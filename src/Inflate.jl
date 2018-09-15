# Pure Julia implementation of decompression of zlib and gzip
# compressed data, as specified by RFC 1950-1952.
#
# Historical note: In memory decompression was implemented in 2013 in
# a gist, https://gist.github.com/GunnarFarneback/8254567. Streaming
# decompression was added in 2018 when it was also turned into a Julia
# package.
module Inflate

export inflate, decompress, gunzip,
       InflateStream, DecompressStream, GunzipStream

# Huffman codes are internally represented by Vector{Vector{Int}},
# where code[k] are a vector of the values with code words of length
# k. Codes are assigned in order from shorter to longer codes and in
# the order listed. E.g.
# [[], [2, 7], [1, 3, 5], [6, 4]]
# would be the code
# 00    - 2
# 01    - 7
# 100   - 1
# 101   - 3
# 110   - 5
# 1110  - 6
# 1111  - 4

const fixed_literal_or_length_table = (Vector{Int})[Int[],
                                                    Int[],
                                                    Int[],
                                                    Int[],
                                                    Int[],
                                                    Int[],
                                                    Int[256:279;],
                                                    Int[0:143; 280:287],
                                                    Int[144:255;]]

const fixed_distance_table = (Vector{Int})[Int[],
                                           Int[],
                                           Int[],
                                           Int[],
                                           Int[0:31;]]

abstract type AbstractInflateData end

mutable struct InflateData <: AbstractInflateData
    bytes::Vector{UInt8}
    current_byte::Int
    bytepos::Int
    bitpos::Int
    literal_or_length_code::Vector{Vector{Int}}
    distance_code::Vector{Vector{Int}}
    crc::UInt32
end

function InflateData(source::Vector{UInt8})
    InflateData(source, 0, 1, 0, fixed_literal_or_length_table,
                fixed_distance_table, init_crc())
end

function get_input_byte(data::InflateData)
    byte = data.bytes[data.bytepos]
    data.bytepos += 1
    data.crc = update_crc(data.crc, byte)
    return byte
end

function getbit(data::AbstractInflateData)
    if data.bitpos == 0
        data.current_byte = Int(get_input_byte(data))
    end
    b = data.current_byte & 1
    data.bitpos += 1
    if data.bitpos == 8
        data.bitpos = 0
    else
        data.current_byte >>= 1
    end
    return b
end

# gets up to N bits or to the next byte boundary. Returns
# the value of the bits as a UInt32 plus the number of bits read.
function get_up_to_n_bits_to_byte_boundary(data::AbstractInflateData, n::Int)
    i = 0
    b = 0
    bitstoread = min((8-data.bitpos) % 8, n)
    while i < bitstoread
        bit = data.current_byte & 1
        data.current_byte >>= 1
        b |= bit << i
        # println("read $bit")
        i += 1
    end
    data.bitpos += bitstoread
    data.bitpos %= 8
    return b, i
end

function getbits(data::AbstractInflateData, n::Int)
    b = 0
    bitpos = 0
    if n < 8
        @inbounds for i = 0 : (n-1)
            b |= getbit(data) << i
        end
        return b
    end

    # step 1: read to the byte boundary
    bs, i = get_up_to_n_bits_to_byte_boundary(data, n)
    b = bs
    n -= i
    bitpos += i

    # step 2: read bytes until n < 8    
    while n >= 8
        bs = Int(get_input_byte(data))
        b |= bs << bitpos
        n -= 8
        bitpos += 8
    end

    # step 3: read remaining bits
    @inbounds for i = 0 : n - 1
        bit = getbit(data)
        b |= bit << bitpos
        bitpos += 1
    end
    return b
end

function skip_bits_to_byte_boundary(data::AbstractInflateData)
    data.bitpos = 0
    return
end

# It is the responsibility of the caller to make sure that bitpos is
# at zero, e.g. by calling skip_bits_to_byte_boundary.
function get_aligned_byte(data::AbstractInflateData)
    return get_input_byte(data)
end

function get_value_from_code(data::AbstractInflateData,
                             code::Vector{Vector{Int}})
    v = 0
    for i = Base.OneTo(length(code))
        v = (v << 1) | getbit(data)
        v < length(code[i]) && return code[i][1 + v]
        v -= length(code[i])
    end
    error("incomplete code table")
end

get_literal_or_length(data::AbstractInflateData) =
    get_value_from_code(data, data.literal_or_length_code)


const base_length = [11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227]
const extra_length_bits = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]

function getlength(data::AbstractInflateData, v::Int)
    v <= 264 && return v - 254
    v <= 284 && return base_length[v - 264] + getbits(data, extra_length_bits[v - 264])
    return 258
end

function getdist(data::AbstractInflateData)
    b = get_value_from_code(data, data.distance_code)
    b <= 3 && return b + 1

    extra_bits = fld(b - 2, 2)
    return 1 + ((2 + b % 2) << extra_bits) + getbits(data, extra_bits)
end

function transform_code_lengths_to_code(code_lengths::Vector{Int})
    code = Vector{Int}[]
    @inbounds for i = Base.OneTo(length(code_lengths))
        n = code_lengths[i]
        n <= 0 && continue

        while n > length(code)
            push!(code, Int[])
        end
        push!(code[n], i - 1)
    end
    return code
end

const order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

function read_code_tables(data::AbstractInflateData)
    hlit = getbits(data, 5) + 257
    hdist = getbits(data, 5) + 1
    hclen = getbits(data, 4) + 4
    code_length_code_lengths = zeros(Int, 19)
    @inbounds for i = Base.OneTo(hclen)
        code_length_code_lengths[1 + order[i]] = getbits(data, 3)
    end
    code_length_code = transform_code_lengths_to_code(code_length_code_lengths)
    code_lengths = zeros(Int, hlit + hdist)
    i = 1
    while i <= hlit + hdist
        c = get_value_from_code(data, code_length_code)
        n = 1
        l = 0
        if c < 16
            l = c
        elseif c == 16
            n = 3 + getbits(data, 2)
            l = code_lengths[i-1]
        elseif c == 17
            n = 3 + getbits(data, 3)
        elseif c == 18
            n = 11 + getbits(data, 7)
        else
            error("invalid code length code")
        end
        code_lengths[i:(i+n-1)] .= l
        i += n
    end
    data.literal_or_length_code = transform_code_lengths_to_code(code_lengths[1:hlit])
    data.distance_code = transform_code_lengths_to_code(code_lengths[(hlit+1):end])
end

function _inflate(data::InflateData)
    out = UInt8[]
    final_block = false
    while !final_block
        final_block = getbits(data, 1) == 1
        compression_mode = getbits(data, 2)
        if compression_mode == 0
            skip_bits_to_byte_boundary(data)
            len = getbits(data, 16)
            nlen = getbits(data, 16)
            if len ⊻ nlen != 0xffff
                error("corrupted data")
            end
            for i = Base.OneTo(len)
                push!(out, get_aligned_byte(data))
            end
            continue
        elseif compression_mode == 1
            data.literal_or_length_code = fixed_literal_or_length_table
            data.distance_code = fixed_distance_table
        elseif compression_mode == 2
            read_code_tables(data)
        else
            error("invalid block compression mode 3")
        end

        while true
            v = get_literal_or_length(data)
            if v < 256
                push!(out, UInt8(v))
            elseif v == 256
                break
            else
                length = getlength(data, v)
                distance = getdist(data)
                for i = Base.OneTo(length)
                    push!(out, out[end - distance + 1])
                end
            end
        end
    end

    return out
end

function init_adler()
    return UInt32(1)
end

function update_adler(adler::UInt32, x::UInt8)
    s1 = Int(adler & 0xffff)
    s2 = Int(adler >> 16)
    s1 += x
    if s1 >= 65521
        s1 -= 65521
    end
    s2 += s1
    if s2 >= 65521
        s2 -= 65521
    end
    return (UInt32(s2) << 16) | UInt32(s1)
end

function finish_adler(adler)
    return adler
end

function compute_adler_checksum(x::Vector{UInt8})
    adler = init_adler()
    @inbounds for b in x
        adler = update_adler(adler, b)
    end
    return finish_adler(adler)
end

const crc_table = zeros(UInt32, 256)

function make_crc_table()
    @inbounds for n = Base.OneTo(256)
        c = UInt32(n - 1)
        @inbounds for k = Base.OneTo(8)
            if (c & 0x00000001) != 0
                c = 0xedb88320 ⊻ (c >> 1)
            else
                c >>= 1
            end
        end
        crc_table[n] = c
    end
end

function init_crc()
    crc_table[1] == 0 && make_crc_table()
    return 0xffffffff
end

update_crc(c::UInt32, x::UInt8) =
    crc_table[1 + ((c ⊻ x) & 0xff)] ⊻ (c >> 8)


finish_crc(c) = c ⊻ 0xffffffff


function crc(x::Vector{UInt8})
    c = init_crc()
    @inbounds for b in x
        c = update_crc(c, b)
    end
    return finish_crc(c)
end

function read_zero_terminated_data(data::AbstractInflateData)
    s = UInt8[]
    while true
        c = get_aligned_byte(data)
        push!(s, c)
        c == 0 && break
    end
    return s
end

function read_decompress_header(data::AbstractInflateData)
    CMF = get_aligned_byte(data)
    FLG = get_aligned_byte(data)
    CM = CMF & 0x0f
    CINFO = CM >> 4
    FLEVEL = FLG >> 6
    FDICT = (FLG >> 5) & 0x01
    CM != 8 && error("unsupported compression method")
    CINFO > 7 && error("invalid LZ77 window size")
    FDICT != 0 && error("preset dictionary not supported")
    mod((UInt(CMF) << 8) | FLG, 31) != 0 && error("header checksum error")
end

function read_gzip_header(data::AbstractInflateData, headers)
    ID1 = get_aligned_byte(data)
    ID2 = get_aligned_byte(data)
    (ID1 != 0x1f || ID2 != 0x8b) && error("not gzipped data")
    
    CM = get_aligned_byte(data)
    CM != 8 && error("unsupported compression method")
    
    FLG = get_aligned_byte(data)
    MTIME = getbits(data, 32)
    XFL = get_aligned_byte(data)
    OS = get_aligned_byte(data)

    if headers != nothing
        headers["mtime"] = MTIME
        headers["os"] = OS
    end

    if (FLG & 0x04) != 0   # FLG.FEXTRA
        xlen = getbits(data, 16)
        if headers != nothing
            headers["fextra"] = zeros(UInt8, xlen)
        end

        @inbounds for i = 1:xlen
            b = get_aligned_byte(data)
            if headers != nothing
                headers["fextra"][i] = b
            end
        end
    end

    if (FLG & 0x08) != 0   # FLG.FNAME
        name = read_zero_terminated_data(data)
        if headers != nothing
            headers["fname"] = String(name[1:end-1])
        end
    end

    if (FLG & 0x10) != 0   # FLG.FCOMMENT
        comment = read_zero_terminated_data(data)
        if headers != nothing
            headers["fcomment"] = String(comment[1:end-1])
        end
    end

    if (FLG & 0x02) != 0   # FLG.FHCRC
        header_crc = finish_crc(data.crc)
        crc16 = getbits(data, 16)
        crc16 != (header_crc & 0xffff) && error("corrupted data, header crc check failed")
    end
end

"""
    inflate(source::Vector{UInt8})

Decompress in memory `source`, in unwrapped deflate format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `InflateStream`.

Reference: [RFC 1951](https://www.ietf.org/rfc/rfc1951.txt)
"""
inflate(source::Vector{UInt8}) = _inflate(InflateData(source))

"""
    decompress(source::Vector{UInt8})

Decompress in memory `source`, in Zlib compressed format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `DecompressStream`.

Reference: [RFC 1950](https://www.ietf.org/rfc/rfc1950.txt)
"""
function decompress(source::Vector{UInt8})
    data = InflateData(source)
    read_decompress_header(data)

    out = _inflate(data)

    skip_bits_to_byte_boundary(data)
    ADLER = 0
    @inbounds for i = [24, 16, 8, 0]
        ADLER |= Int(get_aligned_byte(data)) << i
    end
    compute_adler_checksum(out) != ADLER && error("corrupted data")
    return out
end

"""
    gunzip(source::Vector{UInt8})

Decompress in memory `source`, in Gzip compressed format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `GunzipStream`.

    gzip_headers = Dict{String, Any}()
    gunzip(source::Vector{UInt8}; headers = gzip_headers)

Also retrieve gzip headers in the provided `Dict`.

Reference: [RFC 1952](https://www.ietf.org/rfc/rfc1952.txt)
"""
function gunzip(source::Vector{UInt8}; headers = nothing)
    data = InflateData(source)
    read_gzip_header(data, headers)
    out = _inflate(data)

    skip_bits_to_byte_boundary(data)
    crc32 = getbits(data, 32)
    crc32 != crc(out) && error("corrupted data, crc check failed")

    isize = getbits(data, 32)
    isize != length(out) && error("corrupted data, length check failed")
    
    return out
end

"""
    gunzip(filename::AbstractString)

Convenience wrapper for reading a gzip compressed text file. The
result is returned as a string.
"""
gunzip(filename::AbstractString; kwargs...) = String(gunzip(read(filename); kwargs...))


### Streaming interface. ###

# Max distance for repeated strings is 32768 (RFC 1951). Add two bytes
# margin.
const buffer_size = 32770

mutable struct StreamingInflateData <: AbstractInflateData
    stream::IO
    current_byte::Int
    bitpos::Int
    literal_or_length_code::Vector{Vector{Int}}
    distance_code::Vector{Vector{Int}}
    output_buffer::Vector{UInt8}
    write_pos::Int
    read_pos::Int
    waiting_for_new_block::Bool
    pending_bytes::Int
    distance::Int
    reading_final_block::Bool
    crc::UInt32
end

StreamingInflateData(stream::IO) = 
    StreamingInflateData(stream, 0, 0, fixed_literal_or_length_table,
                            fixed_distance_table,
                            zeros(UInt8, buffer_size), 1, 1,
                            true, 0, -2, false, init_crc())

function get_input_byte(data::StreamingInflateData)
    byte = read(data.stream, UInt8)
    data.crc = update_crc(data.crc, byte)
    return byte
end

abstract type AbstractInflateStream <: IO end

"""
    InflateStream(stream::IO)

Streaming decompression of unwrapped deflate compressed `stream`. For
an in memory counterpart, see `inflate`.

Reference: [RFC 1951](https://www.ietf.org/rfc/rfc1951.txt)
"""
mutable struct InflateStream <: AbstractInflateStream
    data::StreamingInflateData
end

function InflateStream(stream::IO)
    stream = InflateStream(StreamingInflateData(stream))
    getbyte(stream)
    return stream
end

"""
    DecompressStream(stream::IO)

Streaming decompression of Zlib compressed `stream`. For an in memory
counterpart, see `decompress`.

Reference: [RFC 1950](https://www.ietf.org/rfc/rfc1950.txt)
"""
mutable struct DecompressStream <: AbstractInflateStream
    data::StreamingInflateData
    adler::UInt32
end

function DecompressStream(stream::IO)
    stream = DecompressStream(StreamingInflateData(stream), init_adler())
    read_decompress_header(stream.data)
    getbyte(stream)
    return stream
end

"""
    GunzipStream(stream::IO)

Streaming decompression of Gzip compressed `stream`. For an in memory
counterpart, see `gunzip`.

    gzip_headers = Dict{String, Any}()
    GunzipStream(stream::IO; headers = gzip_headers)

Also retrieve gzip headers in the provided `Dict`. The headers are
available directly after the object is constructed.

Reference: [RFC 1952](https://www.ietf.org/rfc/rfc1952.txt)
"""
mutable struct GunzipStream <: AbstractInflateStream
    data::StreamingInflateData
    crc::UInt32
    num_bytes::Int
end

function GunzipStream(stream::IO; headers = nothing)
    stream = GunzipStream(StreamingInflateData(stream), init_crc(), 0)
    read_gzip_header(stream.data, headers)
    getbyte(stream)
    return stream
end

function read_decompress_trailer(stream::DecompressStream)
    computed_adler = finish_adler(stream.adler)
    skip_bits_to_byte_boundary(stream.data)
    stored_adler = 0
    @inbounds for i = [24, 16, 8, 0]
        stored_adler |= Int(get_aligned_byte(stream.data)) << i
    end
    computed_adler != stored_adler && error("corrupted data")
end

function read_gzip_trailer(stream::GunzipStream)
    crc = finish_crc(stream.crc)
    skip_bits_to_byte_boundary(stream.data)
    crc32 = getbits(stream.data, 32)
    crc32 != crc && error("corrupted data, crc check failed")

    isize = getbits(stream.data, 32)
    isize != stream.num_bytes && error("corrupted data, length check failed")
end

function read_output_byte(data::StreamingInflateData)
    byte = data.output_buffer[data.read_pos]
    data.read_pos += 1
    if data.read_pos > buffer_size
        data.read_pos -= buffer_size
    end
    return byte
end

function read_output_byte(stream::AbstractInflateStream)
    byte = read_output_byte(stream.data)
    getbyte(stream)
    return byte
end

function write_to_buffer(data::StreamingInflateData, x::UInt8)
    data.output_buffer[data.write_pos] = x
    data.write_pos += 1
    if data.write_pos > buffer_size
        data.write_pos -= buffer_size
    end
end

write_to_buffer(stream::InflateStream, x::UInt8) = write_to_buffer(stream.data, x)


function write_to_buffer(stream::DecompressStream, x::UInt8)
    write_to_buffer(stream.data, x)
    stream.adler = update_adler(stream.adler, x)
end

function write_to_buffer(stream::GunzipStream, x::UInt8)
    write_to_buffer(stream.data, x)
    stream.crc = update_crc(stream.crc, x)
    stream.num_bytes += 1
end

function getbyte(stream::AbstractInflateStream)
    if stream.data.pending_bytes > 0
        stream.data.pending_bytes -= 1
        if stream.data.distance < 0
            write_to_buffer(stream, get_aligned_byte(stream.data))
        else
            pos = stream.data.write_pos - stream.data.distance
            if pos <= 0
                pos += buffer_size
            end
            write_to_buffer(stream, stream.data.output_buffer[pos])
        end
        return
    end

    if stream.data.waiting_for_new_block
        stream.data.reading_final_block && return

        stream.data.reading_final_block = getbits(stream.data, 1) == 1
        compression_mode = getbits(stream.data, 2)
        if compression_mode == 0
            skip_bits_to_byte_boundary(stream.data)
            len = getbits(stream.data, 16)
            nlen = getbits(stream.data, 16)
            if len ⊻ nlen != 0xffff
                error("corrupted data")
            end
            stream.data.distance = -1
            stream.data.pending_bytes = len
            getbyte(stream)
            return
        elseif compression_mode == 1
            stream.data.literal_or_length_code = fixed_literal_or_length_table
            stream.data.distance_code = fixed_distance_table
        elseif compression_mode == 2
            read_code_tables(stream.data)
        else
            error("invalid block compression mode 3")
        end
        stream.data.waiting_for_new_block = false
    end

    v = get_literal_or_length(stream.data)
    if v < 256
        write_to_buffer(stream, UInt8(v))
    elseif v == 256
        stream.data.waiting_for_new_block = true
        getbyte(stream)
    else
        stream.data.pending_bytes = getlength(stream.data, v)
        stream.data.distance = getdist(stream.data)
        getbyte(stream)
    end
end

Base.eof(stream::AbstractInflateStream) = stream.data.read_pos == stream.data.write_pos

function Base.read(stream::InflateStream, ::Type{UInt8})
    eof(stream) && throw(EOFError())
    return read_output_byte(stream)
end

function Base.read(stream::DecompressStream, ::Type{UInt8})
    eof(stream) && throw(EOFError())
    
    byte = read_output_byte(stream)

    eof(stream) && read_decompress_trailer(stream)
    return byte
end

function Base.read(stream::GunzipStream, ::Type{UInt8})
    eof(stream) && throw(EOFError())
    
    byte = read_output_byte(stream)
    
    eof(stream) && read_gzip_trailer(stream)
    return byte
end

end # module
