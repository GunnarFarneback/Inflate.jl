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

function getbits(data::AbstractInflateData, n::Int)
    b = 0
    for i = 0:(n-1)
        b |= getbit(data) << i
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
    for i = 1:length(code)
        v = (v << 1) | getbit(data)
        if v < length(code[i])
            return code[i][1 + v]
        end
        v -= length(code[i])
    end
    error("incomplete code table")
end

function get_literal_or_length(data::AbstractInflateData)
    return get_value_from_code(data, data.literal_or_length_code)
end

const base_length = [11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227]
const extra_length_bits = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]

function getlength(data::AbstractInflateData, v::Int)
    if v <= 264
        return v - 254
    elseif v <= 284
        return base_length[v - 264] + getbits(data, extra_length_bits[v - 264])
    else
        return 258
    end
end

function getdist(data::AbstractInflateData)
    b = get_value_from_code(data, data.distance_code)
    if b <= 3
        return b + 1
    else
        extra_bits = fld(b - 2, 2)
        return 1 + ((2 + b % 2) << extra_bits) + getbits(data, extra_bits)
    end
end

function transform_code_lengths_to_code(code_lengths::Vector{Int})
    code = Vector{Int}[]
    for i = 1:length(code_lengths)
        n = code_lengths[i]
        if n > 0
            while n > length(code)
                push!(code, Int[])
            end
            push!(code[n], i - 1)
        end
    end
    return code
end

const order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

function read_code_tables(data::AbstractInflateData)
    hlit = getbits(data, 5) + 257
    hdist = getbits(data, 5) + 1
    hclen = getbits(data, 4) + 4
    code_length_code_lengths = zeros(Int, 19)
    for i = 1:hclen
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
            for i = 1:len
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
                for i = 1:length
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
    for b in x
        adler = update_adler(adler, b)
    end
    return finish_adler(adler)
end

const crc_table = zeros(UInt32, 256)

function make_crc_table()
    for n = 1:256
        c = UInt32(n - 1)
        for k = 1:8
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
    if crc_table[1] == 0
        make_crc_table()
    end
    return 0xffffffff
end

function update_crc(c::UInt32, x::UInt8)
    return crc_table[1 + ((c ⊻ x) & 0xff)] ⊻ (c >> 8)
end

function finish_crc(c)
    return c ⊻ 0xffffffff
end

function crc(x::Vector{UInt8})
    c = init_crc()
    for b in x
        c = update_crc(c, b)
    end
    return finish_crc(c)
end

function read_zero_terminated_data(data::AbstractInflateData)
    s = UInt8[]
    while true
        c = get_aligned_byte(data)
        push!(s, c)
        if c == 0
            break
        end
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
    if CM != 8
        error("unsupported compression method")
    end
    if CINFO > 7
        error("invalid LZ77 window size")
    end
    if FDICT != 0
        error("preset dictionary not supported")
    end
    if mod((UInt(CMF) << 8) | FLG, 31) != 0
        error("header checksum error")
    end
end

function read_gzip_header(data::AbstractInflateData, headers)
    ID1 = get_aligned_byte(data)
    ID2 = get_aligned_byte(data)
    if ID1 != 0x1f || ID2 != 0x8b
        error("not gzipped data")
    end
    CM = get_aligned_byte(data)
    if CM != 8
        error("unsupported compression method")
    end
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

        for i = 1:xlen
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
        if crc16 != (header_crc & 0xffff)
            error("corrupted data, header crc check failed")
        end
    end
end

"""
    inflate(source::Vector{UInt8})

Decompress in memory `source`, in unwrapped deflate format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `InflateStream`.

Reference: [RFC 1951](https://www.ietf.org/rfc/rfc1951.txt)
"""
inflate(source::Vector{UInt8}) = return _inflate(InflateData(source))

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
    for i = [24, 16, 8, 0]
        ADLER |= Int(get_aligned_byte(data)) << i
    end
    if compute_adler_checksum(out) != ADLER
        error("corrupted data")
    end

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
    if crc32 != crc(out)
        error("corrupted data, crc check failed")
    end
    isize = getbits(data, 32)
    if isize != length(out)
        error("corrupted data, length check failed")
    end

    return out
end

"""
    gunzip(filename::AbstractString)

Convenience wrapper for reading a gzip compressed text file. The
result is returned as a string.
"""
function gunzip(filename::AbstractString; kwargs...)
    return String(gunzip(read(filename); kwargs...))
end


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

function StreamingInflateData(stream::IO)
    return StreamingInflateData(stream, 0, 0, fixed_literal_or_length_table,
                                fixed_distance_table,
                                zeros(UInt8, buffer_size), 1, 1,
                                true, 0, -2, false, init_crc())
end

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
    for i = [24, 16, 8, 0]
        stored_adler |= Int(get_aligned_byte(stream.data)) << i
    end
    if computed_adler != stored_adler
        error("corrupted data")
    end
end

function read_gzip_trailer(stream::GunzipStream)
    crc = finish_crc(stream.crc)
    skip_bits_to_byte_boundary(stream.data)
    crc32 = getbits(stream.data, 32)
    if crc32 != crc
        error("corrupted data, crc check failed")
    end
    isize = getbits(stream.data, 32)
    if isize != stream.num_bytes
        error("corrupted data, length check failed")
    end
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

function write_to_buffer(stream::InflateStream, x::UInt8)
    write_to_buffer(stream.data, x)
end

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
        if stream.data.reading_final_block
            return
        end
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

function Base.eof(stream::AbstractInflateStream)
    return stream.data.read_pos == stream.data.write_pos
end

function Base.read(stream::InflateStream, ::Type{UInt8})
    if eof(stream)
        throw(EOFError())
    end

    byte = read_output_byte(stream)

    return byte
end

function Base.read(stream::DecompressStream, ::Type{UInt8})
    if eof(stream)
        throw(EOFError())
    end

    byte = read_output_byte(stream)

    if eof(stream)
        read_decompress_trailer(stream)
    end

    return byte
end

function Base.read(stream::GunzipStream, ::Type{UInt8})
    if eof(stream)
        throw(EOFError())
    end

    byte = read_output_byte(stream)

    if eof(stream)
        read_gzip_trailer(stream)
    end

    return byte
end

end # module
