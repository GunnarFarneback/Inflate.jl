# Pure Julia implementation of decompression of zlib and gzip
# compressed data, as specified by RFC 1950-1952.
#
# Historical note: In-memory decompression was implemented in 2013 in
# a gist, https://gist.github.com/GunnarFarneback/8254567. Streaming
# decompression was added in 2018 when it was also turned into a Julia
# package.

"""
    Inflate

Pure Julia implementation of decompression of the Deflate format and
the Zlib and Gzip wrapper formats.

In-memory decompression is done by the following functions:

| function | decompresses |
| -------- | ------------ |
| `inflate(data::Vector{UInt8})` | Deflate data |
| `inflate_zlib(data::Vector{UInt8})` | Zlib data |
| `inflate_gzip(data::Vector{UInt8})` | Gzip data |

Streaming decompression is done using the following types:

| stream | decompresses |
| ------ | ------------ |
| `InflateStream(stream::IO)` | Deflate stream |
| `InflateZlibStream(stream::IO)` | Zlib stream |
| `InflateGzipStream(stream::IO)` | Gzip stream |
"""
module Inflate

export inflate, inflate_zlib, inflate_gzip,
       InflateStream, InflateZlibStream, InflateGzipStream

"integer type for literals"
const LInt = Int

# Huffman codes are internally represented by Vector{Vector{LInt}},
# where code[k] are a vector of the values with code words of length
# k. Codes are assigned in order from shorter to longer codes and in
# the order listed. E.g.
# Vector{LInt}[[], [2, 7], [1, 3, 5], [6, 4]]
# would be the code
# 00    - 2
# 01    - 7
# 100   - 1
# 101   - 3
# 110   - 5
# 1110  - 6
# 1111  - 4

const fixed_literal_or_length_table = (Vector{LInt})[[],
                                                     [],
                                                     [],
                                                     [],
                                                     [],
                                                     [],
                                                     [256:279;],
                                                     [0:143; 280:287],
                                                     [144:255;]]

const fixed_distance_table = (Vector{LInt})[[],
                                            [],
                                            [],
                                            [],
                                            [0:31;]]

abstract type AbstractInflateData end

mutable struct InflateData <: AbstractInflateData
    bytes::Vector{UInt8}
    current_byte::UInt8
    bytepos::Int
    bitpos::Int
    literal_or_length_code::Vector{Vector{LInt}}
    distance_code::Vector{Vector{LInt}}
    update_input_crc::Bool
    crc::UInt32
end

function InflateData(source::Vector{UInt8})
    InflateData(source, 0, 1, 0, fixed_literal_or_length_table,
                fixed_distance_table, false, init_crc())
end

"""
    get_input_byte(data::InflateData) -> UInt8
"""
function get_input_byte(data::InflateData)
    byte = data.bytes[data.bytepos]
    data.bytepos += 1
    if data.update_input_crc
        data.crc = update_crc(data.crc, byte)
    end
    return byte
end

# This isn't called when reading Gzip header, so no need to
# consider updating crc.
"""
    get_input_bytes(data::InflateData, n::Int) -> SubArray{UInt8, 1}
"""
function get_input_bytes(data::InflateData, n::Int)
    bytes = @view data.bytes[data.bytepos:(data.bytepos + n - 1)]
    data.bytepos += n
    return bytes
end

"""
    getbit(data::AbstractInflateData) -> Bool
"""
function getbit(data::AbstractInflateData)
    if data.bitpos == 0
        data.current_byte = get_input_byte(data)
    end
    b = data.current_byte % Bool
    data.current_byte >>= 0x1
    data.bitpos = (data.bitpos + 1) & 0x7
    return b
end

"""
    getbits(data::AbstractInflateData, n::Int) -> UInt32
"""
function getbits(data::AbstractInflateData, n::Int)
    b = UInt32(0)
    for i = 0x00:UInt8(n-1)
        b |= UInt32(getbit(data)) << i
    end
    return b
end

"""
    getbitsint(data::AbstractInflateData, n::Int) -> Int
"""
function getbitsint(data::AbstractInflateData, n::Int)
    n < 32 || error("too long bit length `n`")
    return getbits(data, n) % Int
end

function skip_bits_to_byte_boundary(data::AbstractInflateData)
    data.bitpos = 0
    return
end

# It is the responsibility of the caller to make sure that bitpos is
# at zero, e.g. by calling skip_bits_to_byte_boundary.
"""
    get_aligned_byte(data::AbstractInflateData) -> UInt8
"""
function get_aligned_byte(data::AbstractInflateData)
    return get_input_byte(data)
end

"""
    get_value_from_code(data::AbstractInflateData, code) -> LInt
"""
function get_value_from_code(data::AbstractInflateData,
                             code::Vector{Vector{LInt}})
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

"""
    get_literal_or_length(data::AbstractInflateData) -> LInt
"""
function get_literal_or_length(data::AbstractInflateData)
    return get_value_from_code(data, data.literal_or_length_code)
end

const base_length = [11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227]
const extra_length_bits = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]
"""
    getlength(data::AbstractInflateData, v::Int) -> Int
"""
function getlength(data::AbstractInflateData, v::Int)
    if v <= 264
        return v - 254
    elseif v <= 284
        return base_length[v - 264] + getbitsint(data, extra_length_bits[v - 264])
    else
        return 258
    end
end

"""
    getdist(data::AbstractInflateData) -> Int
"""
function getdist(data::AbstractInflateData)
    b = get_value_from_code(data, data.distance_code) % Int
    if b <= 3
        return b + 1
    else
        extra_bits = fld(b - 2, 2)
        return 1 + ((2 + b % 2) << extra_bits) + getbitsint(data, extra_bits)
    end
end

"""
    transform_code_lengths_to_code(code_lengths::Vector{Int}) -> Vector{Vector{LInt}}
"""
function transform_code_lengths_to_code(code_lengths::Vector{Int})
    code = Vector{LInt}[]
    for i = 1:length(code_lengths)
        n = code_lengths[i]
        if n > 0
            while n > length(code)
                push!(code, LInt[])
            end
            push!(code[n], i - 1)
        end
    end
    return code
end

const order = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]

function read_code_tables(data::AbstractInflateData)
    hlit = getbitsint(data, 5) + 257
    hdist = getbitsint(data, 5) + 1
    hclen = getbitsint(data, 4) + 4
    code_length_code_lengths = zeros(Int, 19)
    for i = 1:hclen
        code_length_code_lengths[1 + order[i]] = getbitsint(data, 3)
    end
    code_length_code = transform_code_lengths_to_code(code_length_code_lengths)
    code_lengths = zeros(Int, hlit + hdist)
    i = 1
    while i <= hlit + hdist
        c = get_value_from_code(data, code_length_code) % Int
        n = 1
        l = 0
        if c < 16
            l = c
        elseif c == 16
            n = 3 + getbitsint(data, 2)
            l = code_lengths[i-1]
        elseif c == 17
            n = 3 + getbitsint(data, 3)
        else
            # A code of length 19 can only yield values between 0 and
            # 18 so we can only get here if c == 18.
            n = 11 + getbitsint(data, 7)
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
        final_block = getbit(data)
        compression_mode = getbits(data, 2) % UInt8
        if compression_mode === 0x0
            skip_bits_to_byte_boundary(data)
            len = getbits(data, 16)
            nlen = getbits(data, 16)
            if len ⊻ nlen != 0xffff
                error("corrupted data")
            end
            append!(out, get_input_bytes(data, len % Int))
            continue
        elseif compression_mode === 0x1
            data.literal_or_length_code = fixed_literal_or_length_table
            data.distance_code = fixed_distance_table
        elseif compression_mode === 0x2
            read_code_tables(data)
        else
            error("invalid block compression mode 3")
        end

        while true
            v = get_literal_or_length(data)
            if v < 0x100
                push!(out, v % UInt8)
            elseif v == 0x100
                break
            else
                length = getlength(data, v % Int)
                distance = getdist(data)
                if length <= distance
                    append!(out, @view out[(end - distance + 1):(end - distance + length)])
                else
                    for i = 1:length
                        push!(out, out[end - distance + 1])
                    end
                end
            end
        end
    end

    return out
end

"""
    init_adler() -> Tuple{Int, Int}
"""
function init_adler()
    return (0, 1)
end

"""
    update_adler(adler::Tuple{Int, Int}, x::UInt8) -> Tuple{Int, Int}
"""
function update_adler(adler::Tuple{Int, Int}, x::UInt8)
    s2, s1 = adler
    s1 += x
    if s1 >= 65521
        s1 -= 65521
    end
    s2 += s1
    if s2 >= 65521
        s2 -= 65521
    end
    return (s2, s1)
end


"""
    finish_adler(adler::Tuple{Int, Int}, x::UInt8) -> UInt32
"""
function finish_adler(adler::Tuple{Int, Int})
    s2, s1 = adler
    return (UInt32(s2) << 16) | UInt32(s1)
end

"""
    compute_adler_checksum(x::Vector{UInt8}) -> UInt32
"""
@inline function compute_adler_checksum(x::Vector{UInt8})
    adler = init_adler()
    for b in x
        adler = update_adler(adler, b)
    end
    return finish_adler(adler)
end

"""
    read_stored_adler(data::AbstractInflateData) -> UInt32
"""
function read_stored_adler(data::AbstractInflateData)
    skip_bits_to_byte_boundary(data)
    stored_adler = UInt32(0)
    for i = 1:4
        stored_adler = stored_adler << 0x8 | get_aligned_byte(data)
    end
    return stored_adler
end

const crc_table = zeros(UInt32, 256)

function make_crc_table()
    for n = 1:256
        c = UInt32(n - 1)
        for k = 1:8
            if (c & 0x00000001) != 0
                c = 0xedb88320 ⊻ (c >> 0x1)
            else
                c >>= 0x1
            end
        end
        crc_table[n] = c
    end
end

"""
    init_crc() -> UInt32
"""
function init_crc()
    if crc_table[1] == 0 # FIXME: `crc_table[1]` is always zero.
        make_crc_table() # `crc_table` should be initialized at compile-time.
    end
    return 0xffffffff
end

"""
    update_crc(c::UInt32, x::UInt8) -> UInt32
"""
@inline function update_crc(c::UInt32, x::UInt8)
    @inbounds return crc_table[1 + ((c % UInt8) ⊻ x)] ⊻ (c >> 0x8)
end

"""
    finish_crc(c::UInt32) -> UInt32
"""
function finish_crc(c::UInt32)
    return c ⊻ 0xffffffff
end

"""
    crc(x::Vector{UInt8}) -> UInt32
"""
function crc(x::Vector{UInt8})
    c = init_crc()
    for b in x
        c = update_crc(c, b)
    end
    return finish_crc(c)
end

"""
    read_zero_terminated_data(data::AbstractInflateData) -> Vector{UInt8}
"""
function read_zero_terminated_data(data::AbstractInflateData)
    s = UInt8[]
    while true
        c = get_aligned_byte(data)
        push!(s, c)
        if c === 0x00
            break
        end
    end
    return s
end

function read_zlib_header(data::AbstractInflateData)
    CMF = get_aligned_byte(data)
    FLG = get_aligned_byte(data)
    CM = CMF & 0x0f
    CINFO = CMF >> 0x4
    FLEVEL = FLG >> 0x6
    FDICT = (FLG >> 0x5) & 0x01
    if CM !== 0x8
        error("unsupported compression method")
    end
    if CINFO > 0x7
        error("invalid LZ77 window size")
    end
    if FDICT !== 0x0
        error("preset dictionary not supported")
    end
    if mod((UInt16(CMF) << 0x8) | FLG, 31) != 0
        error("header checksum error")
    end
end

function read_gzip_header(data::AbstractInflateData, headers, compute_crc)
    data.update_input_crc = compute_crc
    ID1 = get_aligned_byte(data)
    ID2 = get_aligned_byte(data)
    if ID1 !== 0x1f || ID2 !== 0x8b
        error("not gzipped data")
    end
    CM = get_aligned_byte(data)
    if CM !== 0x8
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

    if (FLG & 0x04) !== 0x0   # FLG.FEXTRA
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

    if (FLG & 0x08) !== 0x0   # FLG.FNAME
        name = read_zero_terminated_data(data)
        if headers != nothing
            headers["fname"] = String(name[1:end-1])
        end
    end

    if (FLG & 0x10) !== 0x0   # FLG.FCOMMENT
        comment = read_zero_terminated_data(data)
        if headers != nothing
            headers["fcomment"] = String(comment[1:end-1])
        end
    end

    if (FLG & 0xe0) !== 0x0
        error("reserved FLG bit set")
    end

    data.update_input_crc = false
    if (FLG & 0x02) !== 0x0   # FLG.FHCRC
        crc16 = getbits(data, 16) % UInt16
        if compute_crc
            header_crc = finish_crc(data.crc)
            if crc16 !== (header_crc % UInt16)
                error("corrupted data, header crc check failed")
            end
        end
    end
end

"""
    inflate(source::Vector{UInt8})

Decompress in-memory `source`, in unwrapped deflate format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `InflateStream`.

Reference: [RFC 1951](https://www.ietf.org/rfc/rfc1951.txt)
"""
inflate(source::Vector{UInt8}) = _inflate(InflateData(source))

"""
    inflate_zlib(source::Vector{UInt8})

Decompress in-memory `source`, in Zlib compressed format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `InflateZlibStream`.

    inflate_zlib(source::Vector{UInt8}; ignore_checksum = true)

Skip computing Adler checksum for consistency checking.

Reference: [RFC 1950](https://www.ietf.org/rfc/rfc1950.txt)
"""
function inflate_zlib(source::Vector{UInt8}; ignore_checksum = false)
    data = InflateData(source)
    read_zlib_header(data)

    out = _inflate(data)

    stored_adler = read_stored_adler(data)
    if !ignore_checksum && compute_adler_checksum(out) != stored_adler
        error("corrupted data, adler checksum error")
    end

    return out
end

"""
    inflate_gzip(source::Vector{UInt8})

Decompress in-memory `source`, in Gzip compressed format. The
output will also be a `Vector{UInt8}`. For a streaming counterpart,
see `InflateGzipStream`.

    gzip_headers = Dict{String, Any}()
    inflate_gzip(source::Vector{UInt8}; headers = gzip_headers)

Also retrieve gzip headers in the provided `Dict`.

    inflate_gzip(source::Vector{UInt8}; ignore_checksum = true)

Skip computing CRC of the content, as well as the header, for
consistency checking.

Reference: [RFC 1952](https://www.ietf.org/rfc/rfc1952.txt)
"""
function inflate_gzip(source::Vector{UInt8}; headers = nothing,
                      ignore_checksum = false)
    data = InflateData(source)
    read_gzip_header(data, headers, !ignore_checksum)
    out = _inflate(data)

    skip_bits_to_byte_boundary(data)
    crc32 = getbits(data, 32)
    if !ignore_checksum && crc32 != crc(out)
        error("corrupted data, crc check failed")
    end
    isize = getbits(data, 32)
    if isize != length(out) % UInt32
        error("corrupted data, length check failed")
    end

    return out
end

"""
    inflate_gzip(filename::AbstractString) -> String

Convenience wrapper for reading a gzip compressed text file. The
result is returned as a string.
"""
function inflate_gzip(filename::AbstractString; kwargs...)
    return String(inflate_gzip(read(filename); kwargs...))
end


### Streaming interface. ###

# Max distance for repeated strings is 32768 (RFC 1951). Make it
# twice the size for better performance.
const buffer_size = 65536

mutable struct StreamingInflateData <: AbstractInflateData
    stream::IO
    input_buffer::Vector{UInt8}
    input_buffer_pos::Int
    current_byte::UInt8
    bitpos::Int
    literal_or_length_code::Vector{Vector{LInt}}
    distance_code::Vector{Vector{LInt}}
    output_buffer::Vector{UInt8}
    write_pos::Int
    read_pos::Int
    waiting_for_new_block::Bool
    pending_bytes::Int
    distance::Int
    reading_final_block::Bool
    update_input_crc::Bool
    crc::UInt32
end

function StreamingInflateData(stream::IO)
    return StreamingInflateData(stream, UInt8[], 1, 0, 0,
                                fixed_literal_or_length_table,
                                fixed_distance_table,
                                zeros(UInt8, buffer_size), 1, 1,
                                true, 0, -2, false, false, init_crc())
end

"""
    get_input_byte(data::StreamingInflateData) -> UInt8
"""
function get_input_byte(data::StreamingInflateData)
    if data.input_buffer_pos > length(data.input_buffer)
        data.input_buffer = read(data.stream, 65536) # TODO: Remove the magic number
        data.input_buffer_pos = 1
    end
    byte = data.input_buffer[data.input_buffer_pos]
    data.input_buffer_pos += 1
    if data.update_input_crc
        data.crc = update_crc(data.crc, byte)
    end
    return byte
end

# This isn't called when reading Gzip header, so no need to
# consider updating crc.
"""
    get_input_bytes(data::StreamingInflateData, n) -> SubArray{UInt8, 1}
"""
function get_input_bytes(data::StreamingInflateData, n)
    if data.input_buffer_pos > length(data.input_buffer)
        data.input_buffer = read(data.stream, 65536) # TODO: Remove the magic number
        data.input_buffer_pos = 1
    end
    n = min(n, length(data.input_buffer) - data.input_buffer_pos + 1)
    bytes = @view data.input_buffer[data.input_buffer_pos:(data.input_buffer_pos + n - 1)]
    data.input_buffer_pos += n
    return bytes
end

abstract type AbstractInflateStream <: IO end

"""
    InflateStream(stream::IO)

Streaming decompression of unwrapped deflate compressed `stream`. For
an in-memory counterpart, see `inflate`.

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
    InflateZlibStream(stream::IO)

Streaming decompression of Zlib compressed `stream`. For an in-memory
counterpart, see `inflate_zlib`.

    InflateZlibStream(stream::IO; ignore_checksum = true)

Skip computing Adler checksum for consistency checking.

Reference: [RFC 1950](https://www.ietf.org/rfc/rfc1950.txt)
"""
mutable struct InflateZlibStream <: AbstractInflateStream
    data::StreamingInflateData
    adler::Tuple{Int, Int}
    compute_adler::Bool
end

function InflateZlibStream(stream::IO; ignore_checksum = false)
    stream = InflateZlibStream(StreamingInflateData(stream), init_adler(),
                               !ignore_checksum)
    read_zlib_header(stream.data)
    getbyte(stream)
    if eof(stream)
        read_trailer(stream)
    end
    return stream
end

"""
    InflateGzipStream(stream::IO)

Streaming decompression of Gzip compressed `stream`. For an in-memory
counterpart, see `inflate_gzip`.

    gzip_headers = Dict{String, Any}()
    InflateGzipStream(stream::IO; headers = gzip_headers)

Also retrieve gzip headers in the provided `Dict`. The headers are
available directly after the object is constructed.

    InflateGzipStream(stream::IO; ignore_checksum = true)

Skip computing CRC of the content, as well as the header, for
consistency checking.

Reference: [RFC 1952](https://www.ietf.org/rfc/rfc1952.txt)
"""
mutable struct InflateGzipStream <: AbstractInflateStream
    data::StreamingInflateData
    crc::UInt32
    num_bytes::Int
    compute_crc::Bool
end

function InflateGzipStream(stream::IO; headers = nothing,
                           ignore_checksum = false)
    stream = InflateGzipStream(StreamingInflateData(stream), init_crc(), 0,
                               !ignore_checksum)
    read_gzip_header(stream.data, headers, !ignore_checksum)
    getbyte(stream)
    if eof(stream)
        read_trailer(stream)
    end
    return stream
end

read_trailer(stream::InflateStream) = nothing

function read_trailer(stream::InflateZlibStream)
    computed_adler = finish_adler(stream.adler)
    stored_adler = read_stored_adler(stream.data)
    if stream.compute_adler && computed_adler != stored_adler
        error("corrupted data, adler checksum error")
    end
end

function read_trailer(stream::InflateGzipStream)
    crc = finish_crc(stream.crc)
    skip_bits_to_byte_boundary(stream.data)
    crc32 = getbits(stream.data, 32)
    if stream.compute_crc && crc32 !== crc
        error("corrupted data, crc check failed")
    end
    isize = getbits(stream.data, 32)
    if isize !== stream.num_bytes % UInt32
        error("corrupted data, length check failed")
    end
end

"""
    read_output_byte(data::StreamingInflateData) -> UInt8
"""
@inline function read_output_byte(data::StreamingInflateData)
    @inbounds byte = data.output_buffer[data.read_pos]
    data.read_pos += 1
    if data.read_pos > buffer_size
        data.read_pos -= buffer_size
    end
    return byte
end

function read_output_bytes!(data::StreamingInflateData, out, i)
    if data.read_pos < data.write_pos
        m = data.write_pos - data.read_pos
    else
        m = buffer_size - data.read_pos + 1
    end
    if m == 1
        out[i] = read_output_byte(data)
        return 1
    end
    n = length(out) - i + 1
    n = min(m, n)
    copyto!(out, i, data.output_buffer, data.read_pos, n)
    data.read_pos += n
    if data.read_pos > buffer_size
        data.read_pos -= buffer_size
    end
    return n
end

"""
    read_output_byte(stream::AbstractInflateStream) -> UInt8
"""
function read_output_byte(stream::AbstractInflateStream)
    byte = read_output_byte(stream.data)
    getbyte(stream)
    return byte
end

function read_output_bytes!(stream::AbstractInflateStream, out, n)
    i = 1
    resized = false
    while i <= n && !eof(stream)
        if length(out) < i
            resize!(out, min(length(out) * 2, n))
            resized = true
        end
        i += read_output_bytes!(stream.data, out, i)
        getbyte(stream)
    end
    i -= 1
    if resized
        resize!(out, i)
    end
    return i
end

function write_to_buffer(data::StreamingInflateData, x::UInt8)
    data.output_buffer[data.write_pos] = x
    data.write_pos += 1
    if data.write_pos > buffer_size
        data.write_pos -= buffer_size
    end
end

function write_to_buffer(data::StreamingInflateData, x::AbstractVector{UInt8})
    for i in eachindex(x)
        data.output_buffer[data.write_pos + i - 1] = x[i]
    end
    data.write_pos += length(x)
    if data.write_pos > buffer_size
        data.write_pos -= buffer_size
    end
end

function write_to_buffer(stream::InflateStream, x)
    write_to_buffer(stream.data, x)
end

function write_to_buffer(stream::InflateZlibStream, x::UInt8)
    write_to_buffer(stream.data, x)
    if stream.compute_adler
        stream.adler = update_adler(stream.adler, x)
    end
end

function write_to_buffer(stream::InflateZlibStream, x::AbstractVector{UInt8})
    write_to_buffer(stream.data, x)
    if stream.compute_adler
        for y in x
            stream.adler = update_adler(stream.adler, y)
        end
    end
end

function write_to_buffer(stream::InflateGzipStream, x::UInt8)
    write_to_buffer(stream.data, x)
    if stream.compute_crc
        stream.crc = update_crc(stream.crc, x)
    end
    stream.num_bytes += 1
end

function write_to_buffer(stream::InflateGzipStream, x::AbstractVector{UInt8})
    write_to_buffer(stream.data, x)
    if stream.compute_crc
        for y in x
            stream.crc = update_crc(stream.crc, y)
        end
    end
    stream.num_bytes += length(x)
end

"""
    getbyte(stream::AbstractInflateStream)
"""
function getbyte(stream::AbstractInflateStream)
    if stream.data.write_pos != stream.data.read_pos
        return nothing
    end

    if stream.data.pending_bytes > 0
        if stream.data.read_pos <= stream.data.write_pos
            n = buffer_size - stream.data.write_pos + 1
        else
            n = stream.data.read_pos - stream.data.write_pos - 1
        end
        n = min(stream.data.pending_bytes, n)
        if stream.data.distance < 0
            # Incompressible data.
            bytes = get_input_bytes(stream.data, n)
            stream.data.pending_bytes -= length(bytes)
            write_to_buffer(stream, bytes)
        else
            pos = stream.data.write_pos - stream.data.distance
            if pos <= 0
                n = min(n, 1 - pos)
                pos += buffer_size
            end
            stream.data.pending_bytes -= n
            write_to_buffer(stream, @view stream.data.output_buffer[pos:(pos + n - 1)])
        end
        return nothing
    end

    if stream.data.waiting_for_new_block
        if stream.data.reading_final_block
            return nothing
        end
        stream.data.reading_final_block = getbit(stream.data)
        compression_mode = getbits(stream.data, 2) % UInt8
        if compression_mode === 0x0
            skip_bits_to_byte_boundary(stream.data)
            len = getbits(stream.data, 16)
            nlen = getbits(stream.data, 16)
            if len ⊻ nlen != 0xffff
                error("corrupted data")
            end
            stream.data.distance = -1
            stream.data.pending_bytes = len
            getbyte(stream)
            return nothing
        elseif compression_mode === 0x1
            stream.data.literal_or_length_code = fixed_literal_or_length_table
            stream.data.distance_code = fixed_distance_table
        elseif compression_mode === 0x2
            read_code_tables(stream.data)
        else
            error("invalid block compression mode 3")
        end
        stream.data.waiting_for_new_block = false
    end

    v = get_literal_or_length(stream.data)
    if v < 0x100
        write_to_buffer(stream, UInt8(v))
    elseif v == 0x100
        stream.data.waiting_for_new_block = true
        getbyte(stream)
    else
        stream.data.pending_bytes = getlength(stream.data, v % Int)
        stream.data.distance = getdist(stream.data)
        getbyte(stream)
    end
    return nothing
end

function Base.eof(stream::AbstractInflateStream)
    return stream.data.read_pos == stream.data.write_pos
end

function Base.read(stream::AbstractInflateStream, ::Type{UInt8})
    if eof(stream)
        throw(EOFError())
    end

    byte = read_output_byte(stream)

    if eof(stream)
        read_trailer(stream)
    end

    return byte
end

function Base.readbytes!(stream::AbstractInflateStream,
                         b::AbstractVector{UInt8}, nb = length(b))
    if eof(stream)
        return 0
    end

    n = read_output_bytes!(stream, b, nb)

    if eof(stream)
        read_trailer(stream)
    end

    return n
end

end # module
