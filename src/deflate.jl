export DeflateGzipStream, start_incompressible_block, start_huffman_block

mutable struct BitStream
    out::IO
    data::UInt8
    bitpos::Int
end

BitStream(stream::IO) = BitStream(stream, UInt8(0), 0)

function writebits(stream::BitStream, data, num_bits::Int)
    while num_bits > 0
        n = min(num_bits, 8 - stream.bitpos)
        stream.data |= (data << stream.bitpos) & 0xff
        stream.bitpos += n
        if stream.bitpos == 8
            write(stream.out, stream.data)
            stream.bitpos = 0
            stream.data = 0x00
        end
        num_bits -= n
        data >>= n
    end
    return
end

function flushbits(stream::BitStream)
    if stream.bitpos > 0
        write(stream.out, stream.data)
        stream.bitpos = 0
        stream.data = 0x00
    end
    return
end

function write_aligned_byte(stream::BitStream, x::UInt8)
    flushbits(stream)
    write(stream.out, x)
end

function Base.close(stream::BitStream)
    flushbits(stream)
    close(stream.out)
end

# TODO: Add support for length and distance codes.
function write_deflate_codes(out::BitStream, frequencies::Dict,
                             last_block = false)
    # Normalize the frequency keys to Ints (in case they come as Chars).
    int_frequencies = Dict(Int(k) => v for (k, v) in frequencies)
    # Add one occurence of "end of block".
    int_frequencies[256] = 1
    code_lengths = huffman_code_lengths(int_frequencies)
    @assert maximum(values(code_lengths)) <= 15
    code_length_table = zeros(Int, maximum(keys(code_lengths)) + 1)
    for (k, v) in code_lengths
        code_length_table[Int(k) + 1] = v
    end
    # Hardcoded minimal(?) distance code lengths.
    append!(code_length_table, [1, 1])

    # The code lengths are stored in the Deflate block header with
    # another level of compression.
    i = 1
    compressed_code_lengths = []
    while i <= length(code_length_table)
        l = code_length_table[i]
        j = findfirst(!isequal(l), code_length_table[i+1:end])
        if j == nothing
            j = length(code_length_table[i+1:end])+1
        end
        n = j
        if l == 0
            if n < 3
                for k = 1:n
                    push!(compressed_code_lengths, [l])
                end
            elseif n < 11
                push!(compressed_code_lengths, [17, n - 3])
            else
                n = min(n, 138)
                push!(compressed_code_lengths, [18, n - 11])
            end
        else
            if n < 4
                for k = 1:n
                    push!(compressed_code_lengths, [l])
                end
            else
                n = min(n, 7)
                push!(compressed_code_lengths, [l])
                push!(compressed_code_lengths, [16, n - 4])
            end
        end
        i += n
    end

    # Compute frequencies of the compressed description of the code
    # lengths.
    frequency2 = Dict{Int, Int}()
    for x in first.(compressed_code_lengths)
        frequency2[x] = 1 + get(frequency2, x, 0)
    end

    # Huffman code the entries of the compressed description of the
    # code lengths.
    code_length_code_lengths = huffman_code_lengths(frequency2)
    code_length_codes = huffman_codes(code_length_code_lengths)
    # And reorder the code length code lengths, possibly ending up
    # with trailing zeros, which need not be stored.
    #
    # z could be called reordered_code_length_code_lengths, but around
    # here brevity trumps descriptivity.
    z = [get(code_length_code_lengths, k, 0) for k in order]
    hclen = max(4, findlast(z .!= 0))
    hlit = 257
    hdist = 2

    # Now we can start writing the block header.
    writebits(out, hlit - 257, 5)
    writebits(out, hdist - 1, 5)
    writebits(out, hclen - 4, 4)

    for i = 1:hclen
        writebits(out, z[i], 3)
    end

    for c in compressed_code_lengths
        x = first(c)
        writebits(out, code_length_codes[x], code_length_code_lengths[x])
        if x == 16
            writebits(out, c[2], 2)
        elseif x == 17
            writebits(out, c[2], 3)
        elseif x == 18
            writebits(out, c[2], 7)
        end
    end

    return huffman_codes(code_lengths), code_lengths
end

# Note: This can be made substantially more efficient. The sort in the
# loop only needs to move the last element to its right place.
# Alternatively all of the sorting can be replaced by a heap. It is
# also possible to do the code length housekeeping without forming
# vectors of involved codes.
#
# TODO: Revise so it cannot produce code_lengths longer than 15 bits,
# regardless of the frequencies.
function huffman_code_lengths(frequencies::Dict{<:Any, Int})
    entries = [(v, [k]) for (k, v) in frequencies if v > 0]
    code_lengths = Dict((k, 0) for k in keys(frequencies))
    sort!(entries, by = first, lt = >)
    while length(entries) > 1
        n1, x1 = pop!(entries)
        n2, x2 = entries[end]
        x = vcat(x1, x2)
        entries[end] = (n1 + n2, x)
        sort!(entries, by = first, lt = >)
        for t in x
            code_lengths[t] += 1
        end
    end
    return code_lengths
end

# Based on https://discourse.julialang.org/t/covert-bitarray-to-int64/9193/4.
function reversebits(x, n)
    x = (((x & 0xaaaa) >>  1) | ((x & 0x5555) <<  1))
    x = (((x & 0xcccc) >>  2) | ((x & 0x3333) <<  2))
    x = (((x & 0xf0f0) >>  4) | ((x & 0x0f0f) <<  4))
    x = (((x & 0xff00) >>  8) | ((x & 0x00ff) <<  8))
    return x >> (16 - n)
end

function huffman_codes(code_lengths::Dict{<:Any, Int})
    code = 0
    bits = 1
    codes = Dict((k, 0) for k in keys(code_lengths))
    while code < (1 << bits)
        for k in sort([k for k in keys(code_lengths) if code_lengths[k] == bits])
            codes[k] = reversebits(code, bits)
            code += 1
        end
        bits += 1
        code *= 2
    end

    return codes
end

# TODO: Add support for including optional header fields and setting
# MTIME to anything other than 0.
function write_gzip_header(stream::BitStream)
    write_aligned_byte(stream, 0x1f) # ID1
    write_aligned_byte(stream, 0x8b) # ID2
    write_aligned_byte(stream, 0x08) # CM
    write_aligned_byte(stream, 0x00) # FLG
    write_aligned_byte(stream, 0x00) # MTIME
    write_aligned_byte(stream, 0x00) # MTIME
    write_aligned_byte(stream, 0x00) # MTIME
    write_aligned_byte(stream, 0x00) # MTIME
    write_aligned_byte(stream, 0x00) # XFL
    write_aligned_byte(stream, 0xff) # OS, Unknown
end

function write_gzip_trailer(stream::BitStream, crc, length)
    flushbits(stream)
    writebits(stream, crc, 32)
    writebits(stream, length, 32)
end

mutable struct DeflateGzipStream <: IO
    stream::BitStream
    compression_mode::Int
    pending_incompressible_bytes::Vector{UInt8}
    huffman_codes::Dict{Int, Int}
    huffman_code_lengths::Dict{Int, Int}
    crc::UInt32
    length::Int
end

function DeflateGzipStream(out::IO)
    stream = BitStream(out)
    write_gzip_header(stream)
    return DeflateGzipStream(stream, -1, UInt8[], Dict{Int, Int}(),
                             Dict{Int, Int}(), init_crc(), 0)
end

function Base.write(stream::DeflateGzipStream, x::UInt8)
    stream.crc = update_crc(stream.crc, x)
    stream.length += 1
    if stream.compression_mode == 0
        push!(stream.pending_incompressible_bytes, x)
    elseif stream.compression_mode == 1
        error("Not yet implemented.")
    elseif stream.compression_mode == 2
        writebits(stream.stream, stream.huffman_codes[x],
                  stream.huffman_code_lengths[x])
    else
        error("No compression block has been started.")
    end
    return 1
end

function start_incompressible_block(stream::DeflateGzipStream)
    finish_current_block(stream)
    stream.compression_mode = 0
    stream.pending_incompressible_bytes = UInt8[]
    writebits(stream.stream, 0, 1) # Not final block
    writebits(stream.stream, 0, 2) # Compression mode 0
    flushbits(stream.stream)
end

function finish_incompressible_block(stream::DeflateGzipStream)
    @assert stream.compression_mode == 0
    len = length(stream.pending_incompressible_bytes)
    @assert len <= 0xffff
    writebits(stream.stream, len, 16)
    writebits(stream.stream, len ⊻ 0xffff, 16)
    write(stream.stream.out, stream.pending_incompressible_bytes)
end

function start_huffman_block(stream::DeflateGzipStream, frequencies::Dict)
    finish_current_block(stream)
    stream.compression_mode = 2
    writebits(stream.stream, 0, 1) # Not final block
    writebits(stream.stream, 2, 2) # Compression mode 2
    codes, lengths = write_deflate_codes(stream.stream, frequencies)
    stream.huffman_codes = codes
    stream.huffman_code_lengths = lengths
end

function finish_current_block(stream::DeflateGzipStream)
    if stream.compression_mode == 0
        finish_incompressible_block(stream)
    elseif stream.compression_mode == 1
    elseif stream.compression_mode == 2
        writebits(stream.stream, stream.huffman_codes[256],
                  stream.huffman_code_lengths[256])
    end
end

function Base.close(stream::DeflateGzipStream)
    finish_current_block(stream)
    writebits(stream.stream, 1, 1) # Final block
    writebits(stream.stream, 1, 2) # Compression mode 1
    writebits(stream.stream, 0, 7) # End of block
    write_gzip_trailer(stream.stream, finish_crc(stream.crc),
                       stream.length)
    close(stream.stream)
end
