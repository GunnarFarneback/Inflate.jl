mutable struct GzHeader
    text::Cint           # true if compressed data believed to be text */
    time::Culong         # modification time */
    xflags::Cint         # extra flags (not used when writing a gzip file) */
    os::Cint             # operating system */
    extra::Ptr{Cuchar}   # pointer to extra field or Z_NULL if none */
    extra_len::Cuint     # extra field length (valid if extra != Z_NULL) */
    extra_max::Cuint     # space at extra (only when reading header) */
    name::Ptr{Cuchar}    # pointer to zero-terminated file name or Z_NULL */
    name_max::Cuint      # space at name (only when reading header) */
    comment::Ptr{Cuchar} # pointer to zero-terminated comment or Z_NULL */
    comm_max::Cuint      # space at comment (only when reading header) */
    hcrc::Cint           # true if there was or will be a header crc */
    done::Cint           # true when done reading gzip header (not used
end

function GzHeader(mtime, extra, name, comment, include_header_crc)
    return GzHeader(true, mtime, 0, 255,
                    pointer(extra), length(extra), 0,
                    pointer(name), 0,
                    pointer(comment), 0,
                    include_header_crc, 0)
end

function gzip_with_header(text, mtime, extra, name, comment,
                          include_header_crc)
    zstream = CodecZlib.ZStream()
    CodecZlib.deflate_init!(zstream, 9, 15 + 16)
    data = Vector{UInt8}(text)
    zstream.next_in = pointer(data)
    zstream.avail_in = length(data)
    out = zeros(UInt8, 1000)
    zstream.next_out = pointer(out)
    zstream.avail_out = length(out)
    gzheader = GzHeader(mtime, extra, name, comment, include_header_crc)
    ccall((:deflateSetHeader, CodecZlib.libz), Cint,
          (Ref{CodecZlib.ZStream}, Ref{GzHeader}),
          zstream, gzheader)
    CodecZlib.deflate!(zstream, CodecZlib.Z_FINISH)
    CodecZlib.deflate_end!(zstream)
    return out[1:(length(out) - zstream.avail_out)]
end
