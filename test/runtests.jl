using Test
using Random
using Inflate
using CodecZlib: GzipCompressorStream, ZlibCompressorStream,
                 DeflateCompressorStream, CodecZlib

empty_string = ""
short_string = "This is a short string."
medium_string = read(pathof(Inflate), String)
long_string = join(fill(medium_string, 1000), short_string)

@testset "Text strings" begin
    for s in [empty_string, short_string, medium_string, long_string]
        @test String(inflate(read(DeflateCompressorStream(IOBuffer(s))))) == s
        @test String(decompress(read(ZlibCompressorStream(IOBuffer(s))))) == s
        @test String(gunzip(read(GzipCompressorStream(IOBuffer(s))))) == s
        @test read(InflateStream(DeflateCompressorStream(IOBuffer(s))), String) == s
        @test read(DecompressStream(ZlibCompressorStream(IOBuffer(s))), String) == s
        @test read(GunzipStream(GzipCompressorStream(IOBuffer(s))), String) == s
    end
end

@testset "Incompressible data" begin
    Random.seed!(1)
    for n in [0, 1, 10, 100, 1000, 10000, 100000, 1000000]
        data = rand(UInt8, n)
        @test inflate(read(DeflateCompressorStream(IOBuffer(data)))) == data
        @test decompress(read(ZlibCompressorStream(IOBuffer(data)))) == data
        @test gunzip(read(GzipCompressorStream(IOBuffer(data)))) == data
        @test read(InflateStream(DeflateCompressorStream(IOBuffer(data)))) == data
        @test read(DecompressStream(ZlibCompressorStream(IOBuffer(data)))) == data
        @test read(GunzipStream(GzipCompressorStream(IOBuffer(data)))) == data
    end
end

# Test gzip headers, including header CRC.
include("gzip_with_header.jl")
os = 255   # Unknown
mtime = 1537006040
fextra = UInt8[0x00, 0xbe, 0xef, 0x00]
fname = "foo.txt"
fcomment = "example string"
gz = gzip_with_header("foo", mtime, os, fextra, fname, fcomment, true)

@testset "Gzip headers" begin
    headers = Dict{String, Any}()
    @test String(gunzip(gz, headers = headers)) == "foo"
    @test headers["mtime"] == mtime
    @test headers["os"] == os
    @test headers["fextra"] == fextra
    @test headers["fname"] == fname
    @test headers["fcomment"] == fcomment

    headers = Dict{String, Any}()
    @test read(GunzipStream(IOBuffer(gz), headers = headers), String) == "foo"
    @test headers["mtime"] == mtime
    @test headers["os"] == os
    @test headers["fextra"] == fextra
    @test headers["fname"] == fname
    @test headers["fcomment"] == fcomment
end

# Test failure cases, mostly corrupt data.
include("provoke_errors.jl")
