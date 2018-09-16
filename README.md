# Inflate.jl

Inflate provides a pure Julia implementation of
[zlib](https://zlib.net) *de*compression functionality, with both in
memory and streaming interfaces. This covers decompression of the
Deflate algorithm and the Zlib and Gzip wrapper formats, as specified
in [RFC 1950](https://www.ietf.org/rfc/rfc1950.txt),
[RFC 1951](https://www.ietf.org/rfc/rfc1951.txt), and
[RFC 1952](https://www.ietf.org/rfc/rfc1952.txt).

The main reasons to choose Inflate over
[CodecZlib](https://github.com/bicycle1885/CodecZlib.jl) are:
* 100% Julia code - great for Julia purists.
* No binary dependencies.
* Actually no dependencies at all.
* Can read gzip headers.

You should choose CodecZlib over Inflate if the points above are not
compelling or one or more of the following applies to you:
* Need to compress, not only decompress.
* Want higher speed.
* Want a full-featured streaming interface.
* Want a battle-proven library.

## Documentation

For now, get inspiration from the docstrings and the tests in
`test/runtests.jl`.
