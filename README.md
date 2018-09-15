# Inflate.jl

Inflate provides a pure Julia implementation of
[zlib](https://zlib.net) *de*compression functionality, with both in
memory and streaming interfaces.

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
