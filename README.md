# rmathlib

A Rust port of R's C Library of Special Functions

## Benefits

Some benefits of this port over the native C code are:

- Enable statistics support for the `wasm32-unknown-unknown` target (avoiding emscripten).
- The functions are more clearly documented with the help of `cargo doc`.
- The code is easier to read and inspect thanks to `cargo fmt` and follow definition in IDE's.
