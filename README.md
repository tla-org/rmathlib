# rmathlib

A Rust port of [R's C Library of Special Functions](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).

## Benefits

Some benefits of this port over the native C code are:

- The Rust code is more safe.
- Support `wasm32-unknown-unknown` target (avoiding emscripten).
- Clearer documentation with the help of `cargo doc`.
- Easier to read thanks to `cargo fmt`.
- Enables go to definition.
- More clarity about the used compiler and supported Rust versions.

## Developer Notes

Some tips for debugging the C code in `test/nmath/`:

If you change the C code, run `cargo clean && cargo test` to force recompilation of the code.

When printing inside C, verify that that the numbers are printed correctly.
`REprintf` seems to not always print numbers correctly.
For example, test the representation via:

```c
double tmp = 1.2;
REprintf("from c: x=%g\n", x);
```

or

```c
#include <stdio.h>

double tmp = 1.2;
printf("from c: x=%g\n", x);
```

