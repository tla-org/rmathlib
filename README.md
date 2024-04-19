# rmathlib

A Rust port of [R's C Library of Special Functions](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-preview-blue.svg
[docs-latest-url]: https://rmathlib.poweranalyses.org

## Benefits

Some benefits of this port over the native C code are:

- Support `wasm32-unknown-unknown` target (avoiding emscripten).
- Clearer documentation with the help of `cargo doc`.
- Easier to read thanks to `cargo fmt`.
- Enables Go To Definition (code navigation) when loading this package as a dependency.
- More clarity about the used compiler and supported Rust versions.

## Developer Notes

Some tips for debugging the C code in `test/nmath/`:

To print from C, set

```rust
std::env::set_var("CFLAGS", "-DDEBUG_bratio");
```

and use `--nocapture` like so:

```sh
$ cargo watch -x 'test -- --nocapture'
```

When printing inside C, verify that that the numbers are printed correctly.
`REprintf` seems to not always print numbers correctly.
To fix that, `REprintf` can just be replaced with `printf` (and some `\n`'s).

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

