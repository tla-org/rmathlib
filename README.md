# rmathlib

A Rust port of [R's C Library of Special Functions](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).

## Benefits

Some benefits of this port over the native C code are:

- Support `wasm32-unknown-unknown` target (avoiding emscripten).
- Clearer documentation with the help of `cargo doc`.
- Easier to read thanks to `cargo fmt`.
- Enables Go To Definition (code navigation) when loading this package as a dependency.
- More clarity about the used compiler and supported Rust versions.

## Status

The following functions have been ported:

Distribution | Density | Probability | Quantile | Random Generation
--- | :---: | :---: | :---: | :---:
Normal | `dnorm` | `pnorm` | `qnorm` |
Student's t | `dt` | `pt`, `pnt` | |
Beta | | `pbeta` | |
Poisson | `dpois` | | |
Gamma | `dgamma` | `pgamma` |

## License

The original R code is licensed under the GPL-2.0.
Therefore, this port is also licensed under the GPL-2.0.

The GPL-2.0 is known to be a very restrictive license.
Each project that includes GPL-2.0 code should also be licensed under the GPL-2.0 (
One exception to this restriction is to use this library only for a small part of the functionality of your project (see, e.g., [here](https://opensource.stackexchange.com/questions/1579)).
In such cases, your project would be considered not to be a "derivative work" and therefore not subject to the GPL-2.0.

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

