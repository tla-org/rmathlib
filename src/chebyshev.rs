//! Chebyshev functions.
//!
//! ## DESCRIPTION
//!
//! - `chebyshev_init` determines the number of terms for the
//!   double precision orthogonal series `dos` needed to insure
//!   the error is no larger than `eta`.  Ordinarily eta will be
//!   chosen to be one-tenth machine precision.
//! - `chebyshev_eval` evaluates the n-term Chebyshev series
//!   `a` at `x`.
//!
//! ## NOTES
//!
//! These routines are translations into Rust by R. Huijzer and J. Storopoli
//! from C translations of Fortran routines by W. Fullerton of Los Alamos
//! Scientific Laboratory. Based on the Fortran routine dcsevl by W. Fullerton.
//! Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).

use crate::nmath::*;

/// `chebyshev_init` determines the number of terms for the
/// double precision orthogonal series `dos` needed to insure
/// the error is no larger than `eta`.  Ordinarily eta will be
/// chosen to be one-tenth machine precision.
pub fn chebyshev_init(dos: &[f64], nos: i32, eta: f64) -> i32 {
    if nos < 1 {
        return 0;
    }

    let mut err = 0.0;
    for (i, &d) in dos.iter().enumerate().rev() {
        err += d.abs();
        if err > eta {
            return i as i32;
        }
    }
    0
}

/// `chebyshev_eval` evaluates the n-term Chebyshev series
/// `a` at `x`.
pub fn chebyshev_eval(x: f64, a: &[f64], n: i32) -> f64 {
    if n < 1 || n > 1000 {
        return ml_warn_return_nan();
    }

    if x < -1.1 || x > 1.1 {
        return ml_warn_return_nan();
    }

    let twox = x * 2.0;
    let mut b0 = 0.0;
    let mut b1 = 0.0;
    let mut b2 = 0.0;

    for i in 1..=n {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[(n - i) as usize];
    }

    (b0 - b2) * 0.5
}
