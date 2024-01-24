use std::f64::{consts::LN_2, INFINITY, NAN};

use log::warn;

use crate::sinpi;

extern "C" {
    pub fn gammafn(x: f64) -> f64;
    pub fn lgammacor(x: f64) -> f64;
}

/// Machine dependent constants for IEEE double precision
const XMAX: f64 = 2.5327372760800758e+305;
const DXREL: f64 = 1.490116119384765625e-8;

/// The function lgammafn computes log|gamma(x)|.  The function
/// lgammafn_sign in addition assigns the sign of the gamma function
/// to the address in the second argument if this is not NULL.
///
/// ## NOTES
///
/// This routine is a translation into C of a Fortran subroutine
/// by W. Fullerton of Los Alamos Scientific Laboratory.
///
/// The accuracy of this routine compares (very) favourably
/// with those of the Sun Microsystems portable mathematical
/// library.
///
/// ./toms708.c  has  gamln()
pub fn lgammafn(x: f64) -> f64 {
    lgammafn_sign(x, None)
}

/// Compute the log gamma function and its sign.
///
/// This function computes the logarithm of the absolute value of the gamma function of x
/// and also sets the sign of the gamma function in `sgn`.
fn lgammafn_sign(x: f64, sgn: Option<&mut i32>) -> f64 {
    if let Some(sgn) = sgn {
        *sgn = 1;
        if x < 0.0 && ((-x).floor() % 2.0) == 0.0 {
            *sgn = -1;
        }
    }

    if x.is_nan() {
        return NAN;
    }

    if x <= 0.0 && x == x.trunc() {
        // Negative integer argument
        return INFINITY; // +Inf, since lgamma(x) = log|gamma(x)|
    }

    let y = x.abs();

    if y < 1e-306 {
        return -y.ln();
    }
    if y <= 10.0 {
        return unsafe { gammafn(x).abs().ln() };
    }

    // y = |x| > 10
    if y > XMAX {
        return INFINITY;
    }

    if x > 0.0 {
        // Positive x
        if x > 1e17 {
            return x * (x.ln() - 1.0);
        } else if x > 4934720.0 {
            return LN_2 + (x - 0.5) * x.ln() - x;
        } else {
            return LN_2 + (x - 0.5) * x.ln() - x + unsafe { lgammacor(x) };
        }
    } else {
        // x < -10; y = -x
        let sinpiy = sinpi(y).abs();

        if sinpiy == 0.0 {
            // Handle error: Negative integer argument
            return NAN;
        }

        let ans = LN_2 + (x - 0.5) * y.ln() - x - sinpiy.ln() - unsafe { lgammacor(y) };

        // Check for accuracy
        if ((x - x.trunc() - 0.5) * ans / x).abs() < DXREL {
            // Warning: answer less than half precision
            // because the argument is too near a negative integer
            warn!("** should NEVER happen! *** [lgamma.rs: Neg.int, y={y}]");
            return NAN; // Placeholder for warning
        }
        return ans;
    }
}
