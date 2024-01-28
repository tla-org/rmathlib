//! Rust version by Rik Huijzer and Jose Storopoli

use libm::log;

use crate::dpq::r_d__0;
use crate::dpq::r_d__1;
use crate::dpq::r_d_exp;
use crate::ebd0;
use crate::lgammafn;
use crate::libc::exp;
use crate::libc::sqrt;
use crate::libc::DBL_MIN;
use crate::nmath::ml_warn_return_nan;
use crate::nmath::r_d_nonint_check;
use crate::nmath::r_finite;
use crate::nmath::r_forceint;
use crate::rmath::M_2PI;
use crate::stirlerr;

const M_SQRT_2PI: f64 = 2.506_628_274_631_000_7; // sqrt(2*pi)
const X_LRG: f64 = 2.861_117_485_757_028_3e307; // = 2^1023 / pi

/// Checks argument validity and calls dpois_raw().
///
/// ## AUTHOR
///
/// Catherine Loader, catherine@research.bell-labs.com.
/// October 23, 2000.
pub fn dpois(x: f64, lambda: f64, give_log: bool) -> f64 {
    if x.is_nan() || lambda.is_nan() {
        return x + lambda;
    }
    if lambda < 0.0 {
        return ml_warn_return_nan();
    }
    let x = r_d_nonint_check(x, give_log);
    if x < 0.0 || !r_finite(x) {
        return r_d__0(give_log);
    }

    let x = r_forceint(x);

    dpois_raw(x, lambda, give_log)
}

/// Computes the Poisson probability lb^x exp(-lb) / x!.
///
/// This does not check that x is an integer, since dgamma() may
/// call this with a fractional x argument. Any necessary argument
/// checks should be done in the calling function.
///
/// ## AUTHOR
///
/// Catherine Loader, catherine@research.bell-labs.com.
/// October 23, 2000.
pub fn dpois_raw(x: f64, lambda: f64, give_log: bool) -> f64 {
    if lambda == 0.0 {
        if x == 0.0 {
            return r_d__1(give_log);
        } else {
            return r_d__0(give_log);
        }
    }
    if !r_finite(x) {
        // including for the case where x = lambda = +Inf
        return r_d__0(give_log);
    }
    if x < 0.0 {
        return r_d__0(give_log);
    }
    if x <= lambda * DBL_MIN {
        return r_d_exp(-lambda, give_log);
    }
    if lambda < x * DBL_MIN {
        if !r_finite(x) {
            // lambda < x = +Inf
            return r_d__0(give_log);
        } else {
            return r_d_exp(-lambda + x * log(lambda) - lgammafn(x + 1.0), give_log);
        }
    }
    let (yh, mut yl) = ebd0(x, lambda);
    yl += stirlerr(x);
    let lrg_x = x >= X_LRG; //really large x  <==>  2*pi*x  overflows
    let r = match lrg_x {
        true => M_SQRT_2PI * sqrt(x),
        false => M_2PI * x,
    };
    match give_log {
        true => {
            if lrg_x {
                -yl - yh - log(r)
            } else {
                -yl - yh - (0.5 * log(r))
            }
        }
        false => {
            if lrg_x {
                exp(-yl) / r
            } else {
                exp(-yl) / sqrt(r)
            }
        }
    }
}
