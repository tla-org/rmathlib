//!
//! AUTHOR
//!   Catherine Loader, catherine@research.bell-labs.com.
//!   October 23, 2000.
//!
//!   Merge in to R:
//! Copyright (C) 2000-2021 The R Core Team

use crate::dpq::*;
use crate::libc::*;
use crate::nmath::*;
use crate::rmath::*;
use libm::*;

const M_SQRT_2PI: f64 = 2.50662827463100050241576528481104525301; // sqrt(2 * pi)
const X_LRG: f64 = 2.86111748575702815380240589208115399625e+307; // 2^1023 / pi

fn dpois_raw(x: f64, lambda: f64, give_log: bool) -> f64 {
    if lambda == 0.0 {
        if x == 0.0 {
            return r_d__1(give_log);
        } else {
            return r_d__0(give_log);
        }
    }

    if !r_finite(lambda) {
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
            return r_d__0(give_log);
        }

        // TODO: lgammafn
        // return r_d_exp(-lambda + x*log(lambda) - lgammafn(x+1.0), give_log);
    };

    

    return 0.0;
}

/// Poisson distribution
/// Computes the Poisson probability lb^x exp(-lb) / x!.
/// This does not check that x is an integer, since dgamma() may
/// call this with a fractional x argument. Any necessary argument
/// checks should be done in the calling function.
pub fn dpois() {

}
