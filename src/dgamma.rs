use crate::dpq::r_d__0;
use crate::libc::log;
use crate::nmath::ml_warn_return_nan;
use crate::nmath::r_finite;
use crate::nmath::ML_POSINF;

/// Computes the density of the gamma distribution,
///
///                 1/s (x/s)^{a-1} exp(-x/s)
///      p(x;a,s) = -----------------------
///                          (a-1)!
///
/// where 's' is the scale (= 1/lambda in other parametrizations)
/// and 'a' is the shape parameter ( = alpha in other contexts).
///
/// ## AUTHOR
///
/// - Original C Code:
///   Catherine Loader, catherine@research.bell-labs.com.
///   October 23, 2000.
///
/// - Rust Implementation:
///   Rik Huijzer and Jose Storopoli.
///   January 29, 2024.
pub fn dgamma(x: f64, shape: f64, scale: f64, give_log: bool) -> f64 {
    if x.is_nan() || shape.is_nan() || scale.is_nan() {
        return x + shape + scale;
    }

    if shape < 0.0 || scale <= 0.0 {
        return ml_warn_return_nan();
    }
    if x < 0.0 {
        return r_d__0(give_log);
    }
    if shape == 0.0 {
        // point mass at 0
        return if x == 0.0 {
            ML_POSINF
        } else {
            return r_d__0(give_log);
        };
    }
    if x == 0.0 {
        if shape < 1.0 {
            return ML_POSINF;
        }
        if shape > 1.0 {
            return r_d__0(give_log);
        }
        //else
        return if give_log { -log(scale) } else { 1.0 / scale };
    }
    if shape < 1.0 {
        let pr = dpois_raw(shape, x / scale, give_log);
        return if give_log {
            /* NB: currently *always* shape/x > 0 if shape < 1:
             * -- overflow to Inf happens, but overflow to 0 does NOT : */
            pr + if r_finite(shape / x) {
                log(shape / x)
            } else {
                /* shape/x overflows to +Inf */
                log(shape) - log(x)
            }
        } else {
            pr * shape / x
        };
    }
    let pr = dpois_raw(shape - 1.0, x / scale, give_log);
    return if give_log {
        pr - log(scale)
    } else {
        pr / scale
    };
}
