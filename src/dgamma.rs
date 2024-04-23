use crate::dpois::dpois_raw;
use crate::dpq::r_d__0;
use crate::nmath::ml_warn_return_nan;
use crate::nmath::r_finite;
use crate::nmath::ML_POSINF;
use libm::log;

///
///  AUTHOR
///    Catherine Loader, catherine@research.bell-labs.com.
///    October 23, 2000.
///
///  Merge in to R:
///  Copyright (C) 2000-2019 The R Core Team
///  Copyright (C) 2004-2019 The R Foundation
///
///  This program is free software; you can redistribute it and/or modify
///  it under the terms of the GNU General Public License as published by
///  the Free Software Foundation; either version 2 of the License, or
///  (at your option) any later version.
///
///  This program is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
///  GNU General Public License for more details.
///
///  You should have received a copy of the GNU General Public License
///  along with this program; if not, a copy is available at
///  <https://www.R-project.org/Licenses/>

/// Computes the density of the Gamma distribution.
///
/// _          1/s (x/s)^{a-1} exp(-x/s)
/// p(x;a,s) = -----------------------
/// _                   (a-1)!
///
/// where 's' is the scale (= 1/lambda in other parametrizations)
/// and 'a' is the shape parameter ( = alpha in other contexts).
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
        /* point mass at 0 */
        return if x == 0.0 {
            ML_POSINF
        } else {
            r_d__0(give_log)
        };
    }
    if x == 0.0 {
        if shape < 1.0 {
            return ML_POSINF;
        }
        if shape > 1.0 {
            return r_d__0(give_log);
        }
        /* else */
        return if give_log { -log(scale) } else { 1.0 / scale };
    }

    let pr: f64;
    if shape < 1.0 {
        pr = dpois_raw(shape, x / scale, give_log);
        return if give_log {
            /* NB: currently *always*  shape/x > 0  if shape < 1:
             * -- overflow to Inf happens, but underflow to 0 does NOT : */
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
    /* else  shape >= 1 */
    pr = dpois_raw(shape - 1.0, x / scale, give_log);
    if give_log {
        pr - log(scale)
    } else {
        pr / scale
    }
}
