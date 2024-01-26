use libm::log;
use libm::log1p;

use crate::gammafn;
use crate::lgammacor;
use crate::lgammafn;
use crate::nmath::ML_NAN;
use crate::nmath::ML_NEGINF;
use crate::nmath::ML_POSINF;
use crate::nmath::M_LN_SQRT_2PI;

/// This function returns the value of the log beta function
///
/// log B(a,b) = log G(a) + log G(b) - log G(a+b)
///
/// ## NOTES
///
/// This routine is a translation into C of a Fortran subroutine
/// by W. Fullerton of Los Alamos Scientific Laboratory.
pub fn lbeta(a: f64, b: f64) -> f64 {
    let (mut p, mut q) = (a, a);
    if b < p {
        p = b;
    }
    if b > q {
        q = b;
    }

    // both arguments must be >= 0
    if p < 0.0 {
        return ML_NAN;
    } else if p == 0.0 {
        return ML_POSINF;
    } else if !q.is_finite() {
        return ML_NEGINF;
    }

    if p >= 10.0 {
        // p and q are big.
        let corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
        log(q) * -0.5
            + M_LN_SQRT_2PI
            + corr
            + (p - 0.5) * log(p / (p + q))
            + q * log1p(-p / (p + q))
    } else if q >= 10.0 {
        // p is small, but q is big.
        let corr = lgammacor(q) - lgammacor(p + q);
        lgammafn(p) + corr + p - p * log(p + q) + (q - 0.5) * log1p(-p / (p + q))
    } else {
        // p and q are small: p <= q < 10.
        // R change for very small args
        if p < 1e-306 {
            lgammafn(p) + (lgammafn(q) - lgammafn(p + q))
        } else {
            log(gammafn(p) * (gammafn(q) / gammafn(p + q)))
        }
    }
}
