use libm::exp;
use libm::log;
use libm::log1p;

use crate::dpq::r_d_cval;
use crate::dpq::r_dt_0;
use crate::lbeta;
use crate::libc::fabs;
use crate::nmath::ml_warn_return_nan;
use crate::nmath::r_finite;
use crate::pbeta;
use crate::pnorm;
use crate::rmath::M_LN2;

/// The PDF of the Student's t-distribution.
pub fn pt(x: f64, n: f64, mut lower_tail: bool, log_p: bool) -> f64 {
    if x.is_nan() || n.is_nan() {
        return x + n;
    }

    if n <= 0.0 {
        return ml_warn_return_nan();
    }

    if !r_finite(x) {
        return r_dt_0(lower_tail, log_p);
    }

    if !r_finite(n) {
        return pnorm(x, 0.0, 1.0, lower_tail, log_p);
    }

    let nx = 1.0 + (x / n) * x;

    let val = if nx > 1e100 {
        /* <==>  x*x > 1e100 * n  */
        /* Danger of underflow. So use Abramowitz & Stegun 26.5.4
           pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
           with z = 1/nx,  a = n/2,  b= 1/2 :
        */
        let lval = -0.5 * n * (2.0 * log(fabs(x)) - log(n)) - lbeta(0.5 * n, 0.5) - log(0.5 * n);
        if log_p {
            lval
        } else {
            exp(lval)
        }
    } else if n > x * x {
        pbeta(x * x / (n + x * x), 0.5 * n, 0.5, false, log_p)
    } else {
        pbeta(x * x / (n + x * x), 0.5 * n, 0.5, true, log_p)
    };
    // Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):
    if x <= 0.0 {
        lower_tail = !lower_tail;
    }

    if log_p {
        if lower_tail {
            log1p(-0.5 * exp(val))
        } else {
            val - M_LN2 // = log(.5* pbeta(....))
        }
    } else {
        r_d_cval(val, lower_tail)
    }
}
