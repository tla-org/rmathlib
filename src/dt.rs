use libm::{exp, log1p};

use crate::{
    bd0, dnorm,
    dpq::r_d__0,
    libc::DBL_EPSILON,
    nmath::{ml_warn_return_nan, r_finite, M_LN_SQRT_2PI},
    rmath::M_1_SQRT_2PI,
    stirlerr,
};

/// The density of the Student's t-distribution.
///
/// The t density is evaluated as
///     sqrt(n/2) / ((n+1)/2) * Gamma((n+3)/2) / Gamma((n+2)/2).
///         * (1+x^2/n)^(-n/2)
///         / sqrt( 2 pi (1+x^2/n) )
///
/// This form leads to a stable computation for all
/// values of n, including n -> 0 and n -> infinity.
pub fn dt(x: f64, n: f64, give_log: bool) -> f64 {
    if x.is_nan() || n.is_nan() {
        return x + n;
    }

    if n <= 0.0 {
        return ml_warn_return_nan();
    }

    if !r_finite(x) {
        return r_d__0(give_log);
    }

    if !r_finite(n) {
        return dnorm(x, 0.0, 1.0, give_log);
    }

    #[allow(unused_assignments)]
    let mut u = 0.0;
    let t = -bd0(n / 2.0, (n + 1.0) / 2.0) + stirlerr((n + 1.0) / 2.0) - stirlerr(n / 2.0);
    let x2n = x * x / n; // in [0, Inf]
    let mut ax = 0.0;
    #[allow(unused_assignments)]
    let mut l_x2n = 0.0; // := log(sqrt(1 + x2n)) = log(1 + x2n)/2
    let lrg_x2n = x2n > (1.0 / DBL_EPSILON);
    if lrg_x2n {
        // large x^2/n
        ax = x.abs();
        l_x2n = ax.ln() - (n.ln() / 2.0); // = log(x2n)/2 = 1/2 * log(x^2 / n)
        u = //  log(1 + x2n) * n/2 =  n * log(1 + x2n)/2 =
	    n * l_x2n;
    } else if x2n > 0.2 {
        l_x2n = (1.0 + x2n).ln() / 2.0;
        u = n * l_x2n;
    } else {
        l_x2n = log1p(x2n) / 2.0;
        u = -bd0(n / 2.0, (n + x * x) / 2.0) + x * x / 2.0;
    }

    // R_D_fexp(f,x) :=  (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
    // f = 2pi*(1+x2n)
    //  ==> 0.5*log(f) = log(2pi)/2 + log(1+x2n)/2 = log(2pi)/2 + l_x2n
    //	     1/sqrt(f) = 1/sqrt(2pi * (1+ x^2 / n))
    //		       = 1/sqrt(2pi)/(|x|/sqrt(n)*sqrt(1+1/x2n))
    //		       = M_1_SQRT_2PI * sqrt(n)/ (|x|*sqrt(1+1/x2n))
    if give_log {
        return t - u - (M_LN_SQRT_2PI + l_x2n);
    }

    // else :  if(lrg_x2n) : sqrt(1 + 1/x2n) ='= sqrt(1) = 1
    #[allow(non_snake_case)]
    let I_srqr_ = if lrg_x2n { n.sqrt() / ax } else { exp(-l_x2n) };
    exp(t - u) * M_1_SQRT_2PI * I_srqr_
}
