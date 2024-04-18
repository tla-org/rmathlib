use crate::dpq::r_dt_0;
use crate::dpq::r_dt_1;
use crate::nmath::ml_warn_return_nan;
use crate::rmath::M_LN2;

/// Returns distribution function of the beta distribution.
/// ( = The incomplete beta ratio I_x(p,q) ).
fn pbeta_raw(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    if x >= 1 {
        return r_dt_1(lower_tail, log_p);
    }
    if a == 0.0 || b == 0.0 || !a.is_finite() || !b.is_finite() {
        // NB: 0 < x < 1 :
        if a == 0.0 && b == 0.0 {
            // point mass 1/2 at each of {0, 1} :
            return if log_p { -M_LN2 } else { 0.5 };
        }
        if a == 0.0 || a / b == 0.0 {
            // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
            return r_dt_1(lower_tail, log_p);
        }
        if b == 0.0 || b / a == 0.0 {
            // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
            return r_dt_0(lower_tail, log_p);
        }
        // else, remaining case:  a = b = Inf : point mass 1 at 1/2
        if x < 0.5 {
            return r_dt_0(lower_tail, log_p);
        } else {
            return r_dt_1(lower_tail, log_p);
        }
    }
    if x <= 0.0 {
        return r_dt_0(lower_tail, log_p);
    }

    let x1 = 0.5 - x + 0.5;
    let mut w = 0.0;
    let mut wc = 0.0;
    let mut ierr = 0;
    bratio(a, b, x, x1, &mut w, &mut wc, &mut ierr, log_p);

    if ierr != 0 && ierr != 11 && ierr != 14 {
        println!("pbeta_raw({x}, a={a}, b={b}, lower_tail={lower_tail}, log_p={log_p}) -> bratio() gave error code{ierr}");
    }
    if lower_tail {
        w
    } else {
        wc
    }
}

/// Returns distribution function of the beta distribution.
///( = The incomplete beta ratio I_x(p,q) ).
///
/// ## NOTES
///
/// - A wrapper for TOMS708
/// - 'log_p' partially improved over log(p..)
pub fn pbeta(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    if x.is_nan() || a.is_nan() || b.is_nan() {
        return x + a + b;
    }
    if a < 0.0 || b < 0.0 {
        return ml_warn_return_nan();
    }
    pbeta_raw(x, a, b, lower_tail, log_p)
}
