#![allow(non_snake_case)]

use crate::nmath::*;
use libm::expm1;

pub fn r_d__0(log_p: bool) -> f64 {
    if log_p {
        f64::NEG_INFINITY
    } else {
        0.0
    }
}

pub fn r_d__1(log_p: bool) -> f64 {
    if log_p {
        0.0
    } else {
        1.0
    }
}

pub fn r_dt_0(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        r_d__0(log_p)
    } else {
        r_d__1(log_p)
    }
}

pub fn r_dt_1(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        r_d__1(log_p)
    } else {
        r_d__0(log_p)
    }
}

pub fn r_d_lval(p: f64, lower_tail: bool) -> f64 {
    if lower_tail {
        p
    } else {
        // Using 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy.
        0.5 - p + 0.5
    }
}

pub fn r_d_cval(p: f64, lower_tail: bool) -> f64 {
    if lower_tail {
        // Using 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy.
        0.5 - p + 0.5
    } else {
        p
    }
}

pub fn r_dt_qiv(p: f64, lower_tail: bool, log_p: bool) -> f64 {
    if log_p {
        if lower_tail {
            p.exp()
        } else {
            -expm1(p)
        }
    } else {
        r_d_lval(p, lower_tail)
    }
}

pub fn r_dt_civ(p: f64, lower_tail: bool, log_p: bool) -> f64 {
    if log_p {
        if lower_tail {
            -expm1(p)
        } else {
            p.exp()
        }
    } else {
        r_d_cval(p, !lower_tail)
    }
}

/// Calculate the boundaries exactly for q*() functions.
/// Often left = ML_NEGINF, and very often right = ML_POSINF;
///
/// R_Q_P01_boundaries(p, left, right)  :<==>
///
/// R_Q_P01_check(p);
/// if (p == R_DT_0) return left;
/// if (p == R_DT_1) return right;
///
/// the following implementation should be more efficient (less tests):
///
/// This was originally a macro, but it is now a function.
/// At the caller site, if the return value is not None, then return the
/// result immediately.
pub fn r_q_p01_boundaries(
    p: f64,
    left: f64,
    right: f64,
    lower_tail: bool,
    log_p: bool,
) -> Option<f64> {
    if log_p {
        if p > 0.0 {
            ml_warn_return_nan();
        }
        if p == 0.0 {
            return Some(if lower_tail { right } else { left });
        }
        if p == ML_NEGINF {
            return Some(if lower_tail { left } else { right });
        }
        None
    } else {
        if !(0.0..=1.0).contains(&p) {
            ml_warn_return_nan();
        }
        if p == 0.0 {
            return Some(if lower_tail { left } else { right });
        }
        if p == 1.0 {
            return Some(if lower_tail { right } else { left });
        }
        None
    }
}
