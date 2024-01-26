use std::f64::INFINITY;
use std::f64::NAN;
use std::f64::NEG_INFINITY;

use crate::dpq::r_d__0;
use crate::rmath::ML_LN2;

pub const ML_POSINF: f64 = INFINITY;
pub const ML_NEGINF: f64 = NEG_INFINITY;

pub fn ml_warn_return_nan() -> f64 {
    println!("argument out of domain");
    ML_NAN
}

pub fn r_finite(x: f64) -> bool {
    x.is_finite()
}

pub const ML_NAN: f64 = NAN;

/// log(sqrt(2*pi)) == log(2*pi)/2
pub const M_LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;

/// log(1 - exp(x)) in more stable form than log1p(- r_d_qiv(x)) :
pub fn r_log1_exp(x: f64) -> f64 {
    if x > -ML_LN2 {
        (-x.exp_m1()).ln()
    } else {
        (-x.exp()).ln_1p()
    }
}

pub fn r_nonint(x: f64) -> bool {
    let nearest_int = x.round();
    (x - nearest_int).abs() > 1e-7 * f64::max(1.0, x.abs())
}

pub fn r_forceint(x: f64) -> f64 {
    x.round()
}

pub fn r_d_nonint_check(x: f64, give_log: bool) -> f64 {
    if r_nonint(x) {
        println!("non-integer x = {x}");
    }
    r_d__0(give_log)
}
