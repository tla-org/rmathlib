use std::f64::{MANTISSA_DIGITS, MAX, MIN_EXP};

use libm::ldexp;

use crate::nmath::*;
use crate::rmath::*;

/// Compute the density of the normal distribution.
pub fn dnorm4(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64 {
    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        return x + mu + sigma;
    }
    if sigma < 0.0 {
        return ML_NAN; // Placeholder for warning
    }
    if !sigma.is_finite() {
        return 0.0;
    }
    if !x.is_finite() && mu == x {
        return ML_NAN;
    }
    if sigma == 0.0 {
        return if x == mu { ML_POSINF } else { 0.0 };
    }

    let x = (x - mu) / sigma;

    if !x.is_finite() {
        return 0.0;
    }

    let x = x.abs();
    if x >= 2.0 * MAX.sqrt() {
        return 0.0;
    }
    if give_log {
        return -(M_LN_SQRT_2PI + 0.5 * x * x + sigma.ln());
    }

    // Following block is for the else case of MATHLIB_FAST_dnorm
    if x < 5.0 {
        M_1_SQRT_2PI * (-0.5 * x * x).exp() / sigma
    } else {
        let bound = (-2.0 * ML_LN2 * (MIN_EXP as f64 + 1.0 - MANTISSA_DIGITS as f64)).sqrt();
        if x > bound {
            0.0
        } else {
            let x1 = ldexp(r_forceint(ldexp(x, 16)), -16);
            let x2 = x - x1;
            M_1_SQRT_2PI / sigma * (-0.5 * x1 * x1).exp() * ((-0.5 * x2 - x1) * x2).exp()
        }
    }
}
