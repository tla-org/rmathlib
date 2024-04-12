use crate::dnorm::dnorm4;
use crate::pnorm::pnorm5;
use crate::qnorm::qnorm5;

use std::f64::consts::LN_2;
use std::f64::consts::PI;
use std::f64::consts::SQRT_2;

pub const M_PI: f64 = PI;
pub const M_SQRT2: f64 = SQRT_2;
pub const M_SQRT_32: f64 = 5.656_854_249_492_381;
pub const M_1_SQRT_2PI: f64 = 0.398_942_280_401_432_7; // 1/sqrt(2pi)
pub const M_LN2: f64 = LN_2;
pub const M_LOG10_2: f64 = 0.301_029_995_663_981_2; // log10(2)

pub fn dnorm(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64 {
    dnorm4(x, mu, sigma, give_log)
}

pub fn pnorm(x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    pnorm5(x, mu, sigma, lower_tail, log_p)
}

pub fn qnorm(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    qnorm5(p, mu, sigma, lower_tail, log_p)
}
