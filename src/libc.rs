pub const DBL_MIN: f64 = core::f64::MIN;
pub const DBL_MAX: f64 = core::f64::MAX;
pub const DBL_EPSILON: f64 = core::f64::EPSILON;

pub fn log(x: f64) -> f64 {
    x.ln()
}

pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

pub fn fabs(x: f64) -> f64 {
    x.abs()
}

pub fn exp(x: f64) -> f64 {
    x.exp()
}

pub fn max(a: f64, b: f64) -> f64 {
    a.max(b)
}

pub fn min(a: f64, b: f64) -> f64 {
    a.min(b)
}
