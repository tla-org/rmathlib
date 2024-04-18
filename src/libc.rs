pub const DBL_MIN: f64 = f64::MIN;
pub const DBL_MAX: f64 = f64::MAX;
pub const DBL_EPSILON: f64 = f64::EPSILON;

pub fn log(x: f64) -> f64 {
    x.ln()
}

pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

pub fn fabs(x: f64) -> f64 {
    x.abs()
}
