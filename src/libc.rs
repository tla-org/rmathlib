pub const DBL_MIN: f64 = core::f64::MIN;
pub const DBL_MAX: f64 = core::f64::MAX;

pub fn log(x: f64) -> f64 {
    x.ln()
}

pub fn exp(x: f64) -> f64 {
    x.exp()
}

pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

pub fn fabs(x: f64) -> f64 {
    x.abs()
}
