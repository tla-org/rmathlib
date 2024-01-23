use crate::nmath::ml_warn_return_nan;
use crate::nmath::r_finite;
use crate::rmath::M_PI;

pub fn cospi(mut x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if !r_finite(x) {
        ml_warn_return_nan();
    }

    x = x.abs() % 2.0;
    if x % 1.0 == 0.5 {
        return 0.;
    };
    if x == 1. {
        return -1.;
    };
    if x == 0. {
        return 1.;
    };
    // otherwise
    (M_PI * x).cos()
}
