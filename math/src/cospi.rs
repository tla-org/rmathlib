use crate::nmath::*;
use crate::rmath::*;

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
    (M_PI * x).cos()
}

pub fn sinpi(mut x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if !r_finite(x) {
        ml_warn_return_nan();
    }

    x = x % 2.0;
    if x <= -1.0 {
        x += 2.0;
    } else if x > 1.0 {
        x -= 2.0;
    }
    if x == 0.0 || x == 1.0 {
        return 0.0;
    }
    if x == 0.5 {
        return 1.0;
    }
    if x == -0.5 {
        return -1.0;
    }
    (M_PI * x).sin()
}

pub fn tanpi(mut x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if !r_finite(x) {
        ml_warn_return_nan();
    }

    x = x % 1.0;
    if x <= -0.5 {
        x += 1.0;
    } else if x > 0.5 {
        x -= 1.0;
    }

    if x == 0.0 {
        0.0
    } else if x == 0.5 {
        ML_NAN
    } else if x == 0.25 {
        1.0
    } else if x == -0.25 {
        -1.0
    } else {
        (M_PI * x).tan()
    }
}
