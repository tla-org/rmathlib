pub const ML_POSINF: f64 = f64::INFINITY;
pub const ML_NEGINF: f64 = f64::NEG_INFINITY;

pub fn ml_warn_return_nan() -> f64 {
    println!("argument out of domain");
    f64::NAN
}

pub fn r_finite(x: f64) -> bool {
    x.is_finite()
}

pub const ML_NAN: f64 = f64::NAN;
