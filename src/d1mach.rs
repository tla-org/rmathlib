use crate::libc::DBL_EPSILON;
use crate::libc::DBL_MAX;
use crate::libc::DBL_MIN;
use crate::rmath::M_LOG10_2;

#[allow(dead_code)]
pub fn d1mach(i: i32) -> f64 {
    match i {
        1 => DBL_MIN,
        2 => DBL_MAX,
        3 => 0.5 * DBL_EPSILON,
        4 => DBL_EPSILON,
        5 => M_LOG10_2,
        _ => 0.0,
    }
}
