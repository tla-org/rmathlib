use core::mem::size_of;
use core::{f32, f64, i32};
use std::os::raw::c_int;

// Assuming CHAR_BIT is 8. Adjust as needed for your platform.
const CHAR_BIT: i32 = 8;

pub fn rf_i1mach(i: i32) -> i32 {
    match i {
        1 => 5,
        2 => 6,
        3 => 0,
        4 => 0,
        5 => CHAR_BIT * size_of::<c_int>() as i32, // CHAR_BIT * sizeof(int)
        6 => size_of::<c_int>() as i32 / (CHAR_BIT / 8), // sizeof(int)/sizeof(char)
        7 => 2,
        8 => CHAR_BIT * size_of::<c_int>() as i32 - 1, // CHAR_BIT * sizeof(int) - 1
        9 => i32::MAX,                          // INT_MAX
        10 => f32::RADIX as i32,                       // FLT_RADIX
        11 => f32::MANTISSA_DIGITS as i32,             // FLT_MANT_DIG
        12 => f32::MIN_EXP,                     // FLT_MIN_EXP
        13 => f32::MAX_EXP,                     // FLT_MAX_EXP
        14 => f64::MANTISSA_DIGITS as i32,             // DBL_MANT_DIG
        15 => f64::MIN_EXP,                     // DBL_MIN_EXP
        16 => f64::MAX_EXP,                     // DBL_MAX_EXP
        _ => 0,
    }
}
