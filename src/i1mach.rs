use core::f32;
use core::f64;
use core::i32;
use core::mem::size_of;

const CHAR_BIT: i32 = 8;
const SIZEOF_CHAR: i32 = 1;
// Using `i32` here because that's the integer type used throughout this code.
const SIZEOF_INT: i32 = size_of::<i32>() as i32;

pub fn i1mach(i: i32) -> i32 {
    match i {
        1 => 5,
        2 => 6,
        3 => 0,
        4 => 0,
        5 => CHAR_BIT * SIZEOF_INT,
        6 => SIZEOF_INT / SIZEOF_CHAR,
        7 => 2,
        8 => CHAR_BIT * SIZEOF_INT - 1,
        9 => i32::MAX,
        10 => f32::RADIX as i32,
        11 => f32::MANTISSA_DIGITS as i32,
        12 => f32::MIN_EXP,
        13 => f32::MAX_EXP,
        14 => f64::MANTISSA_DIGITS as i32,
        15 => f64::MIN_EXP,
        16 => f64::MAX_EXP,
        _ => 0,
    }
}
