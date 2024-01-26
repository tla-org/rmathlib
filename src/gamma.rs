use crate::chebyshev::*;
use crate::nmath::*;
use crate::rmath::*;
use crate::sinpi;

/// Chebyshev coefficients for gamma function
const GAMCS: [f64; 42] = [
    8.571_195_590_989_33e-3,
    4.415_381_324_841_007e-3,
    5.685_043_681_599_363e-2,
    -4.219_835_396_418_56e-3,
    1.326_808_181_212_460_3e-3,
    -1.893_024_529_798_880_5e-4,
    3.606_925_327_441_245e-5,
    -6.056_761_904_460_864e-6,
    1.055_829_546_302_283_3e-6,
    -1.811_967_365_542_384e-7,
    3.117_724_964_715_322e-8,
    -5.354_219_639_019_687e-9,
    9.193_275_519_859_589e-10,
    -1.577_941_280_288_339_8e-10,
    2.707_980_622_934_954_4e-11,
    -4.646_818_653_825_73e-12,
    7.973_350_192_007_42e-13,
    -1.368_078_209_830_916e-13,
    2.347_319_486_563_800_7e-14,
    -4.027_432_614_949_067e-15,
    6.910_051_747_372_101e-16,
    -1.185_584_500_221_993e-16,
    2.034_148_542_496_374e-17,
    -3.490_054_341_717_406e-18,
    5.987_993_856_485_306e-19,
    -1.027_378_057_872_228e-19,
    1.762_702_816_060_529_8e-20,
    -3.024_320_653_735_306e-21,
    5.188_914_660_218_398e-22,
    -8.902_770_842_456_576e-23,
    1.527_474_068_493_342_6e-23,
    -2.620_731_256_187_363e-24,
    4.496_464_047_830_539e-25,
    -7.714_712_731_336_878e-26,
    1.323_635_453_126_044e-26,
    -2.270_999_412_942_928_7e-27,
    3.896_418_998_003_991_3e-28,
    -6.685_198_115_125_953e-29,
    1.146_998_663_140_024_4e-29,
    -1.967_938_586_345_134_8e-30,
    3.376_448_816_585_338e-31,
    -5.793_070_335_782_136e-32,
];

/// Machine dependent constants for IEEE double precision
const XMIN: f64 = -170.5674972726612;
const XMAX: f64 = 171.61447887182298;
const XSML: f64 = 2.2474362225598545e-308;
const DXREL: f64 = 1.490_116_119_384_765_6e-8;

/// This function computes the value of the gamma function.
///
/// ## NOTES
///
/// This function is a translation into C of a Fortran subroutine
/// by W. Fullerton of Los Alamos Scientific Laboratory.
/// (e.g. <http://www.netlib.org/slatec/fnlib/gamma.f>)
///
/// The accuracy of this routine compares (very) favourably
/// with those of the Sun Microsystems portable mathematical
/// library.
///
/// MM specialized the case of  n!  for n < 50 - for even better precision
pub fn gammafn(x: f64) -> f64 {
    if x.is_nan() {
        return ML_NAN;
    }

    if x == 0.0 || (x < 0.0 && x == x.round()) {
        return ml_warn_return_nan();
    }

    let y = x.abs();

    if y <= 10.0 {
        let mut n = x as i32;
        if x < 0.0 {
            n -= 1;
        }
        let y = x - n as f64;
        n -= 1;
        let mut value = chebyshev_eval(y * 2.0 - 1.0, &GAMCS, GAMCS.len() as i32) + 0.9375;

        if n == 0 {
            return value;
        }
        if n < 0 {
            if x < -0.5 && ((x - (x - 0.5).round()) / x).abs() < DXREL {
                return ml_warn_return_nan();
            }
            if y < XSML {
                println!("gammafn: Range warning");
                return if x > 0.0 { ML_POSINF } else { ML_NEGINF };
            }
            n = -n;
            for i in 0..n {
                value /= x + i as f64;
            }
            value
        } else {
            for i in 1..=n {
                value *= y + i as f64;
            }
            value
        }
    } else {
        if x > XMAX {
            return ML_POSINF;
        }
        if x < XMIN {
            return 0.0;
        }
        let mut value;
        if y <= 50.0 && y == y.trunc() {
            value = 1.0;
            for i in 2..y as i32 {
                value *= i as f64;
            }
        } else {
            value = ((y - 0.5) * y.ln() - y + ML_LN2 + M_LN_SQRT_2PI).exp();
        }
        if x > 0.0 {
            value
        } else {
            if ((x - (x - 0.5).round()) / x).abs() < DXREL {
                return ml_warn_return_nan();
            }
            let sinpiy = sinpi(y);
            if sinpiy == 0.0 {
                println!("gammafn: Range warning - Negative integer arg - overflow");
                return ML_POSINF;
            }
            -M_PI / (y * sinpiy * value)
        }
    }
}
