use crate::{chebyshev_eval, ml_warn_return_nan};

const ALGMCS: [f64; 15] = [
    1.666_389_480_451_863_4e-1,
    -1.384_948_176_067_564e-5,
    9.810_825_646_924_73e-9,
    -1.809_129_475_572_494e-11,
    6.221_098_041_892_606e-14,
    -3.399_615_005_417_722e-16,
    2.683_181_998_482_699e-18,
    -2.868_042_435_334_643e-20,
    3.962_837_061_046_434_7e-22,
    -6.831_888_753_985_767e-24,
    1.429_227_355_942_498_2e-25,
    -3.547_598_158_101_070_4e-27,
    1.025_680_058_010_471e-28,
    -3.401_102_254_316_749e-30,
    1.276_642_195_630_063e-31,
];

/// Machine dependent constants for IEEE double precision
const NALGM: usize = 5;
const XBIG: f64 = 94906265.62425156;
const XMAX: f64 = 3.745194030963158e306;

/// Compute the log gamma correction factor for x >= 10 so that
///
/// log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
///
/// [lgammacor(x) is called Del(x) in other contexts (e.g. dcdflib)]
///
/// ## NOTES
///
/// This routine is a translation into C of a Fortran subroutine
/// written by W. Fullerton of Los Alamos Scientific Laboratory.
///
/// ## SEE ALSO
///
/// Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
/// is faster and cleaner, but is only defined "fast" for half integers.
pub fn lgammacor(x: f64) -> f64 {
    if x < 10.0 {
        return ml_warn_return_nan();
    } else if x >= XMAX {
        println!("lgammacor: Underflow warning");
        // Allow to underflow
    } else if x < XBIG {
        let tmp = 10.0 / x;
        return chebyshev_eval(tmp * tmp * 2.0 - 1.0, &ALGMCS, NALGM) / x;
    }
    1.0 / (x * 12.0)
}
