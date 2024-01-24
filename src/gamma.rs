use std::f64::{
    consts::{LN_2, PI},
    INFINITY, NAN,
};

use log::warn;

use crate::{chebyshev_eval, nmath::M_LN_SQRT_2PI, sinpi};

/// Chebyshev coefficients for gamma function
const GAMCS: [f64; 42] = [
    0.8571195590989331421920062399942e-2,
    0.4415381324841006757191315771652e-2,
    0.5685043681599363378632664588789e-1,
    -0.4219835396418560501012500186624e-2,
    0.1326808181212460220584006796352e-2,
    -0.1893024529798880432523947023886e-3,
    0.3606925327441245256578082217225e-4,
    -0.6056761904460864218485548290365e-5,
    0.1055829546302283344731823509093e-5,
    -0.1811967365542384048291855891166e-6,
    0.3117724964715322277790254593169e-7,
    -0.5354219639019687140874081024347e-8,
    0.9193275519859588946887786825940e-9,
    -0.1577941280288339761767423273953e-9,
    0.2707980622934954543266540433089e-10,
    -0.4646818653825730144081661058933e-11,
    0.7973350192007419656460767175359e-12,
    -0.1368078209830916025799499172309e-12,
    0.2347319486563800657233471771688e-13,
    -0.4027432614949066932766570534699e-14,
    0.6910051747372100912138336975257e-15,
    -0.1185584500221992907052387126192e-15,
    0.2034148542496373955201026051932e-16,
    -0.3490054341717405849274012949108e-17,
    0.5987993856485305567135051066026e-18,
    -0.1027378057872228074490069778431e-18,
    0.1762702816060529824942759660748e-19,
    -0.3024320653735306260958772112042e-20,
    0.5188914660218397839717833550506e-21,
    -0.8902770842456576692449251601066e-22,
    0.1527474068493342602274596891306e-22,
    -0.2620731256187362900257328332799e-23,
    0.4496464047830538670331046570666e-24,
    -0.7714712731336877911703901525333e-25,
    0.1323635453126044036486572714666e-25,
    -0.2270999412942928816702313813333e-26,
    0.3896418998003991449320816639999e-27,
    -0.6685198115125953327792127999999e-28,
    0.1146998663140024384347613866666e-28,
    -0.1967938586345134677295103999999e-29,
    0.3376448816585338090334890666666e-30,
    -0.5793070335782135784625493333333e-31,
];

/// Machine dependent constants for IEEE double precision
const XMIN: f64 = -170.5674972726612;
const XMAX: f64 = 171.61447887182298;
const XSML: f64 = 2.2474362225598545e-308;
const DXREL: f64 = 1.490116119384765696e-8;

/// This function computes the value of the gamma function.
///
/// ## NOTES
///
/// This function is a translation into C of a Fortran subroutine
/// by W. Fullerton of Los Alamos Scientific Laboratory.
/// (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
///
/// The accuracy of this routine compares (very) favourably
/// with those of the Sun Microsystems portable mathematical
/// library.
///
/// MM specialized the case of  n!  for n < 50 - for even better precision
pub fn gammafn(x: f64) -> f64 {
    if x.is_nan() {
        return NAN;
    }

    if x == 0.0 || (x < 0.0 && x == x.round()) {
        warn!("gammafn: Domain warning");
        return NAN;
    }

    let y = x.abs();

    if y <= 10.0 {
        let mut n = x as i32;
        if x < 0.0 {
            n -= 1;
        }
        let y = x - n as f64;
        n -= 1;
        let mut value = chebyshev_eval(y * 2.0 - 1.0, &GAMCS, GAMCS.len()) + 0.9375;

        if n == 0 {
            return value;
        }
        if n < 0 {
            if x < -0.5 && ((x - (x - 0.5).round()) / x).abs() < DXREL {
                warn!("gammafn: Precision warning");
                return NAN;
            }
            if y < XSML {
                warn!("gammafn: Range warning");
                return if x > 0.0 { INFINITY } else { -INFINITY };
            }
            n = -n;
            for i in 0..n {
                value /= x + i as f64;
            }
            return value;
        } else {
            for i in 1..=n {
                value *= y + i as f64;
            }
            return value;
        }
    } else {
        if x > XMAX {
            return INFINITY;
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
            value = ((y - 0.5) * y.ln() - y + LN_2 + M_LN_SQRT_2PI).exp();
        }
        if x > 0.0 {
            value
        } else {
            if ((x - (x - 0.5).round()) / x).abs() < DXREL {
                warn!("gammafn: Precision warning");
                return NAN;
            }
            let sinpiy = sinpi(y);
            if sinpiy == 0.0 {
                warn!("gammafn: Range warning - Negative integer arg - overflow");
                return INFINITY;
            }
            -PI / (y * sinpiy * value)
        }
    }
}
