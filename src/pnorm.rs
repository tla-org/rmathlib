use crate::dpq::*;
use crate::nmath::*;
use crate::rmath::*;
use libm::*;

/// Computes the cumulative distribution function of the standard normal distribution
///
/// The main computation evaluates near-minimax approximations derived
/// from those in "Rational Chebyshev approximations for the error
/// function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
/// transportable program uses rational functions that theoretically
/// approximate the normal distribution function to at least 18
/// significant decimal digits.  The accuracy achieved depends on the
/// arithmetic system, the compiler, the intrinsic functions, and
/// proper selection of the machine-dependent constants.
///
/// REFERENCE
///
/// Cody, W. D. (1993).
/// ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
/// Special Function Routines and Test Drivers".
/// ACM Transactions on Mathematical Software. 19, 22-32.
pub fn pnorm(mut x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    let mut p: f64;
    let mut cp: f64 = ML_NAN;

    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        return x + mu + sigma;
    }

    if !r_finite(x) && mu == x {
        return ML_NAN;
    }

    if sigma <= 0.0 {
        if sigma < 0.0 {
            return ml_warn_return_nan();
        }
        return if x < mu {
            r_dt_0(lower_tail, log_p)
        } else {
            r_dt_1(lower_tail, log_p)
        };
    }
    p = (x - mu) / sigma;
    if !r_finite(p) {
        return if x < mu {
            r_dt_0(lower_tail, log_p)
        } else {
            r_dt_1(lower_tail, log_p)
        };
    }
    x = p;

    let i_tail = if lower_tail { 0 } else { 1 };
    pnorm_both(x, &mut p, &mut cp, i_tail, log_p);

    if lower_tail {
        p
    } else {
        cp
    }
}

fn d_2(x: f64) -> f64 {
    ldexp(x, -1)
}

#[allow(clippy::too_many_arguments)]
fn do_del(
    xx: f64,
    x: f64,
    cum: &mut f64,
    ccum: &mut f64,
    log_p: bool,
    lower: bool,
    upper: bool,
    temp: f64,
) {
    let xsq = ldexp(ldexp(xx, 4).trunc(), -4);
    let del = (xx - xsq) * (xx + xsq);
    if log_p {
        *cum = (-xsq * d_2(xsq)) - d_2(del) + temp.ln();
        if (lower && x > 0.0) || (upper && x <= 0.0) {
            *ccum = log1p(-exp(-xsq * d_2(xsq)) * exp(-d_2(del)) * temp);
        }
    } else {
        *cum = (-xsq * d_2(xsq)).exp() * (-d_2(del)).exp() * temp;
        *ccum = 1.0 - *cum;
    }
}

fn swap_tail(x: f64, cum: &mut f64, ccum: &mut f64, lower: bool) {
    if x > 0.0 {
        let temp = *cum;
        if lower {
            *cum = *ccum;
        } else {
            *ccum = temp;
        }
    }
}

fn pnorm_both(x: f64, cum: &mut f64, ccum: &mut f64, i_tail: i32, log_p: bool) {
    let a: [f64; 5] = [
        2.235_252_035_460_683_7,
        161.028_231_068_555_87,
        1_067.689_485_460_370_9,
        18_154.981_253_343_56,
        0.065_682_337_918_207_45,
    ];
    let b: [f64; 4] = [
        47.202_581_904_688_245,
        976.098_551_737_776_7,
        10_260.932_208_618_979,
        45_507.789_335_026_73,
    ];
    let c: [f64; 9] = [
        0.398_941_512_088_134_66,
        8.883_149_794_388_377,
        93.506_656_132_177_85,
        597.270_276_394_800_2,
        2_494.537_585_290_372_6,
        6_848.190_450_536_283,
        11_602.651_437_647_35,
        9_842.714_838_383_978,
        1.076_557_677_372_019_2e-8,
    ];
    let d: [f64; 8] = [
        22.266_688_044_328_117,
        235.387_901_782_625,
        1_519.377_599_407_554_7,
        6_485.558_298_266_761,
        18_615.571_640_885_097,
        34_900.952_721_145_98,
        38_912.003_286_093_27,
        19_685.429_676_859_992,
    ];
    let p: [f64; 6] = [
        0.215_898_534_057_957,
        0.127_401_161_160_247_36,
        0.022235277870649807,
        0.001_421_619_193_227_893_4,
        2.9112874951168792e-5,
        0.023_073_441_764_940_174,
    ];
    let q: [f64; 5] = [
        1.284_260_096_144_911,
        0.468_238_212_480_865_1,
        0.065_988_137_868_928_56,
        0.003_782_396_332_027_582_4,
        7.297_515_550_839_662e-5,
    ];

    let mut xden: f64;
    let mut xnum: f64;
    let mut temp: f64;
    let mut _del: f64;

    let xsq: f64;

    let _min = ML_DBL_MAX;

    let mut _i: i32;

    if x.is_nan() {
        *ccum = x;
        *cum = *ccum;
        return;
    }

    let eps: f64 = ML_DBL_EPSILON * 0.5;

    let lower: bool = i_tail != 1;
    let upper: bool = i_tail != 0;

    let y: f64 = x.abs();
    if y <= 0.67448975 {
        if y > eps {
            xsq = x * x;
            xnum = a[4] * xsq;
            xden = xsq;
            for i in 0..3 {
                xnum = (xnum + a[i]) * xsq;
                xden = (xden + b[i]) * xsq;
            }
        } else {
            xnum = 0.0;
            xden = 0.0;
        }

        temp = x * (xnum + a[3]) / (xden + b[3]);
        if lower {
            *cum = 0.5 + temp;
        }
        if upper {
            *ccum = 0.5 - temp;
        }
        if log_p {
            if lower {
                *cum = (*cum).ln();
            }
            if upper {
                *ccum = (*ccum).ln();
            }
        }
    } else if y < M_SQRT_32 {
        /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

        xnum = c[8] * y;
        xden = y;
        for i in 0..7 {
            xnum = (xnum + c[i]) * y;
            xden = (xden + d[i]) * y;
        }
        temp = (xnum + c[7]) / (xden + d[7]);

        do_del(y, x, cum, ccum, log_p, lower, upper, temp);
        swap_tail(x, cum, ccum, lower);
    } else if log_p && y < 1e170
        || lower && -37.5193 < x && x < 8.2924
        || upper && -8.2924 < x && x < 37.5193
    {
        /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
        xsq = 1.0 / (x * x);
        xnum = p[5] * xsq;
        xden = xsq;
        for i in 0..4 {
            xnum = (xnum + p[i]) * xsq;
            xden = (xden + q[i]) * xsq;
        }
        temp = xsq * (xnum + p[4]) / (xden + q[4]);
        temp = (M_1_SQRT_2PI - temp) / y;
        do_del(x, x, cum, ccum, log_p, lower, upper, temp);
        swap_tail(x, cum, ccum, lower);
    } else {
        // large x such that probs are 0 or 1
        if x > 0.0 {
            *cum = r_d__1(log_p);
            *ccum = r_d__0(log_p);
        } else {
            *cum = r_d__0(log_p);
            *ccum = r_d__1(log_p);
        }

        // R normalizes by default so NO_DENORMS is skipped.
        return;
    }
}
