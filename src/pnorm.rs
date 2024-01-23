use crate::dpq::*;
use crate::rmath::*;
use crate::nmath::*;
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
    let mut cp: f64 = f64::NAN;

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

    if lower_tail { p } else { cp }
}

fn d_2(x: f64) -> f64 {
    ldexp(x, -1)
}

fn do_del(xx: f64, x: f64, cum: &mut f64, ccum: &mut f64, log_p: bool, lower: bool, upper: bool, temp: f64) {
    let xsq = ldexp(ldexp(xx, 4).trunc(), -4);
    let del = (xx - xsq) * (xx + xsq);
    if log_p {
        *cum = (-xsq * d_2(xsq)) - d_2(del) + temp.ln();
        if (lower && x > 0.0) || (upper && x <= 0.0) {
            *ccum = log1p(-exp(-xsq * d_2(xsq)) *
                exp(-d_2(del)) * temp);
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
        2.2352520354606839287,
        161.02823106855587881,
        1067.6894854603709582,
        18154.981253343561249,
        0.065682337918207449113
    ];
    let b: [f64; 4] = [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956
    ];
    let c: [f64; 9] = [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8
    ];
    let d: [f64; 8] = [
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727
    ];
    let p: [f64; 6] = [
        0.21589853405795699,
        0.1274011611602473639,
        0.022235277870649807,
        0.001421619193227893466,
        2.9112874951168792e-5,
        0.02307344176494017303
    ];
    let q: [f64; 5] = [
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5
    ];

    let mut xden: f64;
    let mut xnum: f64;
    let mut temp: f64;
    let mut del: f64;
    let mut eps: f64;
    let mut xsq: f64;
    let mut y: f64;

    let min = f64::MIN_POSITIVE;

    let mut i: i32;
    let mut lower: bool;
    let mut upper: bool;

    if x.is_nan() {
        *ccum = x;
        *cum = *ccum;
        return;
    }

    eps = f64::EPSILON * 0.5;

    lower = i_tail != 1;
    upper = i_tail != 0;

    y = x.abs();
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
    } else if log_p && y < 1e170 ||
        lower && -37.5193 < x && x < 8.2924 ||
        upper && -8.2924 < x && x < 37.5193 {

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
