use lazy_static::lazy_static;
use libm::lgamma;
use log::warn;
use std::f64::consts::LN_2;
use std::f64::{EPSILON, INFINITY, MAX, MIN, NAN, NEG_INFINITY};
use std::ops::Neg;

use crate::dpq::{r_dt_0, r_dt_1};
use crate::lgamma::lgammafn;

extern "C" {
    pub fn dpois_raw(x: f64, lambda: f64, give_log: bool) -> f64;
    pub fn pnorm5(x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64;
    pub fn dnorm(x: f64, mu: f64, sigma: f64, log_p: bool) -> f64;
    pub fn R_Log1_Exp(x: f64) -> f64;
}

/// This function computes the distribution function for the
/// gamma distribution with shape parameter alph and scale parameter
/// scale. This is also known as the incomplete gamma function.
/// See Abramowitz and Stegun (6.5.1) for example.
///
/// ## NOTES
///
/// Complete redesign by Morten Welinder, originally for Gnumeric.
/// Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
pub fn pgamma(x: f64, alph: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    // Handling special cases
    if x.is_nan() || alph.is_nan() || scale.is_nan() {
        return NAN;
    }

    if alph <= 0.0 || scale <= 0.0 {
        return NAN; // Undefined for non-positive alpha or scale
    }

    let x = x / scale;

    if x.is_nan() {
        // e.g., original x = scale = +Inf
        return x;
    }

    if alph == 0.0 {
        // Limit case; useful e.g., in pnchisq
        return if x <= 0.0 {
            r_dt_0(lower_tail, log_p)
        } else {
            r_dt_1(lower_tail, log_p)
        };
    }

    pgamma_raw(x, alph, lower_tail, log_p)
}

const SQR: fn(f64) -> f64 = |x| x * x;

lazy_static! {
    /// Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77
    static ref SCALEFACTOR: f64 = SQR(SQR(SQR(4294967296.0))); // (2^32)^8
    // static ref SCALEFACTOR: f64 = 4294967296.0f64.powi(24); // (2^32)^8 = 2^256
}

/// If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x
const M_CUTOFF: f64 = LN_2 * MAX as f64 / MIN as f64; // 3.196577e18

/// Continued fraction for calculation of
/// 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
/// auxiliary in log1pmx() and lgamma1p()
fn logcf(x: f64, i: f64, d: f64, eps: f64) -> f64 {
    let mut c1 = 2.0 * d;
    let mut c2 = i + d;
    let mut c4 = c2 + d;
    let mut a1 = c2;
    let mut b1 = i * (c2 - i * x);
    let mut b2 = d * d * x;
    let mut a2 = c4 * c2 - b2;

    b2 = c4 * b1 - i * b2;

    while (a2 * b1 - a1 * b2).abs() > (eps * b1 * b2).abs() {
        let c3 = c2 * c2 * x;
        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;

        let c3 = c1 * c1 * x;
        c1 += d;
        c4 += d;
        a2 = c4 * a1 - c3 * a2;
        b2 = c4 * b1 - c3 * b2;

        if b2.abs() > *SCALEFACTOR {
            a1 /= *SCALEFACTOR;
            b1 /= *SCALEFACTOR;
            a2 /= *SCALEFACTOR;
            b2 /= *SCALEFACTOR;
        } else if b2.abs() < 1.0 / *SCALEFACTOR {
            a1 *= *SCALEFACTOR;
            b1 *= *SCALEFACTOR;
            a2 *= *SCALEFACTOR;
            b2 *= *SCALEFACTOR;
        }
    }

    a2 / b2
}

/// Accurate calculation of log(1+x)-x, particularly for small x.
fn log1pmx(x: f64) -> f64 {
    const MIN_LOG1_VALUE: f64 = -0.79149064;

    if x > 1.0 || x < MIN_LOG1_VALUE {
        (x + 1.0).ln() - x
    } else {
        // -.791 <= x <= 1 -- expand in [x/(2+x)]^2 =: y
        let r = x / (2.0 + x);
        let y = r * r;

        if x.abs() < 1e-2 {
            // Use series expansion
            const TWO: f64 = 2.0;
            r * (((((TWO / 9.0 * y + TWO / 7.0) * y + TWO / 5.0) * y + TWO / 3.0) * y) - x)
        } else {
            // Use continued fraction expansion
            const TOL_LOGCF: f64 = 1e-14;
            r * (2.0 * y * logcf(y, 3.0, 2.0, TOL_LOGCF) - x)
        }
    }
}

/// Compute log(gamma(a+1)) accurately also for small a (0 < a < 0.5).
pub fn lgamma1p(a: f64) -> f64 {
    if a.abs() >= 0.5 {
        // FIXME: port C function
        unsafe { lgamma(a + 1.0) }
    } else {
        const EULERS_CONST: f64 = 0.5772156649015328606065120900824024;
        const COEFFS: [f64; 40] = [
            0.3224670334241132182362075833230126e-0, // = (zeta(2)-1)/2
            0.6735230105319809513324605383715000e-1, // = (zeta(3)-1)/3
            0.2058080842778454787900092413529198e-1,
            0.7385551028673985266273097291406834e-2,
            0.2890510330741523285752988298486755e-2,
            0.1192753911703260977113935692828109e-2,
            0.5096695247430424223356548135815582e-3,
            0.2231547584535793797614188036013401e-3,
            0.9945751278180853371459589003190170e-4,
            0.4492623673813314170020750240635786e-4,
            0.2050721277567069155316650397830591e-4,
            0.9439488275268395903987425104415055e-5,
            0.4374866789907487804181793223952411e-5,
            0.2039215753801366236781900709670839e-5,
            0.9551412130407419832857179772951265e-6,
            0.4492469198764566043294290331193655e-6,
            0.2120718480555466586923135901077628e-6,
            0.1004322482396809960872083050053344e-6,
            0.4769810169363980565760193417246730e-7,
            0.2271109460894316491031998116062124e-7,
            0.1083865921489695409107491757968159e-7,
            0.5183475041970046655121248647057669e-8,
            0.2483674543802478317185008663991718e-8,
            0.1192140140586091207442548202774640e-8,
            0.5731367241678862013330194857961011e-9,
            0.2759522885124233145178149692816341e-9,
            0.1330476437424448948149715720858008e-9,
            0.6422964563838100022082448087644648e-10,
            0.3104424774732227276239215783404066e-10,
            0.1502138408075414217093301048780668e-10,
            0.7275974480239079662504549924814047e-11,
            0.3527742476575915083615072228655483e-11,
            0.1711991790559617908601084114443031e-11,
            0.8315385841420284819798357793954418e-12,
            0.4042200525289440065536008957032895e-12,
            0.1966475631096616490411045679010286e-12,
            0.9573630387838555763782200936508615e-13,
            0.4664076026428374224576492565974577e-13,
            0.2273736960065972320633279596737272e-13,
            0.1109139947083452201658320007192334e-13, // = (zeta(40+1)-1)/(40+1)
        ];
        const TOL_LOGCF: f64 = 1e-14;
        const C: f64 = 0.2273736845824652515226821577978691e-12; // zeta(N+2)-1

        /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
         * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
         * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
         *
         * Here, another convergence acceleration trick is used to compute
         * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
         */
        let mut lgam = C * logcf(-a / 2.0, 40.0 + 2.0, 1.0, TOL_LOGCF);
        for &coeff in COEFFS.iter().rev() {
            lgam = coeff - a * lgam;
        }

        (a * lgam - EULERS_CONST) * a - log1pmx(a)
    }
}

/// Compute the log of a sum from logs of terms, i.e., log(exp(logx) + exp(logy))
/// without causing overflows and without throwing away large handfuls of accuracy.
pub fn logspace_add(logx: f64, logy: f64) -> f64 {
    logx.max(logy) + (-(logx - logy).abs()).exp().ln_1p()
}

/// Compute the log of a difference from logs of terms, i.e., log(exp(logx) - exp(logy))
/// without causing overflows and without throwing away large handfuls of accuracy.
pub fn logspace_sub(logx: f64, logy: f64) -> f64 {
    logx + (logy - logx).exp().ln_1p().neg()
}

/// Compute the log of a sum from logs of terms, i.e.,
/// log(sum_i exp(logx[i])) in a way that avoids overflows.
pub fn logspace_sum(logx: &[f64]) -> f64 {
    match logx.len() {
        0 => NEG_INFINITY, // log(0) for empty input
        1 => logx[0],
        2 => logspace_add(logx[0], logx[1]),
        _ => {
            // Find the maximum log value to scale other values
            let mx = logx.iter().cloned().fold(NEG_INFINITY, f64::max);
            let sum: f64 = logx.iter().map(|&x| (x - mx).exp()).sum();
            mx + sum.ln()
        }
    }
}

/// dpois_wrap (x_plus_1, lambda) := dpois(x_plus_1 - 1, lambda)
/// where dpois(k, L) := exp(-L) L^k / gamma(k+1) {the usual Poisson probabilities}
/// and dpois*(.., give_log = true) := log(dpois*(..))
fn dpois_wrap(x_plus_1: f64, lambda: f64, give_log: bool) -> f64 {
    if !lambda.is_finite() {
        if give_log {
            NEG_INFINITY
        } else {
            0.0
        }
    } else if x_plus_1 > 1.0 {
        // FIXME: port C function
        unsafe { dpois_raw(x_plus_1 - 1.0, lambda, give_log) }
    } else {
        if lambda > (x_plus_1 - 1.0).abs() * M_CUTOFF {
            // FIXME: port C function
            let res = -lambda - unsafe { lgammafn(x_plus_1) };
            if give_log {
                res
            } else {
                res.exp()
            }
        } else {
            // FIXME: port C function
            let d = unsafe { dpois_raw(x_plus_1, lambda, give_log) };
            if give_log {
                d + (x_plus_1 / lambda).ln()
            } else {
                d * (x_plus_1 / lambda)
            }
        }
    }
}

/// Compute the lower tail of the gamma distribution for small x values.
fn pgamma_smallx(x: f64, alph: f64, lower_tail: bool, log_p: bool) -> f64 {
    let mut sum = 0.0;
    let mut c = alph;
    let mut n = 0.0;
    let mut term;

    // Loop to sum up terms until convergence
    loop {
        n += 1.0;
        c *= -x / n;
        term = c / (alph + n);
        sum += term;

        if term.abs() <= EPSILON * sum.abs() {
            break;
        }
    }

    if lower_tail {
        let f1 = if log_p { (1.0 + sum).ln() } else { 1.0 + sum };
        let f2 = if alph > 1.0 {
            let d = dpois_wrap(alph, x, log_p);
            if log_p {
                d + x
            } else {
                d * x.exp()
            }
        } else {
            if log_p {
                alph * x.ln() - lgamma1p(alph)
            } else {
                x.powf(alph) / lgamma1p(alph).exp()
            }
        };
        if log_p {
            f1 + f2
        } else {
            f1 * f2
        }
    } else {
        let lf2 = alph * x.ln() - lgamma1p(alph);
        if log_p {
            (1.0 + sum).ln_1p() + lf2
        } else {
            let f1m1 = sum;
            let f2m1 = lf2.exp_m1();
            -(f1m1 + f2m1 + f1m1 * f2m1)
        }
    }
}

/// Compute the upper tail of the Poisson distribution using a series expansion.
fn pd_upper_series(x: f64, y: f64, log_p: bool) -> f64 {
    let mut term = x / y;
    let mut sum = term;
    let mut y = y;

    while term > sum * std::f64::EPSILON {
        y += 1.0;
        term *= x / y;
        sum += term;
    }

    if log_p {
        sum.ln()
    } else {
        sum
    }
}

/// Compute the lower tail of the Poisson distribution using a series expansion.
fn pd_lower_series(lambda: f64, y: f64) -> f64 {
    let mut term = 1.0;
    let mut sum = 0.0;
    let mut y = y;

    while y >= 1.0 && term > sum * std::f64::EPSILON {
        term *= y / lambda;
        sum += term;
        y -= 1.0;
    }

    if y != y.floor() {
        // For non-integer y, add another term to the sum
        let f = pd_lower_cf(y, lambda + 1.0 - y);
        sum += term * f;
    }

    sum
}

/// Continued fraction for calculation of scaled lower-tail F distribution.
fn pd_lower_cf(y: f64, d: f64) -> f64 {
    const MAX_IT: u32 = 200000;
    let mut f: f64 = 0.0;
    let mut i = 0;
    let mut of = -1.0; // Initialization far away
    let mut c2 = y;
    let mut c4 = d;
    let mut a1 = 0.0;
    let mut b1 = 1.0;
    let mut a2 = y;
    let mut b2 = d;

    while i < MAX_IT {
        i += 1;
        c2 -= 1.0;
        let c3 = i as f64 * c2;
        c4 += 2.0;
        a1 = c4 * a2 + c3 * a1;
        b1 = c4 * b2 + c3 * b1;

        i += 1;
        c2 -= 1.0;
        let c3 = i as f64 * c2;
        c4 += 2.0;
        a2 = c4 * a1 + c3 * a2;
        b2 = c4 * b1 + c3 * b2;

        if b2.abs() > *SCALEFACTOR {
            a1 /= *SCALEFACTOR;
            b1 /= *SCALEFACTOR;
            a2 /= *SCALEFACTOR;
            b2 /= *SCALEFACTOR;
        } else if b2.abs() < 1.0 / *SCALEFACTOR {
            a1 *= *SCALEFACTOR;
            b1 *= *SCALEFACTOR;
            a2 *= *SCALEFACTOR;
            b2 *= *SCALEFACTOR;
        }

        if b2 != 0.0 {
            f = a2 / b2;
            if (f - of).abs() <= std::f64::EPSILON * f.abs() {
                return f;
            }
            of = f;
        }
    }

    warn!("Non-convergence in pd_lower_cf after {MAX_IT} iterations.");
    f // Returning the last computed value of `f` as a fallback.
}

/// Asymptotic expansion to calculate the probability that Poisson variate
/// has value <= x.
fn ppois_asymp(x: f64, lambda: f64, lower_tail: bool, log_p: bool) -> f64 {
    const COEFS_A: [f64; 8] = [
        -1e99, // placeholder for 1-indexing
        2.0 / 3.0,
        -4.0 / 135.0,
        8.0 / 2835.0,
        16.0 / 8505.0,
        -8992.0 / 12629925.0,
        -334144.0 / 492567075.0,
        698752.0 / 1477701225.0,
    ];

    const COEFS_B: [f64; 8] = [
        -1e99, // placeholder
        1.0 / 12.0,
        1.0 / 288.0,
        -139.0 / 51840.0,
        -571.0 / 2488320.0,
        163879.0 / 209018880.0,
        5246819.0 / 75246796800.0,
        -534703531.0 / 902961561600.0,
    ];

    let dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    let pt_ = -(-dfm / x).ln_1p();
    let s2pt = (2.0 * x * pt_).sqrt();
    let s2pt = if dfm < 0.0 { -s2pt } else { s2pt };

    let elfb = (0..8).map(|i| COEFS_B[i] / x.powi(i as i32)).sum::<f64>();
    let elfb_term = x * elfb;

    let res12 = (0..8)
        .fold((0.0, 1.0 / x.sqrt(), s2pt / x), |(acc, term1, term2), i| {
            let res1_ig = term1 / x + COEFS_A[i] * term1.powi(i as i32 + 1);
            let res2_ig = term2 / x + COEFS_B[i] * term2.powi(i as i32 + 1);
            (acc + res1_ig + res2_ig, res1_ig, res2_ig)
        })
        .0;

    let f = res12 / elfb_term;
    let np = unsafe { pnorm5(s2pt, 0.0, 1.0, !lower_tail, log_p) };

    if log_p {
        np + (1.0 + f * unsafe { dnorm(s2pt, 0.0, 1.0, log_p).exp() }).ln()
    } else {
        np + f * unsafe { dnorm(s2pt, 0.0, 1.0, log_p) }
    }
}

/// Compute the lower tail probability of the gamma distribution.
fn pgamma_raw(x: f64, alph: f64, lower_tail: bool, log_p: bool) -> f64 {
    // Check boundaries
    if x.is_nan() || alph.is_nan() || x < 0.0 || alph <= 0.0 {
        return NEG_INFINITY;
    }

    if x == 0.0 {
        return if lower_tail {
            if log_p {
                NEG_INFINITY
            } else {
                0.0
            }
        } else {
            if log_p {
                0.0
            } else {
                1.0
            }
        };
    }

    if alph == INFINITY {
        return if lower_tail {
            if log_p {
                0.0
            } else {
                1.0
            }
        } else {
            if log_p {
                NEG_INFINITY
            } else {
                0.0
            }
        };
    }

    // Main computation
    if x < 1.0 {
        pgamma_smallx(x, alph, lower_tail, log_p)
    } else if x <= alph - 1.0 && x < 0.8 * (alph + 50.0) {
        let sum = pd_upper_series(x, alph, log_p);
        let d = dpois_wrap(alph, x, log_p);
        if !lower_tail {
            if log_p {
                unsafe { R_Log1_Exp(d + sum) }
            } else {
                1.0 - d * sum.exp()
            }
        } else {
            if log_p {
                sum + d
            } else {
                sum.exp() * d
            }
        }
    } else if alph - 1.0 < x && alph < 0.8 * (x + 50.0) {
        let sum = if alph < 1.0 {
            if x * EPSILON > 1.0 - alph {
                1.0 // To avoid division by zero
            } else {
                let f = pd_lower_cf(alph, x - (alph - 1.0)) * x / alph;
                if log_p {
                    f.ln()
                } else {
                    f
                }
            }
        } else {
            let sum = pd_lower_series(x, alph - 1.0);
            if log_p {
                sum.ln_1p()
            } else {
                1.0 + sum
            }
        };

        let d = dpois_wrap(alph, x, log_p);
        if !lower_tail {
            if log_p {
                sum + d
            } else {
                sum * d
            }
        } else {
            if log_p {
                unsafe { R_Log1_Exp(d + sum) }
            } else {
                1.0 - d * sum
            }
        }
    } else {
        // x >= 1 and x fairly near alph
        ppois_asymp(alph - 1.0, x, !lower_tail, log_p)
    }
}
