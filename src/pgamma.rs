use libm::lgamma;
use std::f64::EPSILON;
use std::f64::MAX;
use std::ops::Neg;

use crate::dnorm::dnorm4;
use crate::dpois::dpois_raw;
use crate::dpq::r_dt_0;
use crate::dpq::r_dt_1;
use crate::lgamma::lgammafn;
use crate::nmath::r_log1_exp;
use crate::nmath::ML_NAN;
use crate::nmath::ML_NEGINF;
use crate::nmath::ML_POSINF;
use crate::pnorm::pnorm5;
use crate::rmath::M_LN2;

/// Computes the distribution function for the gamma distribution
/// with shape parameter alph and scale parameter scale.
///
/// This is also known as the incomplete gamma function.
/// See Abramowitz and Stegun (6.5.1) for example.
///
/// ## NOTES
///
/// Complete redesign by Morten Welinder, originally for Gnumeric.
/// Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
///
/// ## AUTHORS
///
/// 2006-2019 The R Core Team
/// 2005-6 Morten Welinder <terra@gnome.org>
/// 2005-10 The R Foundation
pub fn pgamma(x: f64, alph: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    // Handling special cases
    if x.is_nan() || alph.is_nan() || scale.is_nan() {
        return ML_NAN;
    }

    if alph <= 0.0 || scale <= 0.0 {
        return ML_NAN; // Undefined for non-positive alpha or scale
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

/// If |x| > |k| * M_cutoff,  then  log\[ exp(-x) * k^x \] =~= -x
const M_CUTOFF: f64 = M_LN2 * MAX / EPSILON; // 3.196577e18

/// Continued fraction for calculation of
/// 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
/// auxiliary in log1pmx() and lgamma1p()
fn logcf(x: f64, i: f64, d: f64, eps: f64) -> f64 {
    // Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77
    #![allow(non_snake_case)]
    let SCALEFACTOR: f64 = SQR(SQR(SQR(4294967296.0))); // (2^32)^8

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

        if b2.abs() > SCALEFACTOR {
            a1 /= SCALEFACTOR;
            b1 /= SCALEFACTOR;
            a2 /= SCALEFACTOR;
            b2 /= SCALEFACTOR;
        } else if b2.abs() < 1.0 / SCALEFACTOR {
            a1 *= SCALEFACTOR;
            b1 *= SCALEFACTOR;
            a2 *= SCALEFACTOR;
            b2 *= SCALEFACTOR;
        }
    }

    a2 / b2
}

/// Accurate calculation of log(1+x)-x, particularly for small x.
pub fn log1pmx(x: f64) -> f64 {
    const MIN_LOG1_VALUE: f64 = -0.79149064;

    if !(MIN_LOG1_VALUE..=1.0).contains(&x) {
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
        lgamma(a + 1.0)
    } else {
        const EULERS_CONST: f64 = 0.577_215_664_901_532_9;
        const COEFFS: [f64; 40] = [
            3.224_670_334_241_132e-1, // = (zeta(2)-1)/2
            6.735_230_105_319_81e-2,  // = (zeta(3)-1)/3
            2.058_080_842_778_454_6e-2,
            7.385_551_028_673_986e-3,
            2.890_510_330_741_523_4e-3,
            1.192_753_911_703_261e-3,
            5.096_695_247_430_425e-4,
            2.231_547_584_535_793_9e-4,
            9.945_751_278_180_853e-5,
            4.492_623_673_813_314e-5,
            2.050_721_277_567_069e-5,
            9.439_488_275_268_397e-6,
            4.374_866_789_907_488e-6,
            2.039_215_753_801_366e-6,
            9.551_412_130_407_42e-7,
            4.492_469_198_764_566e-7,
            2.120_718_480_555_466_5e-7,
            1.004_322_482_396_809_9e-7,
            4.769_810_169_363_980_4e-8,
            2.271_109_460_894_316_4e-8,
            1.083_865_921_489_695_5e-8,
            5.183_475_041_970_047e-9,
            2.483_674_543_802_478_5e-9,
            1.192_140_140_586_091_2e-9,
            5.731_367_241_678_862e-10,
            2.759_522_885_124_233_4e-10,
            1.330_476_437_424_449e-10,
            6.422_964_563_838_1e-11,
            3.104_424_774_732_227_6e-11,
            1.502_138_408_075_414_2e-11,
            7.275_974_480_239_079e-12,
            3.527_742_476_575_915e-12,
            1.711_991_790_559_618e-12,
            8.315_385_841_420_285e-13,
            4.042_200_525_289_44e-13,
            1.966_475_631_096_616_5e-13,
            9.573_630_387_838_556e-14,
            4.664_076_026_428_374_4e-14,
            2.273_736_960_065_972_4e-14,
            1.109_139_947_083_452_2e-14, // = (zeta(40+1)-1)/(40+1)
        ];
        const TOL_LOGCF: f64 = 1e-14;
        const C: f64 = 2.273_736_845_824_652_4e-13; // zeta(N+2)-1

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
#[allow(dead_code)]
fn logspace_sub(logx: f64, logy: f64) -> f64 {
    logx + (logy - logx).exp().ln_1p().neg()
}

/// Compute the log of a sum from logs of terms, i.e.,
/// log(sum_i exp(logx\[i\])) in a way that avoids overflows.
#[allow(dead_code)]
fn logspace_sum(logx: &[f64]) -> f64 {
    match logx.len() {
        0 => ML_NEGINF, // log(0) for empty input
        1 => logx[0],
        2 => logspace_add(logx[0], logx[1]),
        _ => {
            // Find the maximum log value to scale other values
            let mx = logx.iter().cloned().fold(ML_NEGINF, f64::max);
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
            ML_NEGINF
        } else {
            0.0
        }
    } else if x_plus_1 > 1.0 {
        dpois_raw(x_plus_1 - 1.0, lambda, give_log)
    } else if lambda > (x_plus_1 - 1.0).abs() * M_CUTOFF {
        let res = -lambda - lgammafn(x_plus_1);
        if give_log {
            res
        } else {
            res.exp()
        }
    } else {
        let d = dpois_raw(x_plus_1, lambda, give_log);
        if give_log {
            d + (x_plus_1 / lambda).ln()
        } else {
            d * (x_plus_1 / lambda)
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
        } else if log_p {
            alph * x.ln() - lgamma1p(alph)
        } else {
            x.powf(alph) / lgamma1p(alph).exp()
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

    while term > sum * EPSILON {
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

    while y >= 1.0 && term > sum * EPSILON {
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
    // Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77
    #![allow(non_snake_case)]
    let SCALEFACTOR: f64 = SQR(SQR(SQR(4294967296.0))); // (2^32)^8
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

        if b2.abs() > SCALEFACTOR {
            a1 /= SCALEFACTOR;
            b1 /= SCALEFACTOR;
            a2 /= SCALEFACTOR;
            b2 /= SCALEFACTOR;
        } else if b2.abs() < 1.0 / SCALEFACTOR {
            a1 *= SCALEFACTOR;
            b1 *= SCALEFACTOR;
            a2 *= SCALEFACTOR;
            b2 *= SCALEFACTOR;
        }

        if b2 != 0.0 {
            f = a2 / b2;
            if (f - of).abs() <= EPSILON * f.abs() {
                return f;
            }
            of = f;
        }
    }

    println!(
        "Non-convergence in pd_lower_cf after {} iterations.",
        MAX_IT
    );
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
    let np = pnorm5(s2pt, 0.0, 1.0, !lower_tail, log_p);

    if log_p {
        np + (1.0 + f * dnorm4(s2pt, 0.0, 1.0, log_p).exp()).ln()
    } else {
        np + f * dnorm4(s2pt, 0.0, 1.0, log_p)
    }
}

/// Compute the lower tail probability of the gamma distribution.
fn pgamma_raw(x: f64, alph: f64, lower_tail: bool, log_p: bool) -> f64 {
    // Check boundaries
    if x.is_nan() || alph.is_nan() || x < 0.0 || alph <= 0.0 {
        return ML_NEGINF;
    }

    if x == 0.0 {
        return if lower_tail {
            if log_p {
                ML_NEGINF
            } else {
                0.0
            }
        } else if log_p {
            0.0
        } else {
            1.0
        };
    }

    if alph == ML_POSINF {
        return if lower_tail {
            if log_p {
                0.0
            } else {
                1.0
            }
        } else if log_p {
            ML_NEGINF
        } else {
            0.0
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
                r_log1_exp(d + sum)
            } else {
                1.0 - d * sum.exp()
            }
        } else if log_p {
            sum + d
        } else {
            sum.exp() * d
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
        } else if log_p {
            r_log1_exp(d + sum)
        } else {
            1.0 - d * sum
        }
    } else {
        // x >= 1 and x fairly near alph
        ppois_asymp(alph - 1.0, x, !lower_tail, log_p)
    }
}
