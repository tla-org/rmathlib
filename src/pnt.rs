use crate::dpq::r_dt_0;
use crate::dpq::r_dt_1;
use crate::lgammafn;
use crate::libc::fabs;
use crate::libc::DBL_EPSILON;
use crate::nmath::ml_warn_return_nan;
use crate::nmath::r_finite;
use crate::nmath::M_SQRT_2dPI;
use crate::nmath::DBL_MIN_EXP;
use crate::nmath::M_LN_SQRT_2PI;
use crate::pbeta;
use crate::pnorm;
use crate::pt;
use crate::rmath::M_LN2;
use libm::exp;
use libm::expm1;
use libm::fmin;
use libm::log1p;
use libm::pow;
use libm::sqrt;

/// Non-central t distribution
///
/// Algorithm AS 243  Lenth,R.V. (1989). Appl. Statist., Vol.38, 185-189.
///
/// Cumulative probability at t of the non-central t-distribution
/// with df degrees of freedom (may be fractional) and non-centrality
/// parameter delta.
pub fn pnt(t: f64, df: f64, ncp: f64, lower_tail: bool, log_p: bool) -> f64 {
    // initialize variables
    let mut albeta = 0.0;
    let mut a = 0.0;
    let mut b = 0.0;
    let mut del = 0.0;
    let mut errbd = 0.0;
    let mut lambda = 0.0;
    let mut rxb = 0.0;
    let mut tt = 0.0;
    let mut x = 0.0;
    let mut geven = 0.0;
    let mut godd = 0.0;
    let mut p = 0.0;
    let mut q = 0.0;
    let mut s = 0.0;
    let mut tnc = 0.0;
    let mut xeven = 0.0;
    let mut xodd = 0.0;
    let mut it = 0;
    let mut negdel = false;
    let mut lower_tail = lower_tail;

    // NOTE: itrmax and errmax may be changed to suit one's needs.
    let itrmax = 1_000;
    let errmax = 1e-12;

    if df <= 0.0 {
        return ml_warn_return_nan();
    }

    if ncp == 0.0 {
        return pt(t, df, lower_tail, log_p);
    }

    if !r_finite(t) {
        if t < 0.0 {
            return r_dt_0(lower_tail, log_p);
        } else {
            return r_dt_1(lower_tail, log_p);
        }
    }

    if t >= 0.0 {
        negdel = false;
        tt = -t;
        del = -ncp;
    } else {
        /* We deal quickly with left tail if extreme,
        since pt(q, df, ncp) <= pt(0, df, ncp) = \Phi(-ncp) */
        if ncp > 40.0 && (!log_p || !lower_tail) {
            return r_dt_0(lower_tail, log_p);
        }
        negdel = true;
        tt = -t;
        del = -ncp;
    }

    if df > 4e5 || del * del > 2.0 * M_LN2 * -DBL_MIN_EXP {
        /*-- 2nd part: if del > 37.62, then p=0 below
        FIXME: test should depend on `df', `tt' AND `del' ! */
        /* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */

        let s = 1.0 / (4.0 * df);

        return pnorm(
            tt * (1.0 - s),
            del,
            sqrt(1.0 + tt * tt * 2.0 * s),
            !negdel,
            log_p,
        );
    }

    /* initialize twin series */
    /* Guenther, J. (1978). Statist. Computn. Simuln. vol.6, 199. */

    x = t * t;
    rxb = df / (x + df); /* := (1 - x) {x below} -- but more accurately */
    x = x / (x + df); /* in [0,1) */

    if x > 0.0 {
        // <==>  t != 0
        lambda = del * del;
        p = 0.5 * exp(-0.5 * lambda);
    }

    if p == 0.0 {
        // underflow
        return r_dt_0(lower_tail, log_p);
    }

    q = M_SQRT_2dPI * p * del;
    s = 0.5 * p;
    /* s = 0.5 - p = 0.5*(1 - exp(-.5 L)) =  -0.5*expm1(-.5 L)) */
    if s < 1e-7 {
        s = 0.5 * expm1(-0.5 * lambda);
    }
    a = 0.5;
    b = 0.5 * df;
    /* rxb = (1 - x) ^ b   [ ~= 1 - b*x for tiny x --> see 'xeven' below]
     *       where '(1 - x)' =: rxb {accurately!} above */
    rxb = pow(rxb, b);
    albeta = M_LN_SQRT_2PI + lgammafn(b) - lgammafn(0.5 + b);
    xodd = pbeta(x, a, b, true, false);
    godd = 2.0 * rxb * exp(a * x.ln() - albeta);
    tnc = b * x;
    xeven = if tnc < DBL_EPSILON { tnc } else { 1.0 - rxb };
    geven = tnc * rxb;
    tnc = p * xodd + q * xeven;

    /* repeat until convergence or iteration limit */
    for it in 1..=itrmax {
        a += 1.0;
        xodd -= godd;
        xeven -= geven;
        godd *= x * (a + b - 1.0) / a;
        geven *= x * (a + b - 0.5) / (a + 0.5);
        p *= lambda / (2.0 * it as f64);
        q *= lambda / (2.0 * it as f64 + 1.0);
        tnc += p * xodd + q * xeven;
        s -= p;

        if s < 1e-10 {
            finis(-del, 0.0, 1.0, true, false, &mut tnc);
        }

        if s <= 0.0 && it < 1 {
            finis(-del, 0.0, 1.0, true, false, &mut tnc);
        }

        errbd = 2.0 * s * (xodd - godd);

        if fabs(errbd) < errmax {
            // convergence
            finis(-del, 0.0, 1.0, true, false, &mut tnc);
            break;
        }
    }
    lower_tail = !negdel;
    r_dt_val(fmin(tnc, 1.0), lower_tail, log_p)
}

// converting the goto statement
fn finis(del: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool, tnc: &mut f64) {
    *tnc = pnorm(-del, 0.0, 1.0, true, false);
}

pub fn r_dt_val(x: f64, lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        r_d_val(x, log_p)
    } else {
        r_d_clog(x, lower_tail)
    }
}

pub fn r_d_val(x: f64, log_p: bool) -> f64 {
    if log_p {
        exp(x)
    } else {
        x
    }
}

pub fn r_d_clog(p: f64, log_p: bool) -> f64 {
    if log_p {
        log1p(-p)
    } else {
        0.5 - p + 0.5
    }
}
