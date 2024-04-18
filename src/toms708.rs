#![allow(dead_code)]
#![allow(clippy::manual_range_contains)]

use crate::i1mach::i1mach;
use libm::exp;
use libm::expm1;
use libm::fabs;
use libm::log;
use libm::log1p;
use libm::pow;
use libm::sqrt;

const ML_NEGINF: f64 = f64::NEG_INFINITY;
const M_LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;
#[allow(clippy::approx_constant)]
const M_LN2: f64 = 0.693_147_180_559_945_3;
#[allow(clippy::approx_constant)]
const M_LOG10_2: f64 = 0.301_029_995_663_981_2;
const M_SQRT_PI: f64 = 1.772_453_850_905_516;

const DBL_MIN: f64 = f64::MIN;
const DBL_MAX: f64 = f64::MAX;
const DBL_EPSILON: f64 = f64::EPSILON;

// Be careful when replacing these min and max with libm functions.
// The definition below is how they are defined in the original code.
// Other implementations may handle NaN values differently.
fn min(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}

fn max(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else {
        b
    }
}

fn r_log1_exp(x: f64) -> f64 {
    if x > -M_LN2 {
        log(-rexpm1(x))
    } else {
        log1p(-exp(x))
    }
}

/// Based on d1mach.f90.
/// https://github.com/certik/fortran-utils/blob/master/src/legacy/amos/d1mach.f90.
fn d1mach(i: i32) -> f64 {
    match i {
        1 => DBL_MIN,
        2 => DBL_MAX,
        3 => 0.5 * DBL_EPSILON,
        4 => DBL_EPSILON,
        5 => M_LOG10_2,
        _ => 0.0,
    }
}

fn r_d_0(log_p: bool) -> f64 {
    if log_p {
        ML_NEGINF
    } else {
        0.0
    }
}

fn r_d_1(log_p: bool) -> f64 {
    if log_p {
        0.0
    } else {
        1.0
    }
}

fn r_d_exp(log_p: bool, x: f64) -> f64 {
    if log_p {
        x
    } else {
        exp(x)
    }
}

fn fmax2(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        x + y
    } else if x < y {
        y
    } else {
        x
    }
}

fn logspace_add(logx: f64, logy: f64) -> f64 {
    fmax2(logx, logy) + log1p(exp(-fabs(logx - logy)))
}

fn l_end(w: &mut f64, w1: &mut f64, do_swap: bool) {
    #[allow(clippy::manual_swap)]
    if do_swap {
        let t = *w;
        *w = *w1;
        *w1 = t;
    }
}

fn l_end_from_w(w: &mut f64, w1: &mut f64, do_swap: bool, log_p: bool) {
    if log_p {
        *w1 = log1p(-*w);
        *w = log1p(-*w);
    } else {
        *w1 = 0.5 - *w + 0.5;
    }
    l_end(w, w1, do_swap);
}

fn l_end_from_w1(w: &mut f64, w1: &mut f64, do_swap: bool, log_p: bool) {
    if log_p {
        *w = log1p(-*w1);
        *w1 = log(-*w1);
    } else {
        *w = 0.5 - *w1 + 0.5;
    }
    l_end(w, w1, do_swap);
}

fn l_end_from_w1_log(w: &mut f64, w1: &mut f64, do_swap: bool, log_p: bool) {
    // *w1 = log(w1) already; w = 1 - w1  ==> log(w) = log(1 - w1) = log(1 -
    // exp(*w1))
    if log_p {
        *w = r_log1_exp(*w1);
    } else {
        *w = /* 1 - exp(*w1) */ -expm1(*w1);
        *w1 = exp(*w1);
    }
    l_end(w, w1, do_swap);
}

/// Power SERies expansion for evaluating I_x(a,b) when b <= 1 or b*x <= 0.7.
///
/// eps is the tolerance used.
/// NB: if log_p is true, also use it if   (b < 40  & lambda > 650)
fn bpser(a: f64, b: f64, x: f64, eps: f64, log_p: bool) -> f64 {
    let m: i32;
    let mut ans: f64;
    let mut c: f64;
    let t: f64;
    let mut u: f64;
    let z: f64;
    let mut b0: f64;
    let apb: f64;

    if x == 0.0 {
        return r_d_0(log_p);
    }
    // Compute the factor x^a/(a*Beta(a,b)).
    let a0: f64 = min(a, b);
    if a0 >= 1.0 {
        /* ------ 1 <= a0 <= b0 ------ */
        z = a * log(x) - betaln(a, b);
        ans = if log_p { z - log(a) } else { exp(z) / a };
    } else {
        b0 = max(a, b);

        if b0 < 8.0 {
            if b0 <= 1.0 {
                /* ------ a0 < 1 and b0 <= 1 ------ */
                if log_p {
                    ans = a * log(x);
                } else {
                    ans = pow(x, a);
                    if ans == 0.0 {
                        /* once underflow, always underflow .. */
                        return ans;
                    }
                }
                apb = a + b;
                if apb > 1.0 {
                    u = a + b - 1.0;
                    z = (gam1(u) + 1.0) / apb;
                } else {
                    z = gam1(apb) + 1.0;
                }
                c = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z;

                if log_p {
                    /* FIXME ? -- improve quite a bit for c ~= 1 */
                    ans += log(c * (b / apb));
                } else {
                    ans *= c * (b / apb);
                }
            } else {
                /* ------ a0 < 1 < b0 < 8 ------ */

                u = gamln1(a0);
                m = (b0 - 1.0) as i32;
                if m >= 1 {
                    c = 1.0;
                    for _i in 1..=m {
                        b0 += -1.0;
                        c *= b0 / (a0 + b0);
                    }
                    u += log(c);
                }

                z = a * log(x) - u;
                b0 += -1.; // => b0 in (0, 7)
                apb = a0 + b0;
                if apb > 1.0 {
                    u = a0 + b0 - 1.0;
                    t = (gam1(u) + 1.0) / apb;
                } else {
                    t = gam1(apb) + 1.0;
                }

                if log_p {
                    /* FIXME? potential for improving log(t) */
                    ans = z + log(a0 / a) + log1p(gam1(b0)) - log(t);
                } else {
                    ans = exp(z) * (a0 / a) * (gam1(b0) + 1.0) / t;
                }
            }
        } else {
            /* ------  a0 < 1 < 8 <= b0  ------ */

            u = gamln1(a0) + algdiv(a0, b0);
            z = a * log(x) - u;

            if log_p {
                ans = z + log(a0 / a);
            } else {
                ans = a0 / a * exp(z);
            }
        }
    }
    // R_ifDEBUG_printf(" bpser(a=%g, b=%g, x=%g, log=%d): prelim.ans = %.14g;\n", a,
    // b, x, log_p, ans);
    if ans == r_d_0(log_p) || (!log_p && a <= eps * 0.1) {
        return ans;
    }

    // Compute the series.
    let tol: f64 = eps / a;
    let mut n = 0.0;
    let mut sum = 0.0;
    let mut w: f64;
    c = 1.0;
    loop {
        // sum is alternating as long as n < b (<==> 1 - b/n < 0)
        n += 1.0;
        c *= (0.5 - b / n + 0.5) * x;
        w = c / (a + n);
        sum += w;
        if !(n < 1e7 && fabs(w) > tol) {
            break;
        }
    }
    if fabs(w) > tol {
        // the series did not converge (in time)
        // warn only when the result seems to matter:
        if (log_p && !(a * sum > -1.0 && fabs(log1p(a * sum)) < eps * fabs(ans)))
            || (!log_p && fabs(a * sum + 1.0) != 1.0)
        {
            // printf(" bpser(a=%g, b=%g, x=%g,...) did not converge (n=1e7, |w|/tol=%g "
            //       "> 1; A=%g)",
            //       a, b, x, fabs(w) / tol, ans);
        }
    }
    // R_ifDEBUG_printf(
    //    "  -> n=%.0f iterations, |w|=%g %s %g=tol:=eps/a ==> a*sum=%g\n", n,
    //    fabs(w), (fabs(w) > tol) ? ">!!>" : "<=", tol, a * sum);
    if log_p {
        if a * sum > -1.0 {
            ans += log1p(a * sum);
        } else {
            if ans > ML_NEGINF {
                // printf("pbeta(*, log.p=true) -> bpser(a=%g, b=%g, x=%g,...) underflow "
                //       "to -Inf",
                //       a, b, x);
            }
            ans = ML_NEGINF;
        }
    } else if a * sum > -1.0 {
        ans *= a * sum + 1.0;
    } else {
        // underflow to
        ans = 0.0;
    }
    ans
}

#[allow(clippy::too_many_arguments)]
fn l_w_bpser(
    a0: f64,
    b0: f64,
    x0: f64,
    w: &mut f64,
    w1: &mut f64,
    eps: f64,
    do_swap: bool,
    log_p: bool,
) {
    *w = bpser(a0, b0, x0, eps, log_p);
    *w1 = if log_p {
        r_log1_exp(*w)
    } else {
        0.5 - *w + 0.5
    };
    // R_ifDEBUG_printf(" L_w_bpser: *w := bpser(*) = %.15g\n", *w);
    l_end(w, w1, do_swap);
}

#[allow(clippy::too_many_arguments)]
fn l_w1_bpser(
    a0: f64,
    b0: f64,
    y0: f64,
    w: &mut f64,
    w1: &mut f64,
    eps: f64,
    do_swap: bool,
    log_p: bool,
) {
    *w1 = bpser(b0, a0, y0, eps, log_p);
    *w = if log_p {
        r_log1_exp(*w1)
    } else {
        0.5 - *w1 + 0.5
    };
    // R_ifDEBUG_printf(" L_w1_bpser: *w1 := bpser(*) = %.15g\n", *w1);
    l_end(w, w1, do_swap)
}

#[allow(clippy::too_many_arguments)]
fn l_bfrac(
    a0: f64,
    b0: f64,
    x0: f64,
    y0: f64,
    lambda: f64,
    eps: f64,
    w: &mut f64,
    w1: &mut f64,
    do_swap: bool,
    log_p: bool,
) {
    *w = bfrac(a0, b0, x0, y0, lambda, eps * 15., log_p);
    *w1 = if log_p {
        r_log1_exp(*w)
    } else {
        0.5 - *w + 0.5
    };
    // R_ifDEBUG_printf(" L_bfrac: *w := bfrac(*) = %g\n", *w);
    l_end(w, w1, do_swap)
}

#[allow(clippy::too_many_arguments)]
/// Evaluates the appropriate algorithm.
fn l140(
    mut a0: f64,
    mut b0: f64,
    x0: f64,
    y0: f64,
    eps: f64,
    w: &mut f64,
    w1: &mut f64,
    do_swap: bool,
    ierr: &mut i32,
    mut ierr1: i32,
    log_p: bool,
) {
    let mut n: i32 = b0 as i32;
    b0 -= n as f64;
    if b0 == 0.0 {
        n -= 1;
        b0 = 1.;
    }

    *w = bup(b0, a0, y0, x0, n, eps, false);

    if *w < DBL_MIN && log_p {
        /* do not believe it; try bpser() : */
        // R_ifDEBUG_printf(" L140: bup(b0=%g,..)=%.15g < DBL_MIN - not used; ", b0,
        //                 *w);
        /*revert: */
        b0 += n as f64;
        /* which is only valid if b0 <= 1 || b0*x0 <= 0.7 */
        return l_w_bpser(a0, b0, x0, w, w1, eps, do_swap, log_p);
    }
    // R_ifDEBUG_printf(" L140: *w := bup(b0=%g,..) = %.15g; ", b0, *w);
    if x0 <= 0.7 {
        /* log_p :  TODO:  w = bup(.) + bpser(.)  -- not so easy to use log-scale */
        *w += bpser(a0, b0, x0, eps, /* log_p = */ false);
        // R_ifDEBUG_printf(" x0 <= 0.7: *w := *w + bpser(*) = %.15g\n", *w);
        return l_end_from_w(w, w1, do_swap, log_p);
    }
    /* L150: */
    if a0 <= 15.0 {
        n = 20;
        *w += bup(a0, b0, x0, y0, n, eps, false);
        // R_ifDEBUG_printf("\n a0 <= 15: *w := *w + bup(*) = %.15g;", *w);
        a0 += n as f64;
    }
    // R_ifDEBUG_printf(" bgrat(*, w=%.15g) ", *w);
    bgrat(a0, b0, x0, y0, w, 15.0 * eps, &mut ierr1, false);
    if ierr1 != 0 {
        *ierr = 10 + ierr1;
    }
    // #ifdef DEBUG_bratio
    //  REprintf("==> new w=%.15g", *w);
    //  if (ierr1) {
    //    REprintf(" Error(code=%d)\n", ierr1);
    //  } else {
    //    REprintf("\n");
    //  }
    // #endif
    l_end_from_w(w, w1, do_swap, log_p)
}

#[allow(clippy::too_many_arguments)]
fn l131(
    // These variables are used for debugging in original code.
    _a: f64,
    _b: f64,
    _x: f64,
    n: i32,
    a0: f64,
    b0: f64,
    x0: f64,
    y0: f64,
    w: &mut f64,
    w1: &mut f64,
    eps: f64,
    ierr: &mut i32,
    mut ierr1: i32,
    did_bup: bool,
    do_swap: bool,
    log_p: bool,
) {
    // R_ifDEBUG_printf(" L131: bgrat(*, w1=%.15g) ", *w1);
    bgrat(b0, a0, y0, x0, w1, 15.0 * eps, &mut ierr1, false);
    // #ifdef DEBUG_bratio
    //   REprintf(" ==> new w1=%.15g", *w1);
    //   if (ierr1) {
    //     REprintf(" ERROR(code=%d)\n", ierr1);
    //   } else {
    //     REprintf("\n");
    //   }
    // #endif
    #[allow(clippy::impossible_comparisons)]
    if *w1 == 0.0 || (0.0 < *w1 && *w1 < DBL_MIN) {
        // w1=0 or very close:
        // "almost surely" from underflow, try more: [2013-03-04]
        // FIXME: it is even better to do this in bgrat *directly* at least for the case
        //  !did_bup, i.e., where *w1 = (0 or -Inf) on entry
        // R_ifDEBUG_printf(" denormalized or underflow (?) -> retrying: ");
        if did_bup {
            // re-do that part on log scale:
            *w1 = bup(b0 - (n as f64), a0, y0, x0, n, eps, true);
        } else {
            *w1 = ML_NEGINF; // = 0 on log-scale
        }
        bgrat(b0, a0, y0, x0, w1, 15.0 * eps, &mut ierr1, true);
        if ierr1 != 0 {
            *ierr = 10 + ierr1;
        }
        // #ifdef DEBUG_bratio
        //     REprintf(" ==> new log(w1)=%.15g", *w1);
        //     if (ierr1) {
        //       REprintf(" Error(code=%d)\n", ierr1);
        //     } else {
        //       REprintf("\n");
        //     }
        // #endif
        return l_end_from_w1_log(w, w1, do_swap, log_p);
    }
    // else
    if ierr1 != 0 {
        *ierr = 10 + ierr1;
    }
    if *w1 < 0.0 {
        // printf("bratio(a=%g, b=%g, x=%g): bgrat() -> w1 = %g", a, b, x, *w1);
    }
    l_end_from_w1(w, w1, do_swap, log_p)
}

#[allow(clippy::too_many_arguments)]
fn bratio_final_else(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    eps: f64,
    w: &mut f64,
    w1: &mut f64,
    ierr: &mut i32,
    ierr1: i32,
    log_p: bool,
) {
    /* L30: -------------------- both  a, b > 1  {a0 > 1  &  b0 > 1}
    ---*/

    /* lambda := a y - b x  =  (a + b)y  =  a - (a+b)x    {using x + y == 1},
     * ------ using the numerically best version : */
    let mut lambda: f64 = if (a + b).is_finite() {
        if a > b {
            (a + b) * y - b
        } else {
            a - (a + b) * x
        }
    } else {
        a * y - b * x
    };
    let do_swap = lambda < 0.0;
    let a0: f64;
    let b0: f64;
    let x0: f64;
    let y0: f64;
    if do_swap {
        lambda = -lambda;
        a0 = b;
        x0 = y;
        b0 = a;
        y0 = x;
    } else {
        a0 = a;
        x0 = x;
        b0 = b;
        y0 = y;
    }

    // R_ifDEBUG_printf("  L30:  both  a, b > 1; |lambda| = %#g, do_swap = %d\n",
    //                 lambda, do_swap);

    if b0 < 40.0 {
        // R_ifDEBUG_printf("  b0 < 40;");
        if b0 * x0 <= 0.7 || (log_p && lambda > 650.) {
            // << added 2010-03; svn r51327
            return l_w_bpser(a0, b0, x0, w, w1, eps, do_swap, log_p);
        } else {
            return l140(a0, b0, x0, y0, eps, w, w1, do_swap, ierr, ierr1, log_p);
        }
    } else if a0 > b0 {
        /* ----  a0 > b0 >= 40  ---- */
        // R_ifDEBUG_printf("  a0 > b0 >= 40;");
        if b0 <= 100. || lambda > b0 * 0.03 {
            return l_bfrac(a0, b0, x0, y0, lambda, eps, w, w1, do_swap, log_p);
        }
    } else if a0 <= 100.0 {
        // R_ifDEBUG_printf("  a0 <= 100; a0 <= b0 >= 40;");
        return l_bfrac(a0, b0, x0, y0, lambda, eps, w, w1, do_swap, log_p);
    } else if lambda > a0 * 0.03 {
        // R_ifDEBUG_printf("  b0 >= a0 > 100; lambda > a0 * 0.03 ");
        return l_bfrac(a0, b0, x0, y0, lambda, eps, w, w1, do_swap, log_p);
    }

    /* else if none of the above    L180: */
    *w = basym(a0, b0, lambda, eps * 100.0, log_p);
    *w1 = if log_p {
        r_log1_exp(*w)
    } else {
        0.5 - *w + 0.5
    };
    // R_ifDEBUG_printf(
    //     "  b0 >= a0 > 100; lambda <= a0 * 0.03: *w:= basym(*) =%.15g\n", *w);
    l_end(w, w1, do_swap)
}

#[allow(clippy::too_many_arguments)]
/// Evaluation of the Incomplete Beta function I_x(a,b)
///
/// It is assumed that a and b are nonnegative, and that x <= 1
/// and y = 1 - x.  Bratio assigns w and w1 the values
///
/// w  = I_x(a,b)
/// w1 = 1 - I_x(a,b)
///
/// ierr is a variable that reports the status of the results.
/// If no input errors are detected then ierr is set to 0 and
/// w and w1 are computed. otherwise, if an error is detected,
/// then w and w1 are assigned the value 0 and ierr is set to
/// one of the following values ...
///
///  ierr = 1  if a or b is negative
///  ierr = 2  if a = b = 0
///  ierr = 3  if x < 0 or x > 1
///  ierr = 4  if y < 0 or y > 1
///  ierr = 5  if x + y != 1
///  ierr = 6  if x = a = 0
///  ierr = 7  if y = b = 0
///  ierr = 8  (not used currently)
///  ierr = 9  NaN in a, b, x, or y
///  ierr = 10     (not used currently)
///  ierr = 11  bgrat() error code 1 [+ warning in bgrat()]
///  ierr = 12  bgrat() error code 2   (no warning here)
///  ierr = 13  bgrat() error code 3   (no warning here)
///  ierr = 14  bgrat() error code 4 [+ WARNING in bgrat()]
///
///
/// -------------------
///    Written by Alfred H. Morris, Jr.
///  Naval Surface Warfare Center
///  Dahlgren, Virginia
///    Revised ... Nov 1991
pub fn bratio(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    w: &mut f64,
    w1: &mut f64,
    ierr: &mut i32,
    log_p: bool,
) {
    let do_swap: bool;
    // n used to be not initialized here, but that meant it was used uninitialized
    // when going through GOTO L131.
    let mut n: i32 = 0;
    let ierr1: i32 = 0;
    let a0: f64;
    let mut b0: f64;
    let x0: f64;
    let y0: f64;

    /*  eps is a machine dependent constant: the smallest
     *      floating point number for which   1. + eps > 1.
     * NOTE: for almost all purposes it is replaced by 1e-15 (~= 4.5 times larger)
     * below */
    let mut eps: f64 = 2.0 * d1mach(3); /* == DBL_EPSILON (in R, Rmath) */

    /* ----------------------------------------------------------------------- */
    *w = r_d_0(log_p);
    *w1 = r_d_0(log_p);

    // safeguard, preventing infinite loops further down
    if x.is_nan() || y.is_nan() || a.is_nan() || b.is_nan() {
        *ierr = 9;
        return;
    }

    if a < 0.0 || b < 0.0 {
        *ierr = 1;
        return;
    }
    if a == 0.0 && b == 0.0 {
        *ierr = 2;
        return;
    }
    if x < 0.0 || x > 1.0 {
        *ierr = 3;
        return;
    }
    if y < 0.0 || y > 1.0 {
        *ierr = 4;
        return;
    }

    /* check that  'y == 1 - x' : */
    let z: f64 = x + y - 0.5 - 0.5;

    if fabs(z) > eps * 3.0 {
        *ierr = 5;
        return;
    }

    // R_ifDEBUG_printf("bratio(a=%g, b=%g, x=%9g, y=%9g, .., log_p=%d): ", a, b, x,
    //                 y, log_p);
    *ierr = 0;
    if x == 0.0 {
        if a == 0.0 {
            *ierr = 6;
            return;
        } else {
            *w = r_d_0(log_p);
            *w1 = r_d_1(log_p);
            return;
        }
    }

    if y == 0.0 {
        if b == 0.0 {
            *ierr = 7;
            return;
        } else {
            *w = r_d_1(log_p);
            *w1 = r_d_0(log_p);
            return;
        }
    }

    if a == 0.0 {
        *w = r_d_1(log_p);
        *w1 = r_d_0(log_p);
        return;
    }
    if b == 0.0 {
        *w = r_d_0(log_p);
        *w1 = r_d_1(log_p);
        return;
    }

    eps = max(eps, 1e-15);
    let a_lt_b: bool = a < b;
    if
    /* max(a,b) */
    (if a_lt_b { b } else { a } < eps * 0.001) {
        /* procedure for a and b < 0.001 * eps */
        // L230:  -- result *independent* of x (!)
        // *w  = a/(a+b)  and  w1 = b/(a+b) :
        if log_p {
            if a_lt_b {
                *w = log1p(-a / (a + b)); // notably if a << b
                *w1 = log(a / (a + b));
            } else {
                // b <= a
                *w = log(b / (a + b));
                *w1 = log1p(-b / (a + b));
            }
        } else {
            *w = b / (a + b);
            *w1 = a / (a + b);
        }

        // R_ifDEBUG_printf("a & b very small -> simple ratios (%g,%g)\n", *w, *w1);
        return;
    }

    if min(a, b) <= 1.0 {
        /*------------------------ a <= 1  or  b <= 1 ---- */

        do_swap = x > 0.5;
        if do_swap {
            a0 = b;
            x0 = y;
            b0 = a;
            y0 = x;
        } else {
            a0 = a;
            x0 = x;
            b0 = b;
            y0 = y;
        }
        /* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */

        // R_ifDEBUG_printf(" min(a,b) <= 1, do_swap=%d;", do_swap);

        if b0 < min(eps, eps * a0) {
            /* L80: */
            *w = fpser(a0, b0, x0, eps, log_p);
            *w1 = if log_p {
                r_log1_exp(*w)
            } else {
                0.5 - *w + 0.5
            };
            // R_ifDEBUG_printf("  b0 small -> w := fpser(*) = %.15g\n", *w);
            return l_end(w, w1, do_swap);
        }

        if a0 < min(eps, eps * b0) && b0 * x0 <= 1.0 {
            /* L90: */
            *w1 = apser(a0, b0, x0, eps);
            // R_ifDEBUG_printf("  a0 small -> w1 := apser(*) = %.15g\n", *w1);
            return l_end_from_w1(w, w1, do_swap, log_p);
        }

        let mut did_bup = false;
        if max(a0, b0) > 1.0 {
            /* L20:  min(a,b) <= 1 < max(a,b)  */
            // R_ifDEBUG_printf("\n L20:  min(a,b) <= 1 < max(a,b); ");
            if b0 <= 1.0 {
                return l_w_bpser(a0, b0, x0, w, w1, eps, do_swap, log_p);
            }

            if x0 >= 0.29 {
                /* was 0.3, PR#13786 */
                return l_w1_bpser(a0, b0, y0, w, w1, eps, do_swap, log_p);
            }

            if x0 < 0.1 && pow(x0 * b0, a0) <= 0.7 {
                return l_w_bpser(a0, b0, x0, w, w1, eps, do_swap, log_p);
            }

            if b0 > 15.0 {
                *w1 = 0.;
                return l131(
                    a, b, x, n, a0, b0, x0, y0, w, w1, eps, ierr, ierr1, did_bup, do_swap, log_p,
                );
            }
        } else {
            /*  a, b <= 1 */
            // R_ifDEBUG_printf("\n      both a,b <= 1; ");
            if a0 >= min(0.2, b0) {
                return l_w_bpser(a0, b0, x0, w, w1, eps, do_swap, log_p);
            }
            if pow(x0, a0) <= 0.9 {
                return l_w_bpser(a0, b0, x0, w, w1, eps, do_swap, log_p);
            }
            if x0 >= 0.3 {
                return l_w1_bpser(a0, b0, y0, w, w1, eps, do_swap, log_p);
            }
        }
        n = 20; /* goto L130; */
        *w1 = bup(b0, a0, y0, x0, n, eps, false);
        did_bup = true;
        // R_ifDEBUG_printf("  ... n=20 and *w1 := bup(*) = %.15g; ", *w1);
        b0 += n as f64;
        l131(
            a, b, x, n, a0, b0, x0, y0, w, w1, eps, ierr, ierr1, did_bup, do_swap, log_p,
        )
    } else {
        /* L30: -------------------- both  a, b > 1  {a0 > 1  &  b0 > 1}
        ---*/
        // Extracted this part into a separate function to reduce complexity.
        bratio_final_else(a, b, x, y, eps, w, w1, ierr, ierr1, log_p)
    } /* else: a, b > 1 */
}

#[allow(clippy::too_many_arguments)]
/// Evaluates I_x(a,b) for b < min(eps, eps*a) and x <= 0.5.
fn fpser(a: f64, b: f64, x: f64, eps: f64, log_p: bool) -> f64 {
    let mut ans: f64;
    let mut c: f64;
    let mut s: f64;
    let mut t: f64;
    let mut an: f64;

    // Set ans := x^a.
    if log_p {
        ans = a * log(x);
    } else if a > eps * 0.001 {
        t = a * log(x);
        if t < exparg(1) {
            /* exp(t) would underflow */
            return 0.0;
        }
        ans = exp(t);
    } else {
        ans = 1.0;
    }

    /* Note that 1/b(a,b) = b */

    if log_p {
        ans += log(b) - log(a);
    } else {
        ans *= b / a;
    }

    let tol: f64 = eps / a;
    an = a + 1.;
    t = x;
    s = t / an;
    loop {
        an += 1.;
        t *= x;
        c = t / an;
        s += c;
        if fabs(c) <= tol {
            break;
        }
    }

    if log_p {
        ans += log1p(a * s);
    } else {
        ans *= a * s + 1.;
    }
    ans
}

#[allow(clippy::too_many_arguments)]
/// Yields the incomplete beta ratio I_{1-x}(b,a) for a <= min(eps,eps*b),
/// b*x <= 1, and x <= 0.5, i.e., a is very small. Use only if above inequalities are satisfied.
pub fn apser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    let g: f64 = 0.577215664901533;

    let bx: f64 = b * x;

    let mut t: f64;
    t = x - bx;
    let c: f64 = if b * eps <= 0.02 {
        log(x) + psi(b) + g + t
    } else {
        // b > 2e13 : psi(b) ~= log(b)
        log(bx) + g + t
    };
    let tol: f64 = eps * 5.0 * fabs(c);
    let mut j: f64 = 1.0;
    let mut s: f64 = 0.0;
    let mut aj: f64;
    loop {
        j += 1.;
        t *= x - bx / j;
        aj = t / j;
        s += aj;
        if fabs(aj) <= tol {
            break;
        }
    }

    -a * (c + s)
}

#[allow(clippy::too_many_arguments)]
/// Evaluates I_x(a,b) - I_x(a+n,b) where n is a positive int.
///
/// eps is the tolerance used.
fn bup(a: f64, b: f64, x: f64, y: f64, n: i32, eps: f64, give_log: bool) -> f64 {
    let apb: f64 = a + b;
    let ap1: f64 = a + 1.0;
    let mut mu: i32;
    let mut k: i32;
    let mut d: f64;
    if n > 1 && a >= 1.0 && apb >= ap1 * 1.10 {
        mu = fabs(exparg(1)) as i32;
        k = exparg(0) as i32;
        if mu > k {
            mu = k;
        }
        d = exp(-mu as f64);
    } else {
        mu = 0;
        d = 1.0;
    }

    /* L10: */
    let mut ret_val = if give_log {
        brcmp1(mu, a, b, x, y, true) - log(a)
    } else {
        brcmp1(mu, a, b, x, y, false) / a
    };
    if n == 1 || (give_log && ret_val == ML_NEGINF) || (!give_log && ret_val == 0.) {
        return ret_val;
    }

    let nm1: i32 = n - 1;
    let mut w = d;
    let mut l: f64;

    // Let k be the index of the maximum term.

    k = 0;
    if b > 1.0 {
        if y > 1e-4 {
            let r = (b - 1.) * x / y - a;
            if r >= 1.0 {
                k = if r < nm1 as f64 { r as i32 } else { nm1 };
            }
        } else {
            k = nm1;
        }

        // Add the increasing terms of the series - if k > 0.
        /* L30: */
        for i in 0..k {
            l = i as f64;
            d *= (apb + l) / (ap1 + l) * x;
            w += d;
        }
    }

    // L40: Add the remaining terms of the series.

    for i in k..nm1 {
        l = i as f64;
        d *= (apb + l) / (ap1 + l) * x;
        w += d;
        if d <= eps * w {
            /* relativ convergence (eps) */
            break;
        }
    }

    // L50: Terminate the procedure.
    if give_log {
        ret_val += log(w);
    } else {
        ret_val *= w;
    }
    ret_val
}

#[allow(clippy::too_many_arguments)]
/// Continued fraction expansion for I_x(a,b) when a, b > 1.
///
/// It is assumed that  lambda = (a + b)*y - b.
fn bfrac(a: f64, b: f64, x: f64, y: f64, lambda: f64, eps: f64, log_p: bool) -> f64 {
    let mut e: f64;
    let mut n: f64;
    let mut p: f64;
    let mut r: f64;
    let mut s: f64;
    let mut t: f64;
    let mut w: f64;
    let mut r0: f64;
    let mut an: f64;
    let mut bn: f64;
    let mut anp1: f64;
    let mut bnp1: f64;
    let mut beta: f64;
    let mut alpha: f64;

    if !lambda.is_finite() {
        return f64::NAN;
    }

    let brc = brcomp(a, b, x, y, log_p);

    if brc.is_nan() {
        // e.g. from   L <- 1e308; pnbinom(L, L, mu = 5)
        return f64::NAN;
    }
    if !log_p && brc == 0.0 {
        return 0.0;
    }

    let c = lambda + 1.0;
    let c0 = b / a;
    let c1 = 1.0 / a + 1.0;
    let yp1 = y + 1.0;

    n = 0.0;
    p = 1.0;
    s = a + 1.0;
    an = 0.0;
    bn = 1.0;
    anp1 = 1.0;
    bnp1 = c / c1;
    r = c1 / c;

    // Continued fraction calculation.

    loop {
        n += 1.0;
        t = n / a;
        w = n * (b - n) * x;
        e = a / s;
        alpha = p * (p + c0) * e * e * (w * x);
        e = (t + 1.0) / (c1 + t + t);
        beta = n + w / s + e * (c + n * yp1);
        p = t + 1.0;
        s += 2.0;

        /* update an, bn, anp1, and bnp1 */

        t = alpha * an + beta * anp1;
        an = anp1;
        anp1 = t;
        t = alpha * bn + beta * bnp1;
        bn = bnp1;
        bnp1 = t;

        r0 = r;
        r = anp1 / bnp1;

        if fabs(r - r0) <= eps * r {
            break;
        }

        // Rescale an, bn, anp1, and bnp1.

        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0;
        if n >= 10000.0 {
            break;
        }
    }

    if log_p {
        brc + log(r)
    } else {
        brc * r
    }
}

#[allow(clippy::too_many_arguments)]
/// Evaluates x^a * y^b / Beta(a,b)
fn brcomp(a: f64, b: f64, x: f64, y: f64, log_p: bool) -> f64 {
    let const__ = 0.398942280401433; /* == 1/sqrt(2*pi); */
    let n: i32;
    let mut c: f64;
    let mut e: f64;
    let mut u: f64;
    let v: f64;
    let mut z: f64;
    let mut b0: f64;
    let apb: f64;

    if x == 0.0 || y == 0.0 {
        return r_d_0(log_p);
    }
    let a0 = min(a, b);
    if a0 < 8.0 {
        let lnx: f64;
        let lny: f64;
        #[allow(clippy::collapsible_else_if)]
        if x <= 0.375 {
            lnx = log(x);
            lny = alnrel(-x);
        } else {
            if y > 0.375 {
                lnx = log(x);
                lny = log(y);
            } else {
                lnx = alnrel(-y);
                lny = log(y);
            }
        }

        z = a * lnx + b * lny;
        if a0 >= 1.0 {
            z -= betaln(a, b);
            return r_d_exp(log_p, z);
        }

        // Procedure for a < 1 OR b < 1.

        b0 = max(a, b);
        if b0 >= 8.0 {
            /* L80: */
            u = gamln1(a0) + algdiv(a0, b0);

            return if log_p {
                log(a0) + (z - u)
            } else {
                a0 * exp(z - u)
            };
        }
        /* else : */

        if b0 <= 1.0 {
            /* algorithm for max(a,b) = b0 <= 1 */

            let e_z = r_d_exp(log_p, z);

            if !log_p && e_z == 0.0 {
                /* exp() underflow */
                return 0.0;
            }

            apb = a + b;
            if apb > 1.0 {
                u = a + b - 1.;
                z = (gam1(u) + 1.) / apb;
            } else {
                z = gam1(apb) + 1.;
            }

            c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;
            /* FIXME? log(a0*c)= log(a0)+ log(c) and that is improvable */
            return if log_p {
                e_z + log(a0 * c) - log1p(a0 / b0)
            } else {
                e_z * (a0 * c) / (a0 / b0 + 1.0)
            };
        }

        // else: algorithm for 1 < b0 < 8.
        u = gamln1(a0);
        n = (b0 - 1.0) as i32;
        if n >= 1 {
            c = 1.0;
            for _i in 1..=n {
                b0 += -1.0;
                c *= b0 / (a0 + b0);
            }
            u += log(c);
        }
        z -= u;
        b0 += -1.0;
        apb = a0 + b0;
        let t = if apb > 1.0 {
            u = a0 + b0 - 1.;
            (gam1(u) + 1.0) / apb
        } else {
            gam1(apb) + 1.0
        };

        if log_p {
            log(a0) + z + log1p(gam1(b0)) - log(t)
        } else {
            a0 * exp(z) * (gam1(b0) + 1.0) / t
        }
    } else {
        // Procedure for a >= 8 and b >= 8.
        let h: f64;
        let x0: f64;
        let y0: f64;
        let lambda: f64;
        if a <= b {
            h = a / b;
            x0 = h / (h + 1.0);
            y0 = 1. / (h + 1.0);
            lambda = a - (a + b) * x;
        } else {
            h = b / a;
            x0 = 1.0 / (h + 1.0);
            y0 = h / (h + 1.0);
            lambda = (a + b) * y - b;
        }

        e = -lambda / a;
        if fabs(e) > 0.6 {
            u = e - log(x / x0);
        } else {
            u = rlog1(e);
        }

        e = lambda / b;
        if fabs(e) <= 0.6 {
            v = rlog1(e);
        } else {
            v = e - log(y / y0);
        }

        z = if log_p {
            -(a * u + b * v)
        } else {
            exp(-(a * u + b * v))
        };

        if log_p {
            -M_LN_SQRT_2PI + 0.5 * log(b * x0) + z - bcorr(a, b)
        } else {
            const__ * sqrt(b * x0) * z * exp(-bcorr(a, b))
        }
    }
}

#[allow(clippy::too_many_arguments)]
// called only once from  bup(), as r = brcmp1(mu, a, b, x, y, false) / a;
fn brcmp1(mu: i32, a: f64, b: f64, x: f64, y: f64, give_log: bool) -> f64 {
    let const__ = 0.398942280401433; /* == 1/sqrt(2*pi); */
    let mut c: f64;
    let t: f64;
    let mut u: f64;
    let v: f64;
    let mut z: f64;
    let apb: f64;

    let a0 = min(a, b);
    if a0 < 8. {
        let lnx: f64;
        let lny: f64;
        if x <= 0.375 {
            lnx = log(x);
            lny = alnrel(-x);
        } else if y > 0.375 {
            // L11:
            lnx = log(x);
            lny = log(y);
        } else {
            lnx = alnrel(-y);
            lny = log(y);
        }

        // L20:
        z = a * lnx + b * lny;
        if a0 >= 1.0 {
            z -= betaln(a, b);
            return esum(mu, z, give_log);
        }
        // else :
        /* -----------------------------------------------------------------------
         */
        /*              PROCEDURE FOR A < 1 OR B < 1 */
        /* -----------------------------------------------------------------------
         */
        // L30:
        let mut b0 = max(a, b);
        if b0 >= 8.0 {
            /* L80:                  ALGORITHM FOR b0 >= 8 */
            u = gamln1(a0) + algdiv(a0, b0);
            // R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1, b0 >= 8;  z=%.15g\n", z);
            return if give_log {
                log(a0) + esum(mu, z - u, true)
            } else {
                a0 * esum(mu, z - u, false)
            };
        } else if b0 <= 1.0 {
            //                   a0 < 1, b0 <= 1
            let ans = esum(mu, z, give_log);
            if ans == (if give_log { ML_NEGINF } else { 0.0 }) {
                return ans;
            }

            apb = a + b;
            if apb > 1.0 {
                // L40:
                u = a + b - 1.0;
                z = (gam1(u) + 1.0) / apb;
            } else {
                z = gam1(apb) + 1.0;
            }
            // L50:
            c = if give_log {
                log1p(gam1(a)) + log1p(gam1(b)) - log(z)
            } else {
                (gam1(a) + 1.) * (gam1(b) + 1.) / z
            };
            // R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1, b0 <= 1;  c=%.15g\n", c);
            return if give_log {
                ans + log(a0) + c - log1p(a0 / b0)
            } else {
                ans * (a0 * c) / (a0 / b0 + 1.0)
            };
        }
        // else:               algorithm for	a0 < 1 < b0 < 8
        // L60:
        u = gamln1(a0);
        let n: i32 = (b0 - 1.0) as i32;
        if n >= 1 {
            c = 1.0;
            for _i in 1..=n {
                b0 += -1.;
                c *= b0 / (a0 + b0);
                /* L61: */
            }
            u += log(c); // TODO?: log(c) = log( prod(...) ) =  sum( log(...) )
        }
        // L70:
        z -= u;
        b0 += -1.;
        apb = a0 + b0;
        if apb > 1.0 {
            // L71:
            t = (gam1(apb - 1.) + 1.) / apb;
        } else {
            t = gam1(apb) + 1.;
        }
        // R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1 < b0 < 8;  t=%.15g\n", t);
        // L72:
        if give_log {
            // TODO? log(t) = log1p(..)
            log(a0) + esum(mu, z, true) + log1p(gam1(b0)) - log(t)
        } else {
            a0 * esum(mu, z, false) * (gam1(b0) + 1.) / t
        }
    } else {
        /* -----------------------------------------------------------------------
         */
        /*              PROCEDURE FOR A >= 8 AND B >= 8 */
        /* -----------------------------------------------------------------------
         */
        // L100:
        let h: f64;
        let x0: f64;
        let y0: f64;
        let lambda: f64;
        if a > b {
            // L101:
            h = b / a;
            x0 = 1. / (h + 1.); // => lx0 := log(x0) = 0 - log1p(h)
            y0 = h / (h + 1.);
            lambda = (a + b) * y - b;
        } else {
            h = a / b;
            x0 = h / (h + 1.); // => lx0 := log(x0) = - log1p(1/h)
            y0 = 1. / (h + 1.);
            lambda = a - (a + b) * x;
        }
        let lx0 = -log1p(b / a); // in both cases

        // R_ifDEBUG_printf(
        //    " brcmp1(mu,a,b,*): a,b >= 8;	x0=%.15g, lx0=log(x0)=%.15g\n", x0,
        //    lx0);
        // L110:
        let mut e = -lambda / a;
        if fabs(e) > 0.6 {
            // L111:
            u = e - log(x / x0);
        } else {
            u = rlog1(e);
        }

        // L120:
        e = lambda / b;
        if fabs(e) > 0.6 {
            // L121:
            v = e - log(y / y0);
        } else {
            v = rlog1(e);
        }

        // L130:
        z = esum(mu, -(a * u + b * v), give_log);
        if give_log {
            log(const__) + (log(b) + lx0) / 2. + z - bcorr(a, b)
        } else {
            const__ * sqrt(b * x0) * z * exp(-bcorr(a, b))
        }
    }
}

#[allow(clippy::too_many_arguments)]
/// Asymptotic Expansion for I_x(a,b)  when a is larger than b.
///
/// Computes w := w + I_x(a,b)
/// It is assumed a >= 15 and b <= 1.
/// `eps` is the tolerance used.
/// `ierr` is a variable that reports the status of the results.
///
/// if(log_w),  *w  itself must be in log-space;
/// compute w := w + I_x(a,b)  but return *w = log(w):
/// ////// *w := log(exp(*w) + I_x(a,b)) = logspace_add(*w, log( I_x(a,b) ))
/// ----------------------------------------------------------------------- */
fn bgrat(a: f64, b: f64, x: f64, y: f64, w: &mut f64, eps: f64, ierr: &mut i32, log_w: bool) {
    const N_TERMS_BGRAT: usize = 30;
    let mut c: [f64; N_TERMS_BGRAT] = [0.0; N_TERMS_BGRAT];
    let mut d: [f64; N_TERMS_BGRAT] = [0.0; N_TERMS_BGRAT];
    let bm1: f64 = b - 0.5 - 0.5;
    let nu: f64 = a + bm1 * 0.5; /* nu = a + (b-1)/2 =: T, in (9.1) of
                                  * Didonato & Morris(1992), p.362 */
    let lnx: f64 = if y > 0.375 { log(x) } else { alnrel(-y) };
    let z: f64 = -nu * lnx; // z =: u in (9.1) of D.&M.(1992)

    if b * z == 0.0 {
        // should not happen, but does, e.g.,
        // for  pbeta(1e-320, 1e-5, 0.5)  i.e., _subnormal_ x,
        // Warning ... bgrat(a=20.5, b=1e-05, x=1, y=9.99989e-321): ..
        //printf("bgrat(a=%g, b=%g, x=%g, y=%g): z=%g, b*z == 0 underflow, hence "
        //       "inaccurate pbeta()",
        //       a, b, x, y, z);
        /* L_Error:    THE EXPANSION CANNOT BE COMPUTED */
        *ierr = 1;
        return;
    }

    /*                 COMPUTATION OF THE EXPANSION */
    /* r1 = b * (gam1(b) + 1.) * exp(b * log(z)),// = b/gamma(b+1) z^b = z^b /
     * gamma(b) set r := exp(-z) * z^b / gamma(b) ; gam1(b) = 1/gamma(b+1) - 1
     * , b in [-1/2, 3/2] */
    // exp(a*lnx) underflows for large (a * lnx); e.g. large a ==> using
    // log_r := log(r): r = r1 * exp(a * lnx) * exp(bm1 * 0.5 * lnx);
    // log(r)=log(b) + log1p(gam1(b)) + b * log(z) + (a * lnx) + (bm1 *
    // 0.5 * lnx),
    let log_r = log(b) + log1p(gam1(b)) + b * log(z) + nu * lnx;
    // FIXME work with  log_u = log(u)  also when log_p=false  (??)
    // u is 'factored out' from the expansion {and multiplied back, at
    // the end}:
    // algdiv(b,a) = log(gamma(a)/gamma(a+b))
    let log_u = log_r - (algdiv(b, a) + b * log(nu));
    /* u = (log_p) ? log_r - u : exp(log_r-u); // =: M  in (9.2) of {reference
    above} */
    /* u = algdiv(b, a) + b * log(nu);// algdiv(b,a) =
    log(gamma(a)/gamma(a+b)) */
    // u = (log_p) ? log_u : exp(log_u); // =: M  in (9.2) of {reference
    // above}
    let u = exp(log_u);

    if log_u == ML_NEGINF {
        // R_ifDEBUG_printf(
        //    " bgrat(*): underflow log_u = -Inf  = log_r -u', log_r = %g ", log_r);
        /* L_Error:    THE EXPANSION CANNOT BE COMPUTED */
        *ierr = 2;
        return;
    }

    let u_0 = u == 0.0; // underflow --> do work with log(u) == log_u !
    let l: f64 = // := *w/u .. but with care: such that it also works when u
             // underflows to 0:
     if log_w { if *w == ML_NEGINF { 0.0 } else { exp(*w - log_u) }
     } else if *w == 0.0 { 0.0 } else { exp(log(*w) - log_u) };

    // R_ifDEBUG_printf(" bgrat(a=%g, b=%g, x=%g, *)\n -> u=%g, l='w/u'=%g, ", a, b,
    //                  x, u, l);
    let q_r = grat_r(b, z, log_r, eps); // = q/r of former grat1(b,z, r, &p, &q)
    let v = 0.25 / (nu * nu);
    let t2 = lnx * 0.25 * lnx;
    let mut j = q_r;
    let mut sum = j;
    let mut t = 1.0;
    let mut cn = 1.0;
    let mut n2 = 0.0;
    for n in 0..=N_TERMS_BGRAT {
        let bp2n = b + n2;
        j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
        n2 += 2.0;
        t *= t2;
        cn /= n2 * (n2 + 1.);
        let nm1 = n - 1;
        c[nm1] = cn;
        let mut s = 0.0;
        if n > 1 {
            let mut coef = b - (n as f64);
            for i in 1..=nm1 {
                s += coef * c[i - 1] * d[nm1 - i];
                coef += b;
            }
        }
        d[nm1] = bm1 * cn + s / (n as f64);
        let dj = d[nm1] * j;
        sum += dj;
        if sum <= 0.0 {
            // R_ifDEBUG_printf(" bgrat(*): sum_n(..) <= 0; should not happen (n=%d)\n",
            //                 n);
            /* L_Error:    THE EXPANSION CANNOT BE COMPUTED */
            *ierr = 3;
            return;
        }
        if fabs(dj) <= eps * (sum + l) {
            *ierr = 0;
            break;
        } else if n == N_TERMS_BGRAT {
            // never? ; please notify R-core if seen:
            *ierr = 4;
            // printf("bgrat(a=%g, b=%g, x=%g) *no* convergence: NOTIFY R-core!\n "
            //       "dj=%g, rel.err=%g\n",
            //       a, b, x, dj, fabs(dj) / (sum + l));
        }
    } // for(n .. n_terms..)

    /*                    ADD THE RESULTS TO W */

    if log_w {
        // *w is in log space already:
        *w = logspace_add(*w, log_u + log(sum));
    } else {
        *w += if u_0 { exp(log_u + log(sum)) } else { u * sum };
    }
}

#[allow(clippy::too_many_arguments)]
/// Scaled complement of incomplete gamma ratio function
/// grat_r(a,x,r) :=  Q(a,x) / r
///  where
/// Q(a,x) = pgamma(x,a, lower.tail=false)
///  and r = e^(-x)* x^a / Gamma(a) ==  exp(log_r)
///
/// It is assumed that a <= 1.  eps is the tolerance to be used.
fn grat_r(a: f64, x: f64, log_r: f64, eps: f64) -> f64 {
    // Called only from bgrat() as q_r = grat_r(b, z, log_r, eps) :
    if a * x == 0.0 {
        /* L130: */
        if x <= a {
            /* L100: */
            exp(-log_r)
        } else {
            /* L110:*/
            0.0
        }
    } else if a == 0.5 {
        // e.g. when called from pt()
        /* L120: */
        if x < 0.25 {
            let p: f64 = erf_(sqrt(x));
            // R_ifDEBUG_printf(" grat_r(a=%g, x=%g ..)): a=1/2 --> p=erf__(.)= %g\n", a,
            //                 x, p);
            return (0.5 - p + 0.5) * exp(-log_r);
        } else {
            // 2013-02-27: improvement for "large" x: direct computation of
            // q/r:
            let sx: f64 = sqrt(x);
            let q_r = erfc1(1, sx) / sx * M_SQRT_PI;
            // R_ifDEBUG_printf(
            //    " grat_r(a=%g, x=%g ..)): a=1/2 --> q_r=erfc1(..)/r= %g\n", a, x,
            //    q_r);
            return q_r;
        }
    } else if x < 1.1 {
        /* L10:  Taylor series for  P(a,x)/x^a */

        let mut an = 3.0;
        let mut c = x;
        let mut sum = x / (a + 3.0);
        let tol = eps * 0.1 / (a + 1.0);
        let mut t;
        loop {
            an += 1.;
            c *= -(x / an);
            t = c / (a + an);
            sum += t;
            if fabs(t) <= tol {
                break;
            }
        }

        // R_ifDEBUG_printf(
        //    " grat_r(a=%g, x=%g, log_r=%g): sum=%g; Taylor w/ %.0f terms", a, x,
        //    log_r, sum, an - 3.);
        let j = a * x * ((sum / 6. - 0.5 / (a + 2.)) * x + 1. / (a + 1.));
        let z = a * log(x);
        let h = gam1(a);
        let g = h + 1.0;

        if (x >= 0.25 && (a < x / 2.59)) || (z > -0.13394) {
            // L40:
            let l = rexpm1(z);
            let q = ((l + 0.5 + 0.5) * j - l) * g - h;
            if q <= 0.0 {
                // R_ifDEBUG_printf(" => q_r= 0.\n");
                /* L110:*/
                return 0.0;
            } else {
                // R_ifDEBUG_printf(" => q_r=%.15g\n", q * exp(-log_r));
                return q * exp(-log_r);
            }
        } else {
            let p = exp(z) * g * (0.5 - j + 0.5);
            // R_ifDEBUG_printf(" => q_r=%.15g\n", (0.5 - p + 0.5) * exp(-log_r));
            return /* q/r = */ (0.5 - p + 0.5) * exp(-log_r);
        }
    } else {
        /* L50: ----  (x >= 1.1)  ---- Continued Fraction Expansion */

        let mut a2n_1 = 1.0;
        let mut a2n = 1.0;
        let mut b2n_1 = x;
        let mut b2n = x + (1. - a);
        let mut c = 1.0;
        let mut am0;
        let mut an0;

        loop {
            a2n_1 = x * a2n + c * a2n_1;
            b2n_1 = x * b2n + c * b2n_1;
            am0 = a2n_1 / b2n_1;
            c += 1.0;
            let c_a = c - a;
            a2n = a2n_1 + c_a * a2n;
            b2n = b2n_1 + c_a * b2n;
            an0 = a2n / b2n;
            if fabs(an0 - am0) < eps * an0 {
                break;
            }
        }

        // R_ifDEBUG_printf(
        //    " grat_r(a=%g, x=%g, log_r=%g): Cont.frac. %.0f terms => q_r=%.15g\n",
        //    a, x, log_r, c - 1., an0);
        /* q/r = (r * an0)/r = */
        an0
    }
}

#[allow(clippy::too_many_arguments)]
/// Asymptotic expansion for I_x(a,b) for large a and b.
///
/// Lambda = (a + b)*y - b and eps is the tolerance used.
/// It is assumed that lambda is nonnegative and that a and b are
/// greater than or equal to 15.
fn basym(a: f64, b: f64, lambda: f64, eps: f64, log_p: bool) -> f64 {
    /* ------------------------ */
    /*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
    /*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
    const NUM_IT: usize = 20;
    /*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

    #[allow(clippy::approx_constant)]
    const E0: f64 = 1.12837916709551; /* e0 == 2/sqrt(pi) */
    const E1: f64 = 0.353553390593274; /* e1 == 2^(-3/2)   */
    const LN_E0: f64 = 0.120782237635245; /* == ln(e0) */

    let mut a0: [f64; NUM_IT + 1] = [0.0; NUM_IT + 1];
    let mut b0: [f64; NUM_IT + 1] = [0.0; NUM_IT + 1];
    let mut c: [f64; NUM_IT + 1] = [0.0; NUM_IT + 1];
    let mut d: [f64; NUM_IT + 1] = [0.0; NUM_IT + 1];

    let f = a * rlog1(-lambda / a) + b * rlog1(lambda / b);
    let t;
    if log_p {
        t = -f;
    } else {
        t = exp(-f);
        if t == 0.0 {
            return 0.0; /* once underflow, always underflow .. */
        }
    }
    let z0 = sqrt(f);
    let z = z0 / E1 * 0.5;
    let z2 = f + f;
    let h: f64;
    let r0: f64;
    let r1: f64;
    let w0: f64;

    if a < b {
        h = a / b;
        r0 = 1.0 / (h + 1.);
        r1 = (b - a) / b;
        w0 = 1.0 / sqrt(a * (h + 1.));
    } else {
        h = b / a;
        r0 = 1. / (h + 1.);
        r1 = (b - a) / a;
        w0 = 1. / sqrt(b * (h + 1.));
    }

    a0[0] = r1 * 0.666_666_666_666_666_6;
    c[0] = a0[0] * -0.5;
    d[0] = -c[0];
    let mut j0 = 0.5 / E0 * erfc1(1, z0);
    let mut j1 = E1;
    let mut sum = j0 + d[0] * w0 * j1;

    let mut s = 1.0;
    let h2 = h * h;
    let mut hn = 1.0;
    let mut w = w0;
    let mut znm1 = z;
    let mut zn = z2;
    for n in (2..=NUM_IT).step_by(2) {
        hn *= h2;
        a0[n - 1] = r0 * 2. * (h * hn + 1.) / (n as f64 + 2.);
        let np1 = n + 1;
        s += hn;
        a0[np1 - 1] = r1 * 2. * s / (n as f64 + 3.0);

        for i in n..=np1 {
            let r = (i as f64 + 1.0) * -0.5;
            b0[0] = r * a0[0];
            for m in 2..=i {
                let mut bsum = 0.0;
                for j in 1..=(m - 1) {
                    let mmj = m - j;
                    bsum += (j as f64 * r - mmj as f64) * a0[j - 1] * b0[mmj - 1];
                }
                b0[m - 1] = r * a0[m - 1] + bsum / m as f64;
            }
            c[i - 1] = b0[i - 1] / (i as f64 + 1.0);

            let mut dsum = 0.0;
            for j in 1..=i {
                dsum += d[i - j - 1] * c[j - 1];
            }
            d[i - 1] = -(dsum + c[i - 1]);
        }

        j0 = E1 * znm1 + (n as f64 - 1.0) * j0;
        j1 = E1 * zn + n as f64 * j1;
        znm1 *= z2;
        zn *= z2;
        w *= w0;
        let t0 = d[n - 1] * w * j0;
        w *= w0;
        let t1 = d[np1 - 1] * w * j1;
        sum += t0 + t1;
        if fabs(t0) + fabs(t1) <= eps * sum {
            break;
        }
    }

    if log_p {
        LN_E0 + t - bcorr(a, b) + log(sum)
    } else {
        let u = exp(-bcorr(a, b));
        E0 * t * u * sum
    }
}

#[allow(clippy::too_many_arguments)]
/// If l = 0 then exparg(l) = The largest possible w for which exp(w) can be computed.
///
/// With 0.99999 fuzz ==> exparg(0) = 709.7756 nowadays.
///
/// if l = 1 (nonzero) then  exparg(l) = the largest negative W for
/// which the computed value of exp(W) is nonzero.
/// With 0.99999 fuzz ==> exparg(1) = -709.0825
/// nowadays
///
/// Note... only an approximate value for exparg(L) is needed.
fn exparg(l: i32) -> f64 {
    #[allow(clippy::approx_constant)]
    const LNB: f64 = 0.69314718055995;
    let m = if l == 0 { i1mach(16) } else { i1mach(15) - 1 };
    (m as f64) * LNB * 0.99999
}

#[allow(clippy::too_many_arguments)]
/// Evaluates exp(mu + x).
fn esum(mu: i32, x: f64, give_log: bool) -> f64 {
    if give_log {
        return x + mu as f64;
    }

    // else :
    let w: f64;
    if x > 0.0 {
        /* L10: */
        if mu > 0 {
            return exp(mu as f64) * exp(x);
        }
        w = (mu as f64) + x;
        if w < 0.0 {
            return exp(mu as f64) * exp(x);
        }
    } else {
        /* x <= 0 */
        if mu < 0 {
            return exp(mu as f64) * exp(x);
        }
        w = (mu as f64) + x;
        if w > 0.0 {
            return exp(mu as f64) * exp(x);
        }
    }
    exp(w)
}

/// Evaluates the function exp(x) - 1.
fn rexpm1(x: f64) -> f64 {
    let p1: f64 = 9.14041914819518e-10;
    let p2: f64 = 0.0238082361044469;
    let q1: f64 = -0.499999999085958;
    let q2: f64 = 0.107141568980644;
    let q3: f64 = -0.0119041179760821;
    let q4: f64 = 5.95130811860248e-4;

    if fabs(x) <= 0.15 {
        x * (((p2 * x + p1) * x + 1.) / ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.))
    } else {
        let w = exp(x);
        if x > 0.0 {
            w * (0.5 - 1.0 / w + 0.5)
        } else {
            w - 0.5 - 0.5
        }
    }
}

/// Calculates the function ln(1 + a).
fn alnrel(a: f64) -> f64 {
    if fabs(a) > 0.375 {
        return log(1.0 + a);
    }
    // else : |a| <= 0.375
    let p1 = -1.29418923021993;
    let p2 = 0.405303492862024;
    let p3 = -0.0178874546012214;
    let q1 = -1.62752256355323;
    let q2 = 0.747811014037616;
    let q3 = -0.0845104217945565;
    let t = a / (a + 2.);
    let t2 = t * t;
    let w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.0) / (((q3 * t2 + q2) * t2 + q1) * t2 + 1.0);
    t * 2. * w
}

/// Evaluates the function x - ln(1 + x).
fn rlog1(x: f64) -> f64 {
    let a = 0.0566749439387324;
    let b = 0.0456512608815524;
    let p0 = 0.333333333333333;
    let p1 = -0.224696413112536;
    let p2 = 0.00620886815375787;
    let q1 = -1.27408923933623;
    let q2 = 0.354508718369557;

    let mut h;
    let w;
    let w1;
    if x < -0.39 || x > 0.57 {
        /* direct evaluation */
        w = x + 0.5 + 0.5;
        return x - log(w);
    }
    /* else */
    if x < -0.18 {
        /* L10: */
        h = x + 0.3;
        h /= 0.7;
        w1 = a - h * 0.3;
    } else if x > 0.18 {
        /* L20: */
        h = x * 0.75 - 0.25;
        w1 = b + h / 3.0;
    } else {
        /*		Argument Reduction */
        h = x;
        w1 = 0.0;
    }

    /* L30:              	Series Expansion */

    let r = h / (h + 2.);
    let t = r * r;
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
    t * 2. * (1. / (1. - r) - r * w) + w1
}

/// Evaluates the real error function.
fn erf_(x: f64) -> f64 {
    const C: f64 = 0.564189583547756;
    const A: [f64; 5] = [
        7.7105849500132e-5,
        -0.00133733772997339,
        0.0323076579225834,
        0.0479137145607681,
        0.128379167095513,
    ];
    const B: [f64; 3] = [0.00301048631703895, 0.0538971687740286, 0.375795757275549];
    const P: [f64; 8] = [
        -1.36864857382717e-7,
        0.564195517478974,
        7.21175825088309,
        43.1622272220567,
        152.98928504694,
        339.320816734344,
        451.918953711873,
        300.459261020162,
    ];
    const Q: [f64; 8] = [
        1.,
        12.7827273196294,
        77.0001529352295,
        277.585444743988,
        638.980264465631,
        931.35409485061,
        790.950925327898,
        300.459260956983,
    ];
    const R: [f64; 5] = [
        2.10144126479064,
        26.2370141675169,
        21.3688200555087,
        4.6580782871847,
        0.282094791773523,
    ];
    const S: [f64; 4] = [
        94.153775055546,
        187.11481179959,
        99.0191814623914,
        18.0124575948747,
    ];

    let mut t: f64;
    
    let bot: f64;
    let top: f64;

    let ax = fabs(x);
    if ax <= 0.5 {
        t = x * x;
        top = (((A[0] * t + A[1]) * t + A[2]) * t + A[3]) * t + A[4] + 1.0;
        bot = ((B[0] * t + B[1]) * t + B[2]) * t + 1.0;

        return x * (top / bot);
    }

    // else:  |x| > 0.5

    if ax <= 4.0 {
        //  |x| in (0.5, 4]
        top = ((((((P[0] * ax + P[1]) * ax + P[2]) * ax + P[3]) * ax + P[4]) * ax + P[5]) * ax
            + P[6])
            * ax
            + P[7];
        bot = ((((((Q[0] * ax + Q[1]) * ax + Q[2]) * ax + Q[3]) * ax + Q[4]) * ax + Q[5]) * ax
            + Q[6])
            * ax
            + Q[7];
        let r = 0.5 - exp(-x * x) * top / bot + 0.5;
        return if x < 0.0 { -r } else { r };
    }

    // else:  |x| > 4

    if ax >= 5.8 {
        return if x > 0.0 { 1.0 } else { -1.0 };
    }

    // else:  4 < |x| < 5.8
    let x2: f64 = x * x;
    t = 1. / x2;
    top = (((R[0] * t + R[1]) * t + R[2]) * t + R[3]) * t + R[4];
    bot = (((S[0] * t + S[1]) * t + S[2]) * t + S[3]) * t + 1.0;
    t = (C - top / (x2 * bot)) / ax;
    let r = 0.5 - exp(-x2) * t + 0.5;
    if x < 0.0 {
        -r
    } else {
        r
    }
}

#[allow(unused_variables)]
/// Evaluates the complementary error function.
///
/// erfc1(ind,x) = erfc(x) if ind = 0,
/// erfc1(ind,x) = exp(x*x)*erfc(x) otherwise.
fn erfc1(ind: i32, x: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Computes 1/gamma(a+1) - 1 for -0.5 <= a <= 1.5.
fn gam1(a: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates ln(gamma(1 + a)) for -0.2 <= a <= 1.25.
fn gamln1(a: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates the Digamma function psi(x).
///
/// Psi(xx) is assigned the value 0 when the digamma function cannot
/// be computed.
///
/// The main computation involves evaluation of rational Chebyshev
/// approximations published in Math. Comp. 27, 123-127(1973) by
/// Cody, Strecok and Thacher.
///
/// Psi was written at Argonne National Laboratory for the FUNPACK
/// package of special function subroutines. Psi was modified by
/// A.H. Morris (NSWC).
fn psi(x: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates the logarithm of the beta function ln(beta(a0,b0)).
fn betaln(a0: f64, b0: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates the function ln(gamma(a + b)) for 1 <= a <= 2 and 1 <= b <= 2.
fn gsumln(a: f64, b: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates del(a0) + del(b0) - del(a0 + b0) where ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
///
/// It is assumed that a0 >= 8 and b0 >= 8.
fn bcorr(a0: f64, b0: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Computation of ln(gamma(b)/gamma(a+b)) when b >= 8.
fn algdiv(a: f64, b: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates ln(gamma(a)) for positive a.
///
/// Written by Alfred H. Morris
/// Naval Surface Warfare Center
/// Dahlgren, Virginia
fn gamln(a: f64) -> f64 {
    panic!("not implemented");
}
