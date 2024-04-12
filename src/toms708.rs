#![allow(dead_code)]
#![allow(clippy::manual_range_contains)]

use libm::exp;
use libm::expm1;
use libm::fabs;
use libm::log;
use libm::log1p;
use libm::pow;

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

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
/// Yields the incomplete beta ratio I_{1-x}(b,a) for a <= min(eps,eps*b),
/// b*x <= 1, and x <= 0.5, i.e., a is very small. Use only if above inequalities are satisfied.
pub fn apser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
/// Evaluates I_x(a,b) - I_x(a+n,b) where n is a positive int.
///
/// eps is the tolerance used.
fn bup(a: f64, b: f64, x: f64, y: f64, n: i32, eps: f64, give_log: bool) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
/// Continued fraction expansion for I_x(a,b) when a, b > 1.
///
/// It is assumed that  lambda = (a + b)*y - b.
fn bfrac(a: f64, b: f64, x: f64, y: f64, lambda: f64, eps: f64, log_p: bool) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
/// Evaluates x^a * y^b / Beta(a,b)
fn brcomp(a: f64, b: f64, x: f64, y: f64, log_p: bool) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
// called only once from  bup(), as r = brcmp1(mu, a, b, x, y, false) / a;
fn brcmp1(mu: i32, a: f64, b: f64, x: f64, y: f64, give_log: bool) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
fn bgrat(a: f64, b: f64, x: f64, y: f64, w: &mut f64, eps: f64, ierr: &mut i32, log_w: bool) {
    panic!("not implemented");
}

#[allow(unused_variables)]
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
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
/// Asymptotic expansion for I_x(a,b) for large a and b.
///
/// Lambda = (a + b)*y - b and eps is the tolerance used.
/// It is assumed that lambda is nonnegative and that a and b are
/// greater than or equal to 15.
fn basym(a: f64, b: f64, lambda: f64, eps: f64, log_p: bool) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
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
    panic!("not implemented");
}

#[allow(unused_variables)]
#[allow(clippy::too_many_arguments)]
/// Evaluates exp(mu + x).
fn esum(mu: i32, x: f64, give_log: i32) -> f64 {
    panic!("not implemented");
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

#[allow(unused_variables)]
/// Calculates the function ln(1 + a).
fn alnrel(x: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates the function x - ln(1 + x).
fn rlog1(x: f64) -> f64 {
    panic!("not implemented");
}

#[allow(unused_variables)]
/// Evaluates the real error function.
fn erf_(x: f64) -> f64 {
    panic!("not implemented");
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
