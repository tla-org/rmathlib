//! Based on C translation of ACM TOMS 708
//!
//! These routines provide very high relative accuracy; about 14 digits.

use crate::libc::*;
use libm::log1p;
use crate::rmath::*;
use crate::d1mach::d1mach;
use crate::i1mach::i1mach;
use crate::dpq::*;

/// Calculates the function exp(x) - 1
fn rexpm1(x: f64) -> f64 {
    let p1: f64 = 9.14041914819518e-10;
    let p2: f64 = 0.0238082361044469;
    let q1: f64 = -0.499999999085958;
    let q2: f64 = 0.107141568980644;
    let q3: f64 = -0.0119041179760821;
    let q4: f64 = 5.95130811860248e-4;

    if x.abs() <= 0.15 {
        return x * (((p2 * x + p1) * x + 1.) /
            ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.));
    } else { /* |x| > 0.15 : */
        let w: f64 = exp(x);
        if x > 0.0 {
            return w * (0.5 - 1. / w + 0.5);
        } else {
            return w - 0.5 - 0.5;
        }
    }
}

// Computes ln(1 + 1).
fn alnrel(a: f64) -> f64 {
    if a.abs() > 0.375 {
        return log(1.0 + a);
    } else {
        let p1: f64 = -1.29418923021993;
        let p2: f64 = 0.405303492862024;
        let p3: f64 = -0.0178874546012214;
        let q1: f64 = -1.62752256355323;
        let q2: f64 = 0.747811014037616;
        let q3: f64 = -0.0845104217945565;
        let t = a / (a + 2.0);
        let t2 = t * t;
        let w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) /
            (((q3 * t2 + q2) * t2 + q1) * t2 + 1.);
        return t * 2. * w;
    }

}

fn r_log1_exp(x: f64) -> f64 {
    if x > -M_LN2 {
        log(-rexpm1(x))
    } else {
        log1p(-exp(x))
    }
}

fn exparg(l: i32) -> f64 {
    // --------------------------------------------------------------------
    // If l = 0 then  exparg(l) = The largest positive W for which
    // exp(W) can be computed. With 0.99999 fuzz  ==> exparg(0) = 709.7756 nowadays.

    // Iff l = 1 (nonzero) then  exparg(l) = the largest negative W for
    // which the computed value of exp(W) is nonzero.
    // With 0.99999 fuzz ==> exparg(1) = -709.0825 nowadays.

    // Note: only an approximate value for exparg(L) is needed.
    // -------------------------------------------------------------------- */
    let lnb: f64 = 0.69314718055995;
    let m: i32 = if l == 0 { i1mach(16) } else { i1mach(15) - 1 };
    m as f64 * lnb * 0.99999
}

/// Evaluates I (A, B) for `b < min(eps, eps * a)` and `x <= 0.5`.
fn fpser(a: f64, b: f64, x: f64, eps: f64, log_p: bool) -> f64 {
    let mut ans: f64;
    let mut c: f64;
    let mut s: f64;
    let mut t: f64;
    let mut an: f64;
    let mut tol: f64;

    if log_p {
        ans = a * log(x);
    } else if a > eps * 0.001 {
        t = a * log(x);
        if t < exparg(1) {
            // exp(t) would underflow.
            return 0.0;
        }
        ans = exp(t);
    } else {
        ans = 1.0;
    }

    return 0.0;
}

/// Power SERies expansion for evaluating I_x(a, b) when b <= 1 or b*x <= 0.7.
/// `eps` is the tolerance used.
/// Note: if `log_p` is true, then also use it if `(b < 40 & lambda > 650)`.
fn bpser(a: f64, b: f64, x: f64, eps: f64, log_p: bool) -> f64 {
    let mut i: i32;
    let mut m: i32;
    let mut ans: f64;
    let mut c: f64;
    let mut t: f64;
    let mut u: f64;
    let mut z: f64;
    let mut a0: f64;
    let mut b0: f64;
    let mut apb: f64;

    if x == 0.0 {
        return r_d__0(log_p);
    }

    // Compute the factor x^a / (a*Beta(a, b)).
    a0 = min(a, b);
    if a0 >= 1.0 {
        // TODO: Implement betaln.
        // z = a * log(x) - betaln(a, b);
    }

    // TODO
    return 0.0;
}

fn l_end(do_swap: bool, w: &mut f64, w1: &mut f64) {
    if do_swap {
        let tmp = *w;
        *w = *w1;
        *w1 = tmp;
    }
}

fn l_w_bpser(do_swap: bool, a0: f64, b0: f64, x0: f64, eps: f64, log_p: bool, w: &mut f64, w1: &mut f64) -> (f64, f64) {
    *w = bpser(a0, b0, x0, eps, log_p);
    *w1 = if log_p { r_log1_exp(*w) } else { 0.5 - *w + 0.5 };
    l_end(do_swap, w, w1);
    return (*w, *w1);
}


/// Calculates the Incomplete Beta Function I_x(a, b)
///
/// ALGORITHM 708, COLLECTED ALGORITHMS FROM ACM.
/// This work published in  Transactions On Mathematical Software,
/// vol. 18, no. 3, September 1992, pp. 360-373.
pub fn bratio(a: f64, b: f64, x: f64, y: f64, log_p: bool) -> (f64, f64, i32) {
    // -----------------------------------------------------------------------
    // 
    //        Evaluation of the Incomplete Beta function I_x(a,b)
    // 
    //  	       --------------------
    // 
    //     It is assumed that a and b are nonnegative, and that x <= 1
    //     and y = 1 - x.  Bratio assigns w and w1 the values
    // 
    //  		w  = I_x(a,b)
    //  		w1 = 1 - I_x(a,b)
    // 
    //     ierr is a variable that reports the status of the results.
    //     If no input errors are detected then ierr is set to 0 and
    //     w and w1 are computed. otherwise, if an error is detected,
    //     then w and w1 are assigned the value 0 and ierr is set to
    //     one of the following values ...
    // 
    //    ierr = 1  if a or b is negative
    //    ierr = 2  if a = b = 0
    //    ierr = 3  if x < 0 or x > 1
    //    ierr = 4  if y < 0 or y > 1
    //    ierr = 5  if x + y != 1
    //    ierr = 6  if x = a = 0
    //    ierr = 7  if y = b = 0
    //    ierr = 8	(not used currently)
    //    ierr = 9  NaN in a, b, x, or y
    //    ierr = 10     (not used currently)
    //    ierr = 11  bgrat() error code 1 [+ warning in bgrat()]
    //    ierr = 12  bgrat() error code 2   (no warning here)
    //    ierr = 13  bgrat() error code 3   (no warning here)
    //    ierr = 14  bgrat() error code 4 [+ WARNING in bgrat()]
    // 
    // 
    // --------------------
    //     Written by Alfred H. Morris, Jr.
    //    Naval Surface Warfare Center
    //    Dahlgren, Virginia
    //     Revised ... Nov 1991
    // ----------------------------------------------------------------------

    let do_swap: bool;
    let n: i32 = 0;
    let mut ierr: i32 = 0;
    let z: f64;
    let a0: f64;
    let b0: f64;
    let x0: f64;
    let y0: f64;
    let lambda: f64;

    // eps is a machine dependent constant: the smallest
    // floating point number for which 1.0 + eps > 1.0
    // NOTE: for almost all purposes it is replaced by 1e-15 (~= 4.5 times larger) below */
    let mut eps: f64 = 2.0 * d1mach(3); // == DBL_EPSILON (in R, Rmath)

    let mut w: f64 = r_d__0(log_p);
    let mut w1: f64 = r_d__0(log_p);

    if x.is_nan() || y.is_nan() || a.is_nan() || b.is_nan() {
        ierr = 9;
        return (w, w1, ierr);
    }
    if a < 0.0 || b < 0.0 {
        ierr = 1;
        return (w, w1, ierr);
    }
    if a == 0.0 && b == 0.0 {
        ierr = 2;
        return (w, w1, ierr);
    }
    if x < 0.0 || x > 1.0 {
        ierr = 3;
        return (w, w1, ierr);
    }
    if y < 0.0 || y > 1.0 {
        ierr = 4;
        return (w, w1, ierr);
    }

    // Check that `y == 1 - x`.
    z = x + y - 0.5 - 0.5;

    if z.abs() > eps * 3.0 {
        ierr = 5;
        return (w, w1, ierr);
    }

    ierr = 0;
    if x == 0.0 {
        if a == 0.0 {
            ierr = 6;
            return (w, w1, ierr);
        }
    }
    if y == 0.0 {
        if b == 0.0 {
            ierr = 7;
            return (w, w1, ierr);
        }
    }
    if a == 0.0 {
        w = r_d__1(log_p);
        w1 = r_d__0(log_p);
        return (w, w1, ierr);
    }
    if b == 0.0 {
        w = r_d__0(log_p);
        w1 = r_d__1(log_p);
        return (w, w1, ierr);
    }

    eps = max(eps, 1e-15);

    let a_lt_b: bool = a < b;
    if (/* max(a, b) */ if a_lt_b { b } else { a }) < eps * 0.001 {
        // Procedure for a and b < 0.001 * eps.
        // w = a / (a + b) and w1 = b / (a + b).
        if log_p {
            if a_lt_b {
                w = log1p(-a / (a + b)); // Notably if a << b.
                w1 = log(a / (a + b));
            } else {
                w = log(b / (a + b));
                w1 = log1p(-b / (a + b));
            }
        } else {
            w = b / (a + b);
            w1 = a / (a + b);
        }
        return (w, w1, ierr);
    }

    if min(a, b) <= 1.0 {
        // a <= 1 or b <= 1.
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
        // Now have `x0 <= 1/2 < y0` (still `x0 + y0 == 1`).

        if b0 < min(eps, eps * a0) {
            w = fpser(a0, b0, x0, eps, log_p);
            w1 = if log_p { r_log1_exp(w) } else { 0.5 - w + 0.5 };
            l_end(do_swap, &mut w, &mut w1);
            return (w, w1, ierr);
        }

        let did_bup: bool = false;
        if max(a0, b0) > 1.0 {
            if b0 <= 1.0 {
                // TODO: implement bpser.
            }
        }

    }

    return (0.0, 0.0, 0);
}

/// Evaluates ln(gamma(1 + a)) for -0.2 <= a <= 1.25.
fn gamln1(a: f64) -> f64 {
    let w: f64;
    if a < 0.6 {
        let p0: f64 = 0.577215664901533;
        let p1: f64 = 0.844203922187225;
        let p2: f64 = -0.168860593646662;
        let p3: f64 = -0.780427615533591;
        let p4: f64 = -0.402055799310489;
        let p5: f64 = -0.0673562214325671;
        let p6: f64 = -0.00271935708322958;
        let q1: f64 = 2.88743195473681;
        let q2: f64 = 3.12755088914843;
        let q3: f64 = 1.56875193295039;
        let q4: f64 = 0.361951990101499;
        let q5: f64 = 0.0325038868253937;
        let q6: f64 = 6.67465618796164e-4;
        w = ((((((p6 * a + p5)* a + p4)* a + p3)* a + p2)* a + p1)* a + p0) /
            ((((((q6 * a + q5)* a + q4)* a + q3)* a + q2)* a + q1)* a + 1.0);
        return -(a) * w;
    } else { // 0.6 <= a <= 1.25.
        let r0: f64 = 0.422784335098467;
        let r1: f64 = 0.848044614534529;
        let r2: f64 = 0.565221050691933;
        let r3: f64 = 0.156513060486551;
        let r4: f64 = 0.017050248402265;
        let r5: f64 = 4.97958207639485e-4;
        let s1: f64 = 1.24313399877507;
        let s2: f64 = 0.548042109832463;
        let s3: f64 = 0.10155218743983;
        let s4: f64 = 0.00713309612391;
        let s5: f64 = 1.16165475989616e-4;
        let x: f64 = a - 0.5 - 0.5;
        w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
            (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.0);
        return x * w;
    }
}

/// Computes ln(gamma(b) / gamma(a + b)) when b >= 8.
fn algdiv(a: f64, b: f64) -> f64 {
    // In this algorithm, del(x) is the function defined
    // by ln(gamma(x)) = (x - 0.5) * ln(x) - x + 0.5 * ln(2 * pi) + del(x).
    let c0: f64 = 0.0833333333333333;
    let c1: f64 = -0.00277777777760991;
    let c2: f64 = 7.9365066682539e-4;
    let c3: f64 = -5.9520293135187e-4;
    let c4: f64 = 8.37308034031215e-4;
    let c5: f64 = -0.00165322962780713;

    let mut c: f64;
    let mut d: f64;
    let mut h: f64;
    let mut t: f64;
    let mut u: f64;
    let mut v: f64;
    let mut w: f64;
    let mut x: f64;
    let mut s3: f64;
    let mut s5: f64;
    let mut x2: f64;
    let mut s7: f64;
    let mut s9: f64;
    let mut s11: f64;

    if a > b {
        h = b / a;
        c = 1.0 / (h + 1.0);
        x = h / (h + 1.0);
        d = a + (b - 0.5);
    } else {
        h = a / b;
        c = h / (h + 1.0);
        x = 1.0 / (h + 1.0);
        d = b + (a - 0.5);
    }

    // Set s<n> = (1 - x^n) / (1 - x).
    
    x2 = x * x;
    s3 = x + x2 + 1.0;
    s5 = x + x2 * s3 + 1.0;
    s7 = x + x2 * s5 + 1.0;
    s9 = x + x2 * s7 + 1.0;
    s11 = x + x2 * s9 + 1.0;

    // w := del(b) - del(a + b).
    
    t = 1.0 / (b*b);
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
        s3) * t + c0;
    w *= c / b;

    // Combine the results.
    
    u = d * alnrel(a / b);
    v = a * (log(b) - 1.);
    if u > v {
        return w - v - u;
    } else {
        return w - u - v;
    }
}

/// Evaluates ln(gamma(a)) for positive a.
///
/// Written by Alfred H. Morris.
/// Naval Surface Warfare Center.
/// Dahlgren, Virginia.
fn gamln(a: f64) -> f64 {
    let d: f64 = 0.418938533204673; // d == 0.5*(LN(2*PI) - 1)
    let c0: f64 = 0.0833333333333333;
    let c1: f64 = -0.002277777777760991;
    let c2: f64 = 7.9365066682539e-4;
    let c3: f64 = -5.9520293135187e-4;
    let c4: f64 = 8.37308034031215e-4;
    let c5: f64 = -0.00165322962780713;

    if a <= 0.8 {
        return gamln1(a) - log(a);
    } else if a <= 2.25 {
        return gamln1(a - 0.5 - 0.5);
    } else if a < 10.0 {
        let n: i32 = (a - 1.25) as i32;
        let mut t = a;
        let mut w = 1.0;
        for _ in 1..=n {
            t += -1.0;
            w *= t;
        }
        return gamln1(t - 1.0) + log(w);
    } else {
        // a >= 10.
        let t = 1.0 / (a * a);
        let w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
        return d + w + (a - 0.5) * (log(a) - 1.0);
    }
}













