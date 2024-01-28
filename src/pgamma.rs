use libm::log1p;

const SCALEFACTOR: f64 = 1.157921e+77;

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
    todo!()
}

/// Continued fraction for calculation of 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
///
/// `eps` denotes the relative tolerance.
fn logcf(x: f64, i: f64, d: f64, eps: f64) -> f64 {
    let mut c1: f64 = 2.0 * d;
    let mut c2: f64 = i + d;
    let mut c4: f64 = c2 + d;
    let mut a1: f64 = c2;
    let mut b1: f64 = i * (c2 - i * x);
    let mut b2: f64 = d * d * x;
    let mut a2: f64 = c4 * c2 - b2;

    b2 = c4 * b1 - i * b2;

    while (a2 * b1 - a1 * b2).abs() > (eps * b1 * b2).abs() {
        let mut c3: f64 = c2 * c2 * x;
        c2 += d;
        c4 += d;
        a1 = c4 * a2 - c3 * a1;
        b1 = c4 * b2 - c3 * b1;

        c3 = c1 * c1 * x;
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

/// Accurately calculates log(1+x)-x, particularly for small x.
pub fn log1pmx(x: f64) -> f64 {
    let min_log_1_value = -0.79149064;

    if x > 1.0 || x < min_log_1_value {
        log1p(x) - x
    } else {
        // * -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
        //  log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
        //  ---------------------------------------------
        //  S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
        let r: f64 = x / (2.0 + x);
        let y: f64 = r * r;
        if x.abs() < 1e-2 {
            let two: f64 = 2.0;
            r * ((((two / 9.0 * y + two / 7.0) * y + two / 5.0) * y + two / 3.0) * y - x)
        } else {
            let tol_logcf = 1e-14;
            r * (2.0 * y * logcf(y, 3.0, 2.0, tol_logcf) - x)
        }
    }
}
