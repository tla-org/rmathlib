//! 
//!  Mathlib : A C Library of Special Functions
//!  Copyright (C) 2006-2019 The R Core Team
//!  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
//!  Copyright (C) 2005-10 The R Foundation
//! 
//!  This program is free software; you can redistribute it and/or modify
//!  it under the terms of the GNU General Public License as published by
//!  the Free Software Foundation; either version 2 of the License, or
//!  (at your option) any later version.
//! 
//!  This program is distributed in the hope that it will be useful,
//!  but WITHOUT ANY WARRANTY; without even the implied warranty of
//!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//!  GNU General Public License for more details.
//! 
//!  You should have received a copy of the GNU General Public License
//!  along with this program; if not, a copy is available at
//!  https://www.R-project.org/Licenses/
//! 
//!  SYNOPSIS
//! 
//!  #include <Rmath.h>
//! 
//!  double pgamma (double x, double alph, double scale,
//!  	       int lower_tail, int log_p)
//! 
//!  double log1pmx	(double x)
//!  double lgamma1p (double a)
//! 
//!  double logspace_add (double logx, double logy)
//!  double logspace_sub (double logx, double logy)
//!  double logspace_sum (double* logx, int n)
//! 
//! 
//!  DESCRIPTION
//! 
//!  This function computes the distribution function for the
//!  gamma distribution with shape parameter alph and scale parameter
//!  scale.	This is also known as the incomplete gamma function.
//!  See Abramowitz and Stegun (6.5.1) for example.
//! 
//!  NOTES
//! 
//!  Complete redesign by Morten Welinder, originally for Gnumeric.
//!  Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler

use libm::log1p;

/// Continued fraction for calculation of 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
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

        let scalefactor: f64 = f64::powf(2.0, 256.0);

        if b2.abs() > scalefactor {
            a1 /= scalefactor;
            b1 /= scalefactor;
            a2 /= scalefactor;
            b2 /= scalefactor;
        } else if b2.abs() < 1.0 / scalefactor {
            a1 *= scalefactor;
            b1 *= scalefactor;
            a2 *= scalefactor;
            b2 *= scalefactor;
        }
    }
    a2 / b2
}

/// Accurately calculates log(1+x)-x, particularly for small x.
pub fn log1pmx(x: f64) -> f64 {
    let min_log_1_value = -0.79149064;

    if x > 1.0 || x < min_log_1_value {
        return log1p(x) - x;
    } else {
        // * -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
	    //  log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
	    //  ---------------------------------------------
	    //  S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
        let r: f64 = x / (2.0 + x);
        let y: f64 = r * r;
        if x.abs() < 1e-2 {
            let two: f64 = 2.0;
            return r * ((((two / 9.0 * y + two / 7.0) * y + two / 5.0) * y + two / 3.0) * y - x);
        } else {
            let tol_logcf = 1e-14;
            return r * (2.0 * y * logcf(y, 3.0, 2.0, tol_logcf) - x);
        }
    }
}
