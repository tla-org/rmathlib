//!
//! AUTHORS
//! Catherine Loader, catherine@research.bell-labs.com, October 23, 2000. [ bd0() ]
//! Morten Welinder, see Bugzilla PR#15628, 2014                          [ebd0() ]
//!
//! Merge in to R (and much more):
//!
//! Copyright (C) 2000-2022 The R Core Team
//!
//! This program is free software; you can redistribute it and/or modify
//! it under the terms of the GNU General Public License as published by
//! the Free Software Foundation; either version 2 of the License, or
//! (at your option) any later version.
//!
//! This program is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//! GNU General Public License for more details.
//!
//! You should have received a copy of the GNU General Public License
//! along with this program; if not, a copy is available at
//! https://www.R-project.org/Licenses/

use crate::nmath::*;
use crate::libc::*;

/// Calculates a stable deviance part using a method that reduces relative error
///
/// Evaluates the "deviance part"
/// bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
///           =  x * log(x/M) + M - x
/// where M = E[X] = n*p (or = lambda), for	  x, M > 0
///
/// in a manner that should be stable (with small relative error)
/// for all x and M=np. In particular for x/np close to 1, direct
/// evaluation fails, and evaluation is based on the Taylor series
/// of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).
///
/// Martyn Plummer had the nice idea to use log1p() and Martin Maechler
/// emphasized the extra need to control cancellation.
///
/// MP:   t := (x-M)/M  ( <==> 1+t = x/M  ==>
///
/// bd0 = M*[ x/M * log(x/M) + 1 - (x/M) ] = M*[ (1+t)*log1p(t) + 1 - (1+t) ]
///     = M*[ (1+t)*log1p(t) - t ] =: M * p1log1pm(t) =: M * p1l1(t)
/// MM: The above is very nice, as the "simple" p1l1() function would be useful
///   to have available in a fast numerical stable way more generally.
pub fn bd0(x: f64, np: f64) -> f64 {
    if !r_finite(x) || !r_finite(np) || np == 0.0 {
        return ml_warn_return_nan();
    }

    if fabs(x-np) < 0.1*(x+np) {
        let mut v: f64 = (x-np)/(x+np);
        let mut s: f64 = (x-np)*v;
        if fabs(s) < DBL_MIN {
            return s;
        }
        let mut ej: f64 = 2.0*x*v;
        v *= v; // v^2
        for j in 1..1000 {
            // Taylor series; 1000: no infinite loop
            // as |v| < 0.1, v^2000 is "zero".
            ej *= v;
            let s_: f64 = s;
            s += ej/((j<<1)+1) as f64;
            if s == s_ {
                // Last term was effectively 0.
                return s;
            }
        }
    }
    println!("bd0: T.series failed to converge in 1000 iterations");
    return x*log(x/np)+np-x;
}
