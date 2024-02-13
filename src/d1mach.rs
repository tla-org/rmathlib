//! Mathlib - A Mathematical Function Library
//! Copyright (C) 1998  Ross Ihaka
//! Copyright (C) 2000-2014 The R Core Team
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

use crate::libc::DBL_EPSILON;
use crate::libc::DBL_MAX;
use crate::libc::DBL_MIN;
use crate::rmath::M_LOG10_2;

pub fn d1mach(i: i32) -> f64 {
    match i {
        1 => DBL_MIN,
        2 => DBL_MAX,
        3 => 0.5 * DBL_EPSILON,
        4 => DBL_EPSILON,
        5 => M_LOG10_2,
        _ => 0.0,
    }
}
