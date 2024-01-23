#![allow(non_snake_case)]

pub fn r_d__0(log_p: bool) -> f64 {
  if log_p { f64::NEG_INFINITY } else { 0.0 }
}

pub fn r_d__1(log_p: bool) -> f64 {
    if log_p { 0.0 } else { 1.0 }
}

pub fn r_dt_0(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail { r_d__0(log_p) } else { r_d__1(log_p) }
}

pub fn r_dt_1(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail { r_d__1(log_p) } else { r_d__0(log_p) }
}
