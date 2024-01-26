//! A Rust port of R's C Library of Special Functions.

// Avoid converting `if n < 1 || n > 1000` to `if !(1..=1000).contains(&n)`.
#![allow(clippy::manual_range_contains)]

mod chebyshev;
mod cospi;
mod dnorm;
// mod dpois;
mod dpq;
mod gamma;
mod lgamma;
mod lgammacor;
mod libc;
mod nmath;
// mod pgamma;
mod pnorm;
mod qnorm;
mod rmath;
mod stirlerr;

pub use chebyshev::*;
pub use cospi::*;
// pub use dpois::*;
pub use gamma::*;
pub use lgamma::*;
pub use lgammacor::*;
pub use libc::*;
pub use nmath::*;
pub use rmath::dnorm;
pub use rmath::pnorm;
pub use rmath::qnorm;
// pub use pgamma::*;
pub use stirlerr::*;
