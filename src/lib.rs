//! A Rust port of R's C Library of Special Functions.
mod cospi;
mod dpq;
mod nmath;
mod pnorm;
mod rmath;
mod qnorm;
mod libc;

pub use cospi::*;
pub use pnorm::*;
