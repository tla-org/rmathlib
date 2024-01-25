//! A Rust port of R's C Library of Special Functions.
mod chebyshev;
mod cospi;
mod dnorm;
mod dpq;
mod gamma;
mod lgamma;
mod lgammacor;
mod libc;
mod nmath;
mod pgamma;
mod pnorm;
mod qnorm;
mod rmath;

pub use chebyshev::*;
pub use cospi::*;
pub use dnorm::*;
pub use gamma::*;
pub use lgamma::*;
pub use lgammacor::*;
pub use libc::*;
pub use nmath::*;
pub use pgamma::*;
pub use pnorm::*;
pub use qnorm::*;
