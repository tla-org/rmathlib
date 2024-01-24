//! A Rust port of R's C Library of Special Functions.
mod cospi;
mod dpq;
mod gamma;
mod lgamma;
mod lgammacor;
mod nmath;
mod pgamma;
mod pnorm;
mod qnorm;
mod rmath;

pub use cospi::*;
pub use gamma::*;
pub use lgamma::*;
pub use lgammacor::*;
pub use pgamma::*;
pub use pnorm::*;
pub use qnorm::*;
