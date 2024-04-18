//! A Rust port of R's C Library of Special Functions.

// Avoid converting `if n < 1 || n > 1000` to `if !(1..=1000).contains(&n)`.
#![allow(clippy::manual_range_contains)]

mod bd0;
mod chebyshev;
mod cospi;
mod d1mach;
mod dnorm;
mod dpois;
mod dpq;
mod gamma;
mod i1mach;
mod lbeta;
mod lgamma;
mod lgammacor;
mod libc;
mod nmath;
mod pgamma;
mod pnorm;
mod qnorm;
mod rmath;
mod stirlerr;
mod toms708;

// Use only explicit exports and no wildcard exports to avoid accidentally
// exporting symbols that should not be exported.
pub use bd0::bd0;
pub use bd0::ebd0;
pub use chebyshev::chebyshev_eval;
pub use chebyshev::chebyshev_init;
pub use cospi::cospi;
pub use cospi::sinpi;
pub use cospi::tanpi;
pub use dpois::dpois;
pub use gamma::gammafn;
pub use i1mach::i1mach;
pub use lbeta::lbeta;
pub use lgamma::lgammafn;
pub use lgamma::lgammafn_sign;
pub use lgammacor::lgammacor;
pub use pgamma::log1pmx;
pub use pgamma::logspace_add;
pub use pgamma::pgamma;
pub use rmath::dnorm;
pub use rmath::pnorm;
pub use rmath::qnorm;
pub use stirlerr::stirlerr;
