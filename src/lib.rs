//! Shared scientific vocabulary and reusable quantum simulation components.
//!
//! The facade always provides [`core`]. Optional Cargo features expose exact,
//! variational, stochastic series expansion, and artifact capabilities as they
//! become supported. Scientific behavior is introduced test-first under the
//! conventions documented by the project.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// Foundational scientific types shared by every qslib algorithm.
pub use qslib_core as core;

#[cfg(feature = "exact")]
/// Exact-basis and exact-solver capabilities.
pub use qslib_exact as exact;

#[cfg(feature = "io")]
/// Versioned scientific configuration and artifact capabilities.
pub use qslib_io as io;

#[cfg(feature = "sse")]
/// Finite-temperature stochastic series expansion capabilities.
pub use qslib_sse as sse;

#[cfg(feature = "variational")]
/// Variational statistics, TDVP, and integration capabilities.
pub use qslib_variational as variational;
