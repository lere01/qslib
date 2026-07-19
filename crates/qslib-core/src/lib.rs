//! Foundational scientific vocabulary shared by qslib algorithms.
//!
//! Checked identifiers, basis meaning, geometry, interactions, operators,
//! models, symmetry primitives, and observable definitions are part of the
//! supported qslib 1.0 API as each milestone completes its conformance gate.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod basis;
mod error;
mod identifiers;
mod scalar;

pub use basis::{
    BasisBit, BasisState, BasisStateView, FullBasis, PackedState, SectorBasis, WordWidth,
};
pub use error::BasisError;
pub use identifiers::{PhysicalAxis, SimulationBasis, SiteCount, SiteId};
pub use scalar::{Complex64, Real, ensure_finite};
