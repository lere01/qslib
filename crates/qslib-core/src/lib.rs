//! Foundational scientific vocabulary shared by qslib algorithms.
//!
//! Checked identifiers, basis meaning, geometry, interactions, operators,
//! models, symmetry primitives, and observable definitions are part of the
//! supported qslib 1.0 API as each milestone completes its conformance gate.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod basis;
mod error;
mod geometry;
mod identifiers;
mod interactions;
mod models;
mod operators;
mod randomness;
mod scalar;
mod symmetry;

pub use basis::{
    BasisBit, BasisState, BasisStateView, FullBasis, PackedState, SectorBasis, WordWidth,
};
pub use error::BasisError;
pub use geometry::{
    Bond, BondMultiplicity, Boundary, Coordinate, CustomGeometry, GeometryError, LatticeKind,
    RectangularGeometry, ShellTolerance, XMajorAdapter,
};
pub use identifiers::{PhysicalAxis, SimulationBasis, SiteCount, SiteId};
pub use interactions::{
    DenseCouplings, DisorderProvenance, DisorderRealization, InteractionChannel, InteractionError,
    InteractionIdentity, InteractionTable, SparseCouplings, WeightedInteraction,
};
pub use models::{
    ModelError, ModelSpecification, ResolvedModel, heisenberg, j1j2, j1j2_disordered, rydberg, tfim,
};
pub use operators::{Hamiltonian, OperatorError, Pauli, PauliString};
pub use randomness::{QSLIB_SEED_SCHEME, derive_seed, expand_master_seed};
pub use scalar::{Complex64, Real, ensure_finite};
pub use symmetry::{
    DiagonalGauge, FiniteGroup, Permutation, SpinInversion, SymmetryCharacter, SymmetryError,
    is_interaction_symmetry, legacy_x_major_permutation, rectangle_point_group, square_point_group,
    sublattice_gauge, translation, translation_group, validate_model_symmetry,
    validate_spin_inversion,
};
