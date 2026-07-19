//! Independent conformance fixtures and test oracles for qslib.
//!
//! This crate is never a production dependency. It validates language-neutral
//! fixture envelopes and their provenance without calling qslib scientific
//! implementations.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod fixture;

pub use fixture::{
    Authorship, CONVENTION_SCHEMA, FIXTURE_SCHEMA, Fixture, FixtureError, FixtureKind, Oracle,
    load_conformance_fixtures, required_fixture_kinds, validate_fixture_set,
};
