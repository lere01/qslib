# ADR-0006: Public stability and feature policy

- Status: accepted
- Date: 2026-07-19
- Owners: qslib maintainers

## Context

The facade must be useful without pulling every backend into every consumer.
At 1.0, users need a clear distinction between stable scientific behavior,
backend adapters, and experiments whose API may still evolve.

## Decision

Core is always present in the root facade. Stable optional features are
`exact`, `variational`, `sse`, and `io`; `full` enables all stable Rust library
features. Python and CLI remain separate packages. Backend interop uses
separately named, non-default adapter features.

Features are additive and must not change the meaning of an existing type or
calculation. Heavy optional dependencies appear only behind the capability that
uses them. Stable facade APIs follow semantic versioning from 1.0. Experimental
items live under an `experimental` module or explicitly named feature and are
not re-exported into stable modules.

All public items state their stability. No experimental algorithm is exported
from the 1.0 facade merely because it is feature-gated. Public enums that may
grow across minor
versions are non-exhaustive where that does not weaken exhaustive scientific
validation. Serialized schemas have their own compatibility rules and never
derive stability merely from Rust type layout.

## Alternatives considered

- Enabling every capability by default was rejected because it makes a basic
  geometry consumer compile solvers, Parquet, SSE, and interfaces.
- A facade with no feature gates was rejected because downstream applications
  would need to know every internal crate.
- Treating all pre-1.0 work as permanently unstable was rejected because the
  release plan requires a finite, dependable 1.0 surface.

## Consequences

Feature combinations add CI cost, but dependency weight and stability promises
become visible. Internal crates may evolve more quickly when the facade contract
is preserved.

## Validation

- CI checks no-default, default, each stable capability, and full builds.
- `cargo tree` verifies that core-only excludes numerical, Arrow, CLI, Python,
  and SSE dependencies.
- Semver checks compare the facade against the previous release candidate.
- Documentation labels experimental modules and generates feature-specific
  examples correctly.
