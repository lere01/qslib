# ADR-0009: Registry package identifiers for the qslib brand

- Status: accepted
- Date: 2026-07-19
- Owners: qslib owner

## Context

The `qslib` name is already an actively published crates.io package for the
unrelated QuantStudio qPCR library, currently version 0.15.2. The same name is
also present on PyPI. crates.io names cannot be reused by a different project.
Local package names and import targets can differ from registry distribution
names, but that distinction must be deliberate and documented.

## Decision

Keep qslib as the project brand and Rust library import name. The owner selected
`qslib-quantum` as the Rust facade and Python distribution name on 2026-07-19.
Read-only crates.io searches on that date returned no matching package. Backend
packages use the `qslib-quantum-` distribution prefix where they are published.
The facade declares `[lib] name = "qslib"`, so Rust source writes
`use qslib::...`.

The Python distribution uses `qslib-quantum` and the collision-safe import
package is `qslib_quantum`. Project branding does not override package-manager
or interpreter namespace safety.

Every workspace package remains `publish = false` until the owner separately
authorizes external publication. Naming acceptance is not publication
authority and does not reserve a registry name.

## Alternatives considered

- Reusing `qslib` on crates.io is impossible without transfer from the existing
  owner and would remain confusing.
- Renaming the entire project avoids package/import distinction but discards the
  owner's chosen brand.
- Keeping `import qslib` preserves the shortest spelling but conflicts with the
  existing qPCR package in a shared Python environment.
- Declaring that qslib will never use public registries removes the collision,
  but changes the intended distribution route and requires an explicit owner
  decision.

## Consequences

Documentation must distinguish project, distribution, package, and import
names. The selected unique distribution improves discoverability and prevents
dependency confusion. Capability packages must use the same prefix so a future
facade can depend on published registry packages rather than path-only crates.

## Validation

- Read-only registry searches find no existing package under the selected
  identifiers immediately before acceptance.
- `cargo metadata`, built wheel metadata, Rust imports, and Python imports use
  the documented names.
- Package smoke tests install from local artifacts into clean environments.
- Publication remains disabled until a separate owner authorization.
