# Local release-candidate evidence

This document describes the local evidence bundle for the current `0.1.0`
implementation snapshot. It is not a publication or a release authorization.

## Test matrix

- Rust 1.85.0: locked workspace tests, formatting, and feature checks.
- Current stable Rust: workspace tests, Clippy with warnings denied, rustdoc
  with warnings denied, and feature-boundary checks.
- Python: ABI3 wheel contract tests on the locally available interpreter, plus
  the exact four-site example.
- Documentation: mdBook, Rust API docs, doctests, Markdown links, and combined
  site API-copy validation.
- Policy: cargo-deny advisories, licenses, and sources.
- Hardening: bounded structured CLI parser/resolution fuzz smoke and expanded
  kernel benchmark target.

## Known external gates

The remote Linux/macOS/Windows CI matrix, nightly Miri execution, and Linux
semver-checks run are authored but not executed in this local-only workflow.
The ncli backend-selection/parity adapter remains a separate ownership unit.
The Python cdylib is packaged through Maturin; a workspace-wide Cargo release
link is not a supported way to build that extension on macOS.

## Migration status

qslib uses canonical row-major site order, little-endian site-zero bit packing,
explicit simulation bases, and resolved pair coefficients. Legacy ncli and
standalone-SSE inputs require the named adapters documented in the migration
chapters. No destructive migration is performed by the release bundle.
