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
- Coverage: `cargo llvm-cov --locked --workspace --all-features --summary-only`
  passed with 78.28% line coverage and 71.63% region coverage. The tool emitted
  no branch counters for this run; Python FFI and the CLI binary entry point
  were not exercised by this Rust-only coverage command.
- API stability: `cargo-semver-checks 0.42.0` compared the current facade with
  baseline commit `2584261` as an assumed patch release and passed 165 checks
  with 12 skips.

## Known external gates

The remote Linux/macOS/Windows CI matrix and nightly Miri execution are authored
but not executed in this local-only workflow. The local semver comparison above
does not replace a Linux registry or release-baseline check.
The ncli backend-selection/parity adapter remains a separate ownership unit.
The Python cdylib is packaged through Maturin; a workspace-wide Cargo release
link is not a supported way to build that extension on macOS.

## Migration status

qslib uses canonical row-major site order, little-endian site-zero bit packing,
explicit simulation bases, and resolved pair coefficients. Legacy ncli and
standalone-SSE inputs require the named adapters documented in the migration
chapters. No destructive migration is performed by the release bundle.
