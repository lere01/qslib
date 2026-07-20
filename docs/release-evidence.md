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
- Generated invariants: bounded packed-state byte/word round trips and
  permutation inverse/compose properties pass on Rust 1.85 and under Miri.
- Fuzzing: the isolated `fuzz/state_conversion` libFuzzer target passes
  `cargo fuzz run state_conversion --sanitizer none -- -runs=1000` with no
  crash; CI has the same bounded invocation.
- Safety: nightly Rust 1.99.0 Miri passes every qslib-core target, including
  the generated conversion tests.
- Dependency audit: `cargo audit` exits zero with one allowed unmaintained
  `paste 1.0.15` warning in the Parquet dependency tree; cargo-deny reports
  the same documented exception.
- Coverage: nightly `cargo llvm-cov --locked --workspace --all-features
  --branch --summary-only` passed with 77.99% line coverage, 77.66% region
  coverage, and 57.76% branch coverage. The file-level report was reviewed;
  the CLI binary entry point and Python FFI remain unexercised by Rust tests,
  while core scientific branches have direct conformance or invalid-input
  coverage.
- API stability: `cargo-semver-checks 0.42.0` compared the current facade with
  baseline commit `2584261` as an assumed patch release and passed 165 checks
  with 12 skips; the broader disposition is recorded in
  [`docs/api-stability.md`](api-stability.md).
- Packaging: the ABI3 wheel and Maturin source distribution both install in a
  temporary environment; each runs the ten Python contract tests and the exact
  four-site example. The current bundle is
  `/private/tmp/qslib-release-candidate-20260720k`.
- Checksums: `SHA256SUMS` is generated from inside the bundle with `./...`
  relative paths and verifies both in place and after relocation.

## Known external gates

The remote Linux/macOS/Windows CI matrix is authored but not executed in this
local-only workflow. The local semver comparison above does not replace a
Linux registry or release-baseline check. The ncli backend-selection/parity
adapter remains a separate ownership unit.
The Python cdylib is packaged through Maturin; a workspace-wide Cargo release
link is not a supported way to build that extension on macOS.

## Migration status

qslib uses canonical row-major site order, little-endian site-zero bit packing,
explicit simulation bases, and resolved pair coefficients. Legacy ncli and
standalone-SSE inputs require the named adapters documented in the migration
chapters. No destructive migration is performed by the release bundle.
