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
- API stability: `cargo semver-checks check-release -p qslib-quantum
  --baseline-rev 2584261 --release-type patch --only-explicit-features`
  (cargo-semver-checks 0.42.0) compared the current facade with baseline
  commit `2584261` and passed 165 checks with 12 skips; the broader disposition
  is recorded in
  [`docs/api-stability.md`](api-stability.md).
- Packaging: the ABI3 wheel and Maturin source distribution both install in a
  temporary environment; each runs the ten Python contract tests and the exact
  four-site example. The current bundle is
  `/private/tmp/qslib-release-final-MmZAR0`, built from commit `61c5138`.
  The guarded release workflow is the reproducible path for rebuilding a
  bundle from the current revision.
- Distribution preparation (2026-07-20): Python packaging metadata now carries
  the physicist-facing description, Apache-2.0 license, repository and
  documentation URLs, classifiers, and a crate-local README/license boundary
  that works in both wheels and source distributions. The release workflow has
  a guarded PyPI trusted-publishing job and builds Linux, macOS, and Windows
  CLI archives with the same archive contract tests.
- Revalidation (2026-07-20 09:58Z): the full Rust 1.85 workspace test matrix,
  stable Clippy and rustdoc with warnings denied, formatting, all facade
  feature boundaries, conformance/workspace harnesses, Markdown links, and CI
  YAML parsing passed on the clean source revision used for the current bundle.
- Hosted CI (2026-07-20): GitHub Actions run
  [`29752717024`](https://github.com/lere01/qslib/actions/runs/29752717024) for
  commit `e63466f` completed successfully. All 17 jobs passed, including the
  Linux/macOS/Windows Rust stable and MSRV matrix, Python 3.12/3.13 wheels on
  all three operating systems, Miri, branch coverage, bounded fuzzing, facade
  feature checks, formatting, documentation, and dependency policy.
- Hosted revalidation (2026-07-20): GitHub Actions run
  [`29755394438`](https://github.com/lere01/qslib/actions/runs/29755394438) for
  commit `61c5138` completed successfully with the same 17-job matrix. This
  run includes the corrected guarded release workflow in the pushed source.
- ncli parity (2026-07-20 10:20Z): the separately owned parent repository has
  an explicit optional qslib backend with TFIM, signed J1-J2, Rydberg, and exact
  spectrum parity tests. Its native backend remains the default. The adapter
  preflights a 256 MiB dense budget, rejects nonzero Rydberg diagonals and
  non-Hermitian outputs, and the parent workflow installs qslib at pinned
  revision `e63466f` and marks the parity tests required on all three hosted
  operating systems. Local execution against the current ABI3 wheel passed
  all 11 parity tests with `NCLI_REQUIRE_QSLIB_PARITY=1`.
- Publication preparation (2026-07-20 10:20Z): `.github/workflows/release.yml`
  defines a guarded artifact build and a manual, tag-checked GitHub release job.
  The build clean-installs and exercises wheel/sdist artifacts, smoke-tests the
  CLI, and verifies the generated checksum manifest before upload.
- Release rehearsal (2026-07-20): the local workflow-equivalent build exposed
  and corrected an invalid `--locked` flag on `maturin sdist`. Wheel and sdist
  installs, ten Python contract tests, the exact ground-state example, CLI
  help, mdBook, the combined Rust API site, and relocated checksum verification
  all passed for commit `61c5138`.
- Checksums: `SHA256SUMS` is generated from inside the bundle with `./...`
  relative paths and verifies both in place and after relocation.

## Known external gates

The hosted qslib Linux/macOS/Windows matrix is green for commit `e63466f` as
recorded above. Its Rust jobs explicitly exclude the Python cdylib, which is
built and tested through Maturin jobs on all three platforms. The local semver
comparison above does not replace a Linux registry or release-baseline check.
The ncli adapter remains a separate ownership unit and must be validated in
its own repository's CI.
The Python cdylib is packaged through Maturin; a workspace-wide Cargo release
link is not a supported way to build that extension on macOS.

## Migration status

qslib uses canonical row-major site order, little-endian site-zero bit packing,
explicit simulation bases, and resolved pair coefficients. Legacy ncli and
standalone-SSE inputs require the named adapters documented in the migration
chapters. No destructive migration is performed by the release bundle.
