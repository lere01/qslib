# Rust toolchain and MSRV policy

qslib uses Rust edition 2024. The minimum supported Rust version for qslib 1.0
is Rust 1.85.0, the first stable release supporting edition 2024. Every stable
workspace package declares `rust-version = "1.85"` through workspace metadata.

The checked-in `rust-toolchain.toml` selects exactly Rust 1.85.0 and includes
rustfmt and Clippy, so ordinary local commands continuously exercise the MSRV.
CI must test both this pinned MSRV and current stable. A dependency that cannot
compile on the MSRV is not accepted merely because it works on the development
compiler.

Raising the MSRV requires an accepted ADR that records the dependency or
language requirement, user impact, and migration date. Within the 1.x series,
the project should provide at least one minor-release notice before an MSRV
increase unless a security fix makes that impractical.

The lockfile is committed for reproducible workspace, CLI, Python, and release
candidate validation. Library dependency requirements remain semver ranges, but
CI validates the locked graph and a fresh dependency resolution separately.

The standard checks are listed in `AGENTS.md` and `CONTRIBUTING.md`. Nightly
Rust may be used by optional analysis tools such as fuzzers, but supported
library code and release artifacts must not require nightly features.
