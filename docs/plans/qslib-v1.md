# Deliver qslib 1.0 as a reusable quantum simulation library

This ExecPlan is a living document governed by `PLANS.md` at the qslib project
root. The sections `Progress`, `Surprises and discoveries`, `Decision log`, and
`Outcomes and retrospective` must be maintained as implementation proceeds.

## Purpose and user-visible outcome

qslib 1.0 will provide one coherent Rust library for the reusable scientific
parts of the existing ncli and SSE programs. A physicist will be able to define
a lattice and pair-dependent couplings, construct a TFIM, isotropic Heisenberg,
J1-J2, disordered exchange, or Rydberg Hamiltonian, use symmetries and conserved
sectors, compute exact ground states and small-system dynamics, evaluate common
observables, solve reusable TDVP numerical problems, and run supported
finite-temperature SSE simulations. The same conventions will be available
through Rust, a documented command line, and Python bindings suitable for ncli.

A complete local demonstration will construct a four-site model with
heterogeneous couplings, print its exact ground-state energy and observables,
evolve a state for a short time, run a small thermal SSE job where the model is
supported, and reproduce the same exact result through Python. The release
candidate will build from a clean checkout and pass all Rust, Python,
documentation, conformance, portability, and packaging gates.

This plan defines completion as qslib 1.0. Quantum simulation is an open-ended
field, so capabilities listed under non-goals are not allowed to delay 1.0.

## Progress

- [x] (2026-07-19 16:46Z) Defined normative scientific conventions in
  `docs/conventions.md`, including pair-dependent weighted interactions.
- [x] (2026-07-19 16:46Z) Added physics-first repository guidance, TDD policy,
  architecture notes, ADR machinery, one architect agent, and three qslib
  skills.
- [x] (2026-07-19 16:46Z) Inventoried the existing ncli physics, exact, TDVP,
  symmetry, observable, and test modules and the standalone Rust SSE crate.
- [x] (2026-07-19 16:46Z) Wrote this initial self-contained qslib 1.0 ExecPlan.
- [x] (2026-07-19 17:00Z) Added complete Apache-2.0 licensing, Cargo metadata,
  security and contribution policies, changelog, Rust 1.85 toolchain policy,
  and a code-of-conduct decision without exposing an inferred personal email.
- [x] (2026-07-19 17:00Z) Drafted and architect-reviewed the seven required
  technical ADRs plus explicit repository-ownership and registry-name ADRs.
  Dependency-sensitive records were accepted after their exact graphs passed
  MSRV, stable, license, source, and advisory validation.
- [x] (2026-07-19 17:00Z) Verified through crates.io that `qslib` is occupied by
  an unrelated active package and that the two initial alternative facade names
  returned no search matches. Publication remains disabled.
- [x] (2026-07-19 17:08Z) Exhausted safe offline Milestone 0 work, reverified
  the ownership state, and paused at the explicit owner gate after three
  consecutive goal turns without the required decisions. Resume from the
  current gate-status list under `External authority and owner gates`.
- [x] (2026-07-19 17:13Z) With explicit owner authorization, initialized
  `qslib/` as a dedicated Git repository on branch `main` and copied the
  architect agent plus all three qslib skills into the new project root. No
  commit, remote, push, or parent-file removal was performed.
- [x] (2026-07-19 17:15Z) Reverified the dedicated repository and paused after
  three consecutive resumed goal turns at the reduced owner gate: remote URL,
  registry names, dependency downloads, conduct reporting, and commits remain
  unauthorized or undecided.
- [x] (2026-07-19 17:16Z) Owner approved remote
  `https://github.com/lere01/qslib.git`, `qslib-quantum` distributions, Rust
  import `qslib`, Python import `qslib_quantum`, repository-scoped dependency
  downloads, closed external contributions, and local commits. Configured the
  remote without network contact; push and pull remain prohibited.
- [x] (2026-07-19 17:34Z) Installed Rust 1.85.0 and validated the selected
  linear-algebra, RNG, checksum, serialization, Parquet/NPY, and Python FFI
  dependency graph on Rust 1.85 and current stable. Cargo-deny license, source,
  and advisory checks pass with one documented unmaintained `paste 1.0.15`
  exception for which no safe transitive upgrade exists.
- [x] (2026-07-19 18:03Z) Began Milestone 1 with the approved nine-package
  workspace and empty default facade feature set. Added eight independent,
  checksummed neutral fixtures, captured the intended compile failure for the
  missing harness API, then implemented the strict parser and made five fixture
  tests plus three workspace-boundary tests green.
- [x] (2026-07-19 17:43Z) Completed Milestone 0: established the dedicated
  repository, ownership boundary, accepted architecture decisions, verified
  the exact staged tree from an isolated reconstruction, and prepared the
  authorized local milestone commit without contacting the remote.
- [x] (2026-07-19 18:09Z) Completed Milestone 1: created the nine-package Rust
  workspace, additive facade features, pinned cross-platform CI definitions,
  strict checksummed neutral fixtures, and an independent conformance and
  dependency-boundary harness. Rust 1.85 and stable checks, docs, feature
  combinations, cargo-deny, and all eight harness tests pass locally.
- [x] (2026-07-19 18:20Z) Completed Milestone 2: implemented checked site
  identifiers and counts, distinct physical and simulation axes, explicit
  binary states and views, multiword little-endian packed states, width-
  labelled byte serialization, reference scalar aliases and finite checks,
  full-basis iteration, and fixed-Hamming-weight sectors. Seven focused core
  tests plus the complete workspace quality suite pass on Rust 1.85 and stable.
- [x] (2026-07-19 18:52Z) Completed Milestone 3: implemented checked
  rectangular and triangular geometry, boundaries, minimum-image shells,
  simple and periodic-image bonds, custom coordinates, explicit x-major
  conversion, weighted interactions, dense/sparse coupling validation, named
  terms, and ADR-0003 disorder provenance. Architect closure review approved;
  Rust 1.85 and stable workspace quality checks pass.
- Observation: the exact staged M3 tree at commit `d0d68f8` was reconstructed
  under `/private/tmp/qslib-m3-candidate.*` and its complete workspace test
  suite passed on Rust 1.85 without relying on the working tree.
- [x] (2026-07-19) Completed Milestone 4: implemented operators, Hamiltonians,
  and canonical model constructors. Architect closure review approved; the
  complete workspace quality suite passes.
- [x] (2026-07-19) Completed Milestone 5: implemented symmetry groups, actions,
  orbit and sector projection utilities, gauges, and resolved-model symmetry
  validation. Architect closure review approved; the complete workspace
  quality suite passes.
- [x] (2026-07-19) Completed Milestone 6: implemented exact full and fixed-sector
  bases, dense and direct CSR matrices, Hermitian spectra and residuals,
  restarted sparse ground states, stable thermal sums, projectors, and signed
  real- and imaginary-time evolution. Architect closure review approved; the
  complete workspace quality suite and CLI ground-state gate pass.
- [x] (2026-07-19) Completed Milestone 7: implemented observables, online statistics,
  autocorrelation diagnostics, and disorder aggregation.
- [x] (2026-07-19) Completed Milestone 8: implemented weighted caller-supplied
  variational local-energy and TDVP kernels, including exact mode signs, dense and
  streamed QGT operations, fixed/GCV/spectral regularization, clipping, residual
  diagnostics, typed solver provenance, and BLAKE3 parameter-layout fingerprints.
  Architect closure review approved; the complete workspace quality suite passes.
- [x] (2026-07-19) Completed Milestone 9: implemented transactional Euler and
  Heun evolution, adaptive Euclidean/QGT acceptance, deterministic retry-safe
  stage seeds, accepted-only observations, overflow-safe commits, and complete
  checkpoint metadata/configuration validation. Architect closure review
  approved; the complete workspace quality suite passes.
- [x] (2026-07-19) Completed Milestone 10: migrated and reconciled the supported
  SSE algorithms on canonical qslib types. Architect review gates cover
  versioned chain seeds, explicit legacy adapters, independent-chain thermal
  confidence, worker-count determinism, cutoff safety, and parity fixtures.
- [x] (2026-07-19 23:51Z) Completed Milestone 11: implemented strict versioned
  configuration, resolved-model reconstruction, manifests, summaries,
  checkpoint envelopes with typed evolution/RNG state and named NPY arrays,
  immutable Parquet trajectory parts with recoverable completion markers,
  legacy migration errors, and an independent Python artifact verifier. Rust
  1.85/stable workspace tests, Clippy, rustdoc, cargo-deny, and focused IO,
  checkpoint-resume, and Python-reader gates pass.
- [ ] Local qslib portion of M12 complete (2026-07-20 00:29Z): the ABI3
  PyO3/NumPy binding has structured exceptions, Maturin metadata, owned-array
  contracts, row-major geometry and coupling resolution, exact basis and
  TFIM ground-state/matrix bindings, Heisenberg/Rydberg matrix bindings,
  convention-labelled observables, and TDVP estimate/solve bindings. Ten
  Python contract tests cover independent four-site and model fixtures,
  real/imaginary-time signs, negative strides, dtype/rank/read-only inputs,
  garbage collection, and concurrent entry. Local Python 3.14 wheel tests,
  Rust 1.85/stable tests, Clippy, rustdoc, formatting, and cargo-deny pass.
  M12 remains open only for the cross-platform CI evidence and the separately
  owned ncli backend-selection/parity adapter.
- [ ] Complete Milestone 12: implement Python bindings and ncli parity adapters.
- [x] (2026-07-20 01:25Z) Completed qslib-local Milestone 13: the CLI now
  uses qslib-io's strict Parquet dataset loader for artifact inspection, rejects
  unknown or unversioned configuration fields, validates model-specific fields
  and command options, reports resolved interaction identities and coefficients,
  emits provenance, and preflights dense exact-memory budgets. Architect closure
  review approved local closure after the strict artifact, schema, option,
  resource, and resolved-model audit.
- [x] (2026-07-20 01:25Z) Completed qslib-local Milestone 14: added an
  installation chapter, combined-site API links, Python example execution in
  CI, a doctested CLI Rust example, corrected CLI mathematics, and a guarded
  documentation-site output path. Architect closure review approved local
  closure after the relative API-link and generated-site checks passed.
- [ ] In progress (2026-07-20 02:05Z): M15 hardening now benchmarks geometry,
  interaction resolution, matrix construction/action, exact diagonalization,
  expectation, TDVP statistics/solve, and short SSE sweeps, with a recorded
  same-host baseline. Structured parser/resolution fuzz smoke coverage and a
  bounded CI invocation are authored. A local previous-commit semver check now
  passes all 165 applicable checks, and workspace llvm-cov reports 78.28% line
  coverage (71.63% regions; branch data is unavailable from this tool output).
  The remote CI matrix and nightly Miri execution remain external gates.
- [ ] Local M16 release-candidate preparation (2026-07-20 01:25Z): optimized
  Rust binaries excluding the Python cdylib, an ABI3 Python wheel, a combined
  Markdown/Rust API site, a portable workspace source archive plus the core
  Cargo package, license/readme/changelog/release-evidence files, and relative
  checksums are available under
  `/private/tmp/qslib-release-candidate-20260720c`. Clean Python-wheel and CLI
  smoke runs pass. This is a 0.1.0 candidate; version 1.0.0 remains prohibited
  until every acceptance gate passes.
- [ ] Complete Milestone 15: complete performance, fuzzing, portability,
  dependency, and API-stability hardening.
- [ ] Complete Milestone 16: build and validate the qslib 1.0 release candidate.

## Surprises and discoveries

- Observation: qslib began as a binary `Hello, world!` template and
  documentation, not a library implementation.
  Evidence: the initial audit found only `src/main.rs`; Milestone 0 has since
  added an intentionally empty documented library target without scientific
  behavior.
- Observation: the modern ncli physics and dynamics paths already have a broad
  independent test corpus that can seed neutral parity fixtures.
  Evidence: the repository currently contains more than two hundred Python test
  functions, including exact Hamiltonian, row-major geometry, observable, and
  TDVP parity tests.
- Observation: the standalone Rust SSE crate already contains validated
  geometry, Hamiltonian terms, TFIM and Rydberg decompositions, samplers,
  artifacts, a CLI, documentation, and parallel-chain support.
  Evidence: `sse/src/lib.rs` exports these public surfaces and `sse/tests/cli.rs`
  exercises the command line.
- Observation: legacy SSE `Spin::Up` cannot be converted universally to a qslib
  bit because TFIM and Rydberg attach different physical meanings to it.
  Evidence: `docs/conventions.md` requires a model-aware compatibility adapter.
- Observation: native Rust neural-network training would duplicate PyTorch GPU
  functionality without addressing the primary computational cost.
  Evidence: `crates/README.md` records that the dominant ncli work is dense GPU
  matrix multiplication and requires profiling before native acceleration.
- Observation: before repository separation, the parent ncli worktree was
  extensively dirty and qslib resolved to the ncli root.
  Evidence: the 2026-07-19 Milestone 0 audit showed `qslib/` and all qslib agent
  assets as untracked parent paths alongside unrelated work. After explicit
  owner approval, `git rev-parse --show-toplevel` from qslib resolves to the
  dedicated qslib root and the required agent assets are self-contained there.
- Observation: `qslib` is not an available Rust or Python distribution name.
  Evidence: `cargo search qslib --limit 10` returned the unrelated crates.io
  package `qslib = "0.15.2"`, and PyPI identifies the same brand as a
  QuantStudio qPCR library. Rust source can retain library target `qslib`, but
  public distributions and the Python import need collision-safe names.
- Observation: the host has current stable Rust 1.96.0 but not the declared
  Rust 1.85.0 MSRV toolchain.
  Evidence: `rustup run 1.85.0 rustc --version` reported that toolchain 1.85.0
  is not installed. Installing it is part of the unresolved dependency-download
  authority gate.
- Observation: the likely remote `https://github.com/lere01/qslib.git` is not
  publicly accessible from the current unauthenticated environment.
  Evidence: a read-only `git ls-remote` returned `Repository not found`; this
  does not distinguish an absent repository from a private one and therefore
  cannot establish remote metadata.
- Observation: only a subset of proposed dependencies is already present in the
  local Cargo source cache.
  Evidence: cached `rand_chacha 0.9.0`, `serde 1.0.229`, `serde_json 1.0.150`,
  and `serde_yaml_ng 0.10.0` declare MSRVs at or below Rust 1.85 and compatible
  licenses. `faer`, BLAKE3, Parquet, NPY, PyO3, and NumPy are not cached, so the
  corresponding ADRs cannot be validated offline.
- Observation: the original Milestone 1 gate required intentionally failing
  tests while `PLANS.md` requires every completed milestone to remain green.
  Evidence: the architect review identified the contradiction. Milestone 1 now
  validates neutral fixtures and harnesses while later milestones record their
  red step immediately before implementation and close green.
- Observation: both selected `faer 0.24.4` and Parquet `59.1.0` pull in
  unmaintained `paste 1.0.15` through their current transitive graphs.
  Evidence: cargo-deny reported RUSTSEC-2024-0436 and showed paths through
  `gemm` and Parquet. The advisory describes maintenance status rather than an
  exploitable vulnerability and states that no safe upgrade is available.
- Observation: Milestone 1 initially could not download Serde and BLAKE3 in the
  restricted sandbox even though the owner had authorized repository-scoped
  dependencies.
  Evidence: Cargo reported DNS resolution failure; the approved elevated rerun
  downloaded the locked Rust 1.85-compatible graph and then exposed the
  intended missing-harness compile error.
- Observation: using the same target name `qslib` for the public library and
  the installed CLI command produces a Cargo rustdoc output collision.
  Evidence: workspace rustdoc warned that the library and binary both wanted
  `target/doc/qslib/index.html`; marking the CLI binary `doc = false` preserves
  the command name and removes the collision while keeping library docs.
- Observation: the committed M1 tree is reproducible without relying on the
  source worktree or untracked generated output.
  Evidence: the staged index was reconstructed under
  `/private/tmp/qslib-m1-candidate.gEDI5F/qslib`; Rust 1.85 and stable checks,
  tests, rustdoc, metadata, license/source policy, links, agent metadata, and
  text-policy checks passed. A fresh local clone at commit `9c7cda5` was clean
  and passed both workspace test suites.
- Observation: packed-state serialization can remain independent of machine
  storage width when the serialized word width is explicit.
  Evidence: M2 round-trips the same four-site state through U8, U16, U32, and
  U64 words and round-trips a 65-site state through two U64 words, rejecting
  nonzero padding and exact-length mismatches.
- Observation: a generic fixed-weight iterator preserves increasing packed
  integer order without eagerly materializing a sector.
  Evidence: the iterator advances a checked ascending position list and passes
  the independent N=4,K=2 vector `[3, 5, 6, 9, 10, 12]`, including K=0 and K=N.
- Observation: the finalized M2 commit is reproducible from a fresh local
  checkout after adding the scalar dependency and direct facade exports.
  Evidence: clone at `931fc67` was clean and passed complete Rust 1.85 and
  stable workspace test suites, including seven core basis tests, the facade
  smoke test, five fixture tests, and three workspace-boundary tests.
- Observation: the initial M2 implementation review exposed portability and
  panic risks that focused tests alone did not reveal.
  Evidence: the architect review required checked full-basis dimensions,
  overflow-safe serialized offsets, a non-lossy large-index error, a dense
  borrowed view, and removal of an internal iterator `expect`; all are now
  covered by the implementation and quality checks.
- Observation: a tree reconstructed solely from the Git index is sufficient to
  expose missing untracked governance or agent files before the first commit.
  Evidence: `git checkout-index --all` produced an isolated candidate under
  `/private/tmp`; both Rust toolchains and all integrity checks passed there.
- Observation: a checkpoint envelope that merely labels opaque bytes is not a
  restart contract.
  Evidence: architect review required typed accepted-state and RNG metadata,
  complete evolution controls, exact array sets, bounded next steps, and a
  real `EvolutionDriver::from_parts` continuation test before M11 closure.
- Observation: Parquet completion is a two-file transaction boundary.
  Evidence: the completed manifest can be durable before its marker; `load`
  now revalidates immutable parts and repairs a missing or incorrect marker,
  and the recovery test exercises that interrupted state.

- Observation: the Python cdylib must be built through Maturin rather than a
  workspace-wide optimized Cargo link.
  Evidence: `cargo +1.85.0 build --release --workspace --all-features` tried
  to link `qslib-quantum-python` as a regular macOS dynamic library and failed
  on unresolved Python symbols; excluding that binding for the Rust binary
  build and running `maturin build` produced the validated ABI3 wheel.
- Observation: the first registry-based `cargo-semver-checks` 0.42.0 attempt
  on macOS exposed a host-specific `system-configuration` null-object panic,
  and a subsequent baseline attempt initially needed authorized registry
  access. A local previous-commit comparison now succeeds without that panic:
  `cargo semver-checks check-release --package qslib-quantum
  --baseline-rev 2584261 --release-type patch` reports 165 passes and 12
  skips. Linux/release-baseline evidence remains useful because this local
  comparison is not a registry publication check.
- Observation: the local Rust 1.85 toolchain does not ship the Miri component.
  Evidence: `cargo miri --version` reports that `cargo-miri` is unavailable;
  the repository now authors a nightly Ubuntu Miri job for the core tests so
  the required safety gate can run on a supported host.
- Observation: compact CLI model input cannot reuse qslib-io's
  `qslib-config-v1` identifier because that schema already names a complete
  `ScientificConfig` document. The two documents have different fields and
  reconstruction guarantees.
  Evidence: the architect audit caught the collision; the CLI now uses
  `qslib-model-input-v1` while retaining shared convention metadata.
- Observation: a documentation generator's destination marker must survive
  mdBook's output replacement. The ownership marker is therefore written
  after mdBook finishes and is required before a later recursive replacement.
  Evidence: a same-destination rebuild now succeeds while an unowned directory
  is refused.
- Observation: Cargo cannot package the entire unpublished workspace as
  registry-ready `.crate` files because path dependencies such as
  `qslib-quantum-core` do not exist on crates.io. A portable source archive and
  the independently packageable core crate are included instead; publication
  remains disabled.
  Evidence: `cargo package --workspace --no-verify` packages core then fails
  resolving the unpublished core dependency for the exact crate.
- Observation: CLI artifact inspection must use the same transactional loader
  as scientific consumers. A filename named `COMPLETE` is not evidence of a
  complete dataset; marker bytes, manifest schema, configuration checksum, and
  every immutable part must be validated together.
  Evidence: the architect's M13 audit found that the previous CLI accepted
  `complete\n`; the red test now rejects it and the green path calls
  `ParquetDatasetManifest::load`.
- Observation: a physics-facing model validation result is incomplete if it
  reports only counts. Pair identities, signed coefficients, site parameters,
  geometry, boundaries, convention schema, and resolved provenance are needed
  to audit a disordered or frustrated Hamiltonian.
  Evidence: the CLI now emits `resolved_interactions`,
  `resolved_specification`, and a provenance object in JSON output.
- Observation: exact CLI commands need a resource contract even when the Rust
  kernels themselves are checked. Dense Hilbert-space matrices can otherwise
  request infeasible allocations before an error is returned.
  Evidence: the CLI now computes checked dimensions and enforces a 256 MiB
  dense-matrix budget before construction, with focused command tests.

## Decision log

- Decision: qslib 1.0 is a finite release contract, not every conceivable
  quantum simulation method.
  Rationale: a bounded definition of done is required for reliable autonomous
  execution and semantic versioning.
  Date/author: 2026-07-19, initial plan.
- Decision: use a layered Cargo workspace with a small facade crate and
  capability crates separated by dependency weight.
  Rationale: core geometry and physics types should remain reusable without
  pulling in eigensolvers, Python, CLI, parallel SSE, or artifact dependencies.
  Date/author: 2026-07-19, initial plan.
- Decision: keep neural-network architectures, autograd, and GPU training in
  ncli for 1.0. qslib variational APIs consume caller-supplied amplitudes,
  ratios, samples, local energies, or log derivatives.
  Rationale: this creates reusable numerical machinery without recreating a
  tensor framework in Rust.
  Date/author: 2026-07-19, initial plan.
- Decision: port SSE behavior into a qslib workspace crate and retain an
  explicit compatibility path for the standalone `quantum-sse` crate during
  migration.
  Rationale: qslib must own canonical basis and interaction types while current
  SSE users need a safe transition.
  Date/author: 2026-07-19, initial plan.
- Decision: external publication is not part of autonomous completion.
  Rationale: pushing, tagging, crates.io publication, package signing, and
  release creation require explicit owner authority. Local verified artifacts
  are sufficient for the release-candidate gate.
  Date/author: 2026-07-19, initial plan.
- Decision: pin ordinary qslib development to Rust 1.85.0 and test current
  stable separately.
  Rationale: edition 2024 begins at Rust 1.85, and continuously exercising the
  MSRV prevents accidental reliance on newer language or standard-library APIs.
  Date/author: 2026-07-19, primary agent after architect review.
- Decision: keep dependency-sensitive ADRs proposed until actual locked
  dependencies pass license, MSRV, and platform checks.
  Rationale: documentation research supports the selected architecture, but an
  accepted dependency decision needs build evidence rather than an assumption.
  Date/author: 2026-07-19, primary agent and qslib-architect.
- Decision: accept ADR-0002 through ADR-0005 after their exact candidate graphs
  compiled on Rust 1.85 and stable and passed license and source policy checks.
  Linux, Windows, and wheel validation remain Milestone 1 and release gates,
  not prerequisites for accepting these architecture choices.
  Rationale: the locked dependency evidence resolves the architecture risk;
  cross-platform execution validates the implementation and packaging later.
  Date/author: 2026-07-19, primary agent after qslib-architect final review.
- Decision: use immutable Parquet trajectory parts with an atomic manifest and
  JSON plus named NPY checkpoint payloads.
  Rationale: append-only parts provide transactional recovery, while NPY keeps
  checkpoint arrays simple, typed, portable, and independently readable from
  Python.
  Date/author: 2026-07-19, primary agent after architect review.
- Decision: use ChaCha20 for qslib 1.x streams with BLAKE3 domain-separated,
  versioned seed derivation.
  Rationale: the algorithm and exact derivation are durable metadata, and
  logical stream identities remain independent of scheduling.
  Date/author: 2026-07-19, primary agent after architect review.
- Decision: do not silently claim the `qslib` Python import while an unrelated
  active distribution owns that namespace.
  Rationale: a distinct distribution name alone does not prevent both packages
  from installing files into the same interpreter namespace.
  Date/author: 2026-07-19, primary agent after registry verification.
- Decision: qslib is a dedicated Git repository rather than an untracked ncli
  subtree.
  Rationale: qslib has independent versioning, CI, consumers, agent guidance,
  and release artifacts, while the parent ncli worktree contains unrelated
  active changes. A separate root gives every qslib operation an unambiguous
  ownership boundary.
  Date/author: 2026-07-19, owner approval implemented by the primary agent.
- Decision: use `qslib-quantum` for Rust and Python distributions, `qslib` for
  the Rust library target, and `qslib_quantum` for the Python import package.
  Rationale: the selected distribution is collision-safe while preserving the
  qslib brand. The Python import must also avoid the unrelated active `qslib`
  namespace.
  Date/author: 2026-07-19, owner approval implemented by the primary agent.
- Decision: keep unsolicited external contributions closed and do not adopt a
  public code of conduct at this project stage.
  Rationale: opening contributions requires a private reporting contact
  approved for publication and a later explicit governance decision.
  Date/author: 2026-07-19, owner.
- Decision: local initial and milestone commits are authorized, while push and
  pull are prohibited.
  Rationale: local commits provide recoverable milestone boundaries without
  causing external state changes.
  Date/author: 2026-07-19, owner.
- Decision: accept a monitored cargo-deny exception for RUSTSEC-2024-0436 until
  `faer` and Arrow/Parquet remove `paste 1.0.15`.
  Rationale: `paste` is used as a compile-time macro, the advisory reports
  unmaintained status rather than a vulnerability, both selected foundational
  graphs require it, and no safe upgrade exists. The exception is explicit and
  must be removed when upstream migrates.
  Date/author: 2026-07-19, primary agent after architect confirmation.
- Decision: use prefixed Cargo package identifiers with ergonomic Rust targets,
  retain the unpublished name `qslib-test-support`, and keep the facade default
  feature set empty.
  Rationale: this reconciles ADR-0001's conceptual crate names with ADR-0009,
  preserves concise imports, and makes core-only dependency cost explicit.
  Date/author: 2026-07-19, primary agent after qslib-architect review.
- Decision: store neutral conformance evidence in
  `fixtures/conformance/v1/` with strict case-specific payload validation and a
  BLAKE3 manifest.
  Rationale: a language-neutral, checksummed oracle can be reused by Rust,
  Python, and future backends without depending on production implementation.
  Date/author: 2026-07-19, primary agent after qslib-architect and conformance
  review.
- Decision: keep core basis serialization explicit and width-labelled without
  deriving durable Serde layouts, keep model-specific occupation aliases out
  of generic `BasisBit`, and expose concrete `Real`, `Complex64`, and finite
  scalar validation helpers.
  Rationale: ADR-0004 assigns durable schemas to `qslib-io`; explicit bytes
  preserve little-endian semantics without accidental schema promises. Generic
  zero/one values remain distinct from Rydberg occupation and legacy spin
  labels, while numerical code has one reference scalar policy.
  Date/author: 2026-07-19, primary agent after qslib-architect review.
- Decision: make resolved IO reconstruction lossless for serialized physical
  semantics by applying interaction multiplicity, restoring physical constants,
  and rejecting duplicate onsite roles at one site.
  Rationale: a seed or positional coefficient list cannot reconstruct a unique
  Hamiltonian; the durable term table must determine the same operator.
  Date/author: 2026-07-19, primary agent after qslib-architect review.
- Decision: checkpoint accepted-state metadata stores the complete trajectory
  controls and treats RNG position as complete 64-byte ChaCha20 blocks, with
  explicit conversion to the 16-word stream position.
  Rationale: `EvolutionDriver::from_parts` compares every trajectory-changing
  control, and durable restart must reject unsupported seed derivation versions
  rather than silently reinterpret them.
  Date/author: 2026-07-19, primary agent after qslib-architect review.
- Decision: the first Python binding surface is coarse grained and returns
  owned NumPy arrays rather than mutable Rust handles or per-sample callbacks.
  Rationale: this keeps scientific semantics in Rust, makes buffer lifetime
  explicit, and gives ncli a stable backend boundary before TDVP and broader
  parity adapters are added.
  Date/author: 2026-07-20, primary agent after ADR-0005 and M12 inventory.
- Decision: CLI configuration uses YAML or JSON with one explicit dense
  coupling matrix and model-specific site arrays; every upper-triangle entry
  remains a resolved pair coefficient, including zero values. This keeps
  pair-dependent disorder and normalization visible at the command boundary.
  Date/author: 2026-07-20, primary agent after M13 contract tests.
- Decision: release-candidate Rust builds exclude the Python cdylib target and
  invoke Maturin separately. This is the supported PyO3 packaging boundary and
  avoids pretending that an extension module is an ordinary standalone Rust
  dynamic library on macOS.
  Date/author: 2026-07-20, primary agent after the optimized-link audit.
- Decision: CLI configurations use the distinct physicist-facing
  `qslib-model-input-v1` envelope plus the shared `qslib-conventions-v1`
  metadata, require canonical `row_major` and
  `little_endian` metadata, reject unknown and model-inapplicable fields, and
  include a resolved provenance object in stable JSON results.
  Rationale: a typo or omitted basis convention must never silently become a
  different scientific run, and downstream artifact consumers need to audit
  the resolved model without reconstructing it from counts.
  Date/author: 2026-07-20, primary agent after architect M13 audit.
- Decision: `qslib artifacts inspect` is a strict Parquet dataset inspection
  command, not generic filesystem metadata. It delegates marker, schema,
  checksum, and part validation to `qslib-io` and returns manifest provenance.
  Rationale: duplicate artifact semantics in a CLI would drift from the
  durable IO contract and make incomplete results look publishable.
  Date/author: 2026-07-20, primary agent after architect M13 audit.

## Outcomes and retrospective

Milestone 0 is complete. The project has a
dedicated repository root, self-contained agent assets, complete licensing and
baseline governance, an explicit library target, a pinned MSRV policy, and
architect-reviewed architecture records. No scientific production behavior has
been implemented. All owner choices are resolved, Rust 1.85 and stable checks
pass, all nine initial ADRs are accepted, and the dependency policy is
machine-readable. The architect's final review found no scientific or owner
gate. The exact staged tree was reconstructed outside the working tree and
passed the MSRV, stable, metadata, documentation-link, agent-metadata, license,
source, and text-policy checks.

Milestone 1 is now complete. A physicist-facing contributor can inspect the
approved package map and feature costs, run the same local checks used by CI,
and examine eight independently derived small-system evidence files before
any model implementation exists. The harness rejects unknown schemas,
unsupported conventions, missing provenance, malformed matrices, non-finite
values, duplicate identities, wrong checksums, and incomplete fixture sets.
Linux, macOS, and Windows workflow definitions are authored and action-pinned;
their remote execution remains a later owner-authorized CI observation because
push and pull remain prohibited. Milestone 2 then started with failing public
tests for checked identifiers, basis states, and deterministic packing.

Milestone 2 is now complete. Core users can construct and validate non-empty
systems, distinguish physical axes from the simulation basis, inspect dense
binary states, pack arbitrary supported site counts little-endian by site,
serialize with an explicit word width, enumerate the full basis, and enumerate
fixed-occupation sectors deterministically. Invalid bits, empty systems,
out-of-range sites, non-canonical high bits, wrong byte lengths, overflowed
dimensions, non-finite scalars, and invalid sector weights return structured
errors. Milestone 3 is complete: its red-to-green tests cover geometry,
boundaries, canonical bonds, pair-dependent interactions, disorder provenance,
parity adapters, and anisotropic closest-vector cases. Milestone 4 now starts
with operator channels, Hamiltonian term assembly, and model constructors.

Milestones 4 through 10 are also complete, covering canonical Hamiltonians,
symmetry sectors, exact spectra and dynamics, observables and statistics,
variational/TDVP kernels, transactional integration, and the migrated SSE
backend. Milestone 11 is complete: physicists can now persist strict resolved
configurations and summaries, bind artifacts to conventions and checksums,
write portable checkpoint arrays and typed restart metadata, append immutable
Parquet trajectories, recover interrupted dataset publication, and verify the
result independently from Python. Milestone 12 is the next active boundary;
publication and remote repository changes remain intentionally outside this
autonomous run.

M12's qslib-local implementation is complete. The ABI3 wheel imports as
`qslib_quantum`, validates strided and Fortran-order inputs, resolves
pair-dependent couplings, reports convention-labelled observables, and
reproduces exact TFIM, Heisenberg, and Rydberg matrix behavior through owned
NumPy outputs. Cross-platform wheel jobs are authored but cannot run until CI
executes, and ncli backend adoption remains a separate ownership boundary
that must not modify the parent repository without explicit authority. The
qslib-local CLI milestone is now complete: documented four-site and tiny SSE
commands execute through public kernels, and JSON output is tested as a stable
machine-facing surface. M15 hardening and M16 local artifact preparation are
in progress. The local release candidate has a reproducible checksum manifest
and clean-environment smoke evidence; local semver and coverage evidence are
now recorded, while cross-platform CI, Miri execution, and ncli backend parity
remain open gates.

## Context and orientation

The qslib project root is its dedicated Git repository at `qslib/`. Its
authoritative files are `AGENTS.md`, `PLANS.md`, `docs/conventions.md`,
`docs/architecture/README.md`, and this ExecPlan. Project-scoped agent assets
live at `.codex/agents/qslib-architect.toml` and `.agents/skills/qslib-*` inside
that root. The surrounding ncli repository remains a behavioral reference and
contains unrelated active work that qslib operations must not modify.

The existing Python implementation is a behavioral reference, not the qslib
architecture. Important source areas are:

- `src/ncli/physics/hamiltonians/` for TFIM, Heisenberg J1-J2, Rydberg,
  lattice-pair construction, and local-energy kernels;
- `src/ncli/physics/observables/` for energy moments, magnetization,
  correlations, structure factors, Fisher information, Shannon entropy, and
  entanglement estimators;
- `src/ncli/dynamics/exact.py` and `src/ncli/dynamics/exact_tfim.py` for exact
  bases, sparse matrices, ground-state references, and exact time evolution;
- `src/ncli/variational/tdvp/` for parameter layouts, QGT statistics, dense and
  matrix-free solvers, regularization, replay, and diagnostics;
- `src/ncli/dynamics/integrators.py`, `state.py`, and `checkpoint.py` for
  accepted-state integration and recovery;
- `src/ncli/models/common/symmetry/` and
  `src/ncli/wavefunction/symmetry.py` for existing group actions and gauges;
- `tests/` for independent behavioral anchors and regression cases.

The existing Rust reference is `sse/`. Its core source lies in
`sse/src/geometry.rs`, `lattice.rs`, `hamiltonian.rs`, and `sse/`. It currently
supports TFIM and Rydberg SSE, parallel deterministic chains, configuration,
artifacts, a CLI, and documentation. Its physical conventions must be converted
explicitly rather than copied blindly.

The canonical qslib site order is row-major, `site = x + Lx*y`. Packed states
are little-endian by site. In a selected Pauli simulation basis, bit zero has
eigenvalue plus one and bit one has eigenvalue minus one. Rydberg occupation is
the bit value, `n = (1-Z)/2`. Couplings belong to resolved interaction terms and
may vary by pair. These meanings are release-critical.

## Scope and non-goals

qslib 1.0 includes the following supported capabilities:

- validated chain, rectangular, square, triangular, and custom geometries with
  open and periodic boundaries where the geometry defines them;
- simple and explicit periodic-image bond multiplicity;
- dense and packed two-level basis states, fixed-Hamming-weight sectors, and
  explicit simulation bases;
- weighted per-site and per-pair operator terms, constants, deterministic term
  combination, and model provenance;
- inhomogeneous TFIM, isotropic Heisenberg exchange, homogeneous and disordered
  J1-J2 geometry, and driven Rydberg models;
- pair-supplied, matrix-supplied, shell-generated, and seeded disorder coupling
  inputs, with realized values persisted;
- lattice transformations, finite symmetry groups, spin inversion, characters,
  orbit representatives, invariant sectors, and projection utilities;
- full and conserved-sector exact bases, dense and sparse Hermitian matrices,
  ground and low-lying eigensystems, exact thermal sums for small systems, and
  real- and imaginary-time exact evolution;
- energy, variance, Pauli and spin magnetization, raw and connected
  correlations, structure factors, sublattice order, total spin, Fisher
  information, Shannon entropy, and exact pure-state bipartite entropy;
- online weighted moments, multiple-chain aggregation, autocorrelation and
  effective-sample diagnostics, and separate disorder-ensemble aggregation;
- local-energy evaluation from caller-provided amplitudes or ratios, QGT and
  force statistics from caller-provided log derivatives, dense and
  matrix-vector TDVP solves, documented regularizers, and transactional
  integrators;
- sign-safe SSE for the models and parameter regimes validated by the migrated
  implementation, initially TFIM and Rydberg, with explicit rejection of
  unsupported sign structures;
- versioned YAML and JSON configuration, JSON summaries, a columnar trajectory
  format, checksummed artifacts, and restartable accepted-boundary checkpoints;
- Rust APIs, Python bindings for stable scientific kernels, a CLI, examples,
  rustdoc, a physics-first guide, and local release artifacts.

The following are post-1.0 unless an accepted ADR replaces a listed 1.0 item
without expanding the schedule:

- native Rust neural-network architectures, autograd, optimizers, or GPU
  training;
- a general tensor library or direct replacement for PyTorch;
- matrix-product states, PEPS, tensor-network contraction, circuit simulation,
  quantum hardware control, or chemistry-specific integral pipelines;
- a general cure for the SSE sign problem or an assertion that arbitrary
  disordered Heisenberg models are sign safe;
- Ewald summation and arbitrary periodic simulation cells beyond the validated
  1.0 geometry contract;
- distributed multi-GPU TDVP or production-scale exact diagonalization;
- every observable and model present in the literature;
- remote publication or deployment without owner authorization.

## Architecture and interfaces

Convert the qslib package into a Cargo workspace while retaining the root
package as the public `qslib` facade. Conceptual crate labels are followed by
their Cargo package identifiers where those differ:

- `qslib-core` (`qslib-quantum-core`): identifiers, scalar and error policy,
  basis states, geometry,
  weighted interactions, operators, canonical models, symmetry primitives, and
  observable definitions;
- `qslib-exact` (`qslib-quantum-exact`): basis enumeration, matrix construction,
  Hermitian eigensolvers,
  exact thermodynamics, and exact evolution;
- `qslib-variational` (`qslib-quantum-variational`): local-energy aggregation,
  QGT and force statistics,
  TDVP linear solves, diagnostics, and integrators;
- `qslib-sse` (`qslib-quantum-sse`): SSE decomposition, state, updates,
  measurements, deterministic
  parallel chains, and thermal diagnostics;
- `qslib-io` (`qslib-quantum-io`): versioned configuration, manifests,
  artifacts, checksums, and
  checkpoints;
- `qslib-python` (`qslib-quantum-python`): PyO3 and NumPy bindings for
  intentionally stable kernels;
- `qslib-cli` (`qslib-quantum-cli`): physicist-first commands built only on
  public library APIs;
- `qslib-test-support`: unpublished neutral fixtures, independent tiny-matrix
  builders, and legacy parity loaders used only by tests.

The root `qslib` facade reexports `qslib-core` and feature-gates heavier
capabilities. Core-only use must not compile eigensolver, Rayon, Python, CLI, or
columnar-I/O dependencies. Provide a documented `full` feature for users who
want all Rust capabilities. Python and CLI remain separate packages rather than
facade features.

Use `f64` and `Complex64` for reference public numerics. Lower precision may be
accepted only by APIs that record it explicitly. Keep unsafe Rust forbidden in
core scientific crates. Any unavoidable FFI unsafe code stays inside
`qslib-python`, is documented, and has boundary tests.

The initial stable vocabulary should include checked equivalents of `SiteId`,
`Axis`, `BasisBit`, `PackedState`, `SectorBasis`, `BoundaryCondition`,
`Geometry`, `Bond`, `InteractionId`, `WeightedInteraction`, `OperatorChannel`,
`HamiltonianTerm`, `Hamiltonian`, and model-specific validated builders. Exact,
variational, and SSE crates consume these shared meanings rather than defining
parallel spin or bond types.

The implementation must settle the open convention decisions through ADRs
before the affected production code. Required early ADRs cover workspace
boundaries, linear algebra and sparse storage, deterministic RNG and seed
derivation, serialization and columnar artifacts, Python FFI ownership, public
stability and feature flags, and migration of the standalone SSE crate. Each
ADR must record evaluated alternatives and a testable consequence.

If no maintained sparse Hermitian eigensolver satisfies the accepted ADR,
implement a qslib-owned Lanczos reference with deterministic initialization,
full or selective reorthogonalization, residual diagnostics, degeneracy tests,
and a dense fallback. Do not silently call a non-Hermitian solver or discard
complex phases.

## Milestones

### Milestone 0: repository and architecture gate

Establish qslib as a clean ownership and versioning unit before substantial
code is written. The preferred outcome is a dedicated qslib Git repository with
the current `.agents` skills and `.codex` agent copied into its root. If the
owner elects to keep qslib in the ncli repository through 1.0, record that
decision and ensure every command and commit path is scoped to `qslib/` plus the
project agent assets.

Add the Apache-2.0 license, package metadata, security policy, contribution
guide, code-of-conduct decision, changelog, rust-toolchain policy, and initial
ADRs named above. Choose an MSRV no lower than Rust 1.85 because edition 2024 is
in use; raise it only when an accepted dependency decision requires that.
Resolve whether the qslib name is available for intended publication, but do
not publish. The gate passes when a clean checkout has one unambiguous project
root, all guidance is discoverable, licensing is complete, ADRs are accepted,
and no production crate decision remains implicit.

Current evidence at 2026-07-19 17:34Z:

- Rust 1.85 and current stable formatting, locked checks, Clippy with warnings
  denied, all-target tests, and rustdoc with warnings denied pass for the
  architecture scaffold. There are no scientific tests or behaviors yet.
- Cargo metadata reports package `qslib-quantum`, library target `qslib`,
  Apache-2.0, edition 2024, Rust 1.85, the approved remote URLs, and publication
  disabled.
- `LICENSE` is byte-identical to the complete Apache-2.0 text already used by
  the standalone SSE project.
- Architect TOML, all three qslib skill front matters, and all three skill
  interface metadata files parse successfully.
  An initial validation command used paths relative to the wrong directory and
  assumed unavailable PyYAML; the corrected check used parent-root paths,
  Python `tomllib`, and Ruby's standard YAML parser.
- The no-en-dash and no-em-dash scan passes for qslib and its project agent
  assets.
- Every relative Markdown link in the qslib documentation resolves to an
  existing local target.
- The initial and final independent `qslib-architect` reviews are complete.
  Their corrections are incorporated in ADRs and this plan. The final review
  found no scientific-convention conflict or additional owner gate.
- The dedicated repository, self-contained agent assets, remote metadata, final
  package identifiers, closed-contribution decision, and commit authority are
  verified.
- A temporary dependency probe pinned every selected dependency family and
  compiled its used surface on Rust 1.85 and stable. `cargo deny check
  advisories licenses sources` passes with the documented RUSTSEC-2024-0436
  maintenance exception. This evidence supports acceptance of ADR-0002 through
  ADR-0005. Linux, Windows, and wheel execution remain later implementation and
  release validation gates.
- Immediately before the final gate review, `cargo search qslib-quantum`
  returned no crates.io match and the PyPI JSON endpoint returned HTTP 404.
  This verifies current availability but does not reserve either name.
- The final self-containment check initially used `path` as a zsh loop variable,
  which temporarily replaced that subprocess's command search path. The
  corrected command used `skill_name`; all six skill files match their source,
  parse successfully, and no repository state was affected by the failed check.
- The exact staged candidate was reconstructed with `git checkout-index --all`
  into `/private/tmp/qslib-m0-candidate.Nr8BU0/qslib`. Rust 1.85 and stable
  formatting, locked check, Clippy, tests, and rustdoc passed there. Cargo-deny
  license and source checks, Cargo metadata assertions, agent TOML/YAML parsing,
  relative Markdown-link validation, and the en-dash/em-dash scan also passed.

### Milestone 1: workspace and conformance harness

Replace the binary template with the workspace and facade structure described
above. Add CI definitions for Linux, macOS, Windows, stable Rust, and MSRV.
Configure formatting, clippy with warnings denied, rustdoc with missing public
documentation denied, dependency license checks, and a test-support crate.

Before implementing scientific types, add neutral JSON fixtures for rectangular
indexing, bit packing, one-bond TFIM and Heisenberg matrices, heterogeneous
Heisenberg couplings, two-site Rydberg energies, observable normalization, and
basis-rotation spectrum parity. Milestone 1 tests validate fixture schemas,
independent-oracle provenance, and harness behavior while remaining green.
Fixtures live under `fixtures/conformance/v1/`, use strict case-specific
payload validation, record complex `f64` matrix shape and row-major layout, and
are bound to a sorted BLAKE3 manifest. The harness checks schema versions,
provenance, comparison policy, dimensions, finite values, Hermiticity where
claimed, completeness, uniqueness, and checksums without calling a production
qslib crate.
Each later implementation milestone adds its convention tests first, records
their intended failing run, and makes them green before that milestone closes.
The gate passes when the workspace scaffolds cleanly, fixture and harness tests
pass, and CI definitions run the same commands locally. No completed milestone
retains a required failing test.

### Milestone 2: basis and foundational types

Implement checked identifiers, axes, basis bits, dense binary state views,
packed states, checked word-width serialization, Hamming weight, full basis
iteration, and fixed-weight sector enumeration. Design APIs so Rydberg
occupation is not an alias for an ambiguous `Spin::Up` value. Add overflow,
invalid-bit, maximum-size, ordering, round-trip, property, and compile-time API
tests first.

The gate passes when every basis conformance vector is green, full and sector
bases enumerate in canonical order, random valid round trips pass property
tests, invalid states return structured errors, and no model-specific meaning
leaks into the generic bit type.

### Milestone 3: geometry, interactions, and disorder

Implement the canonical geometries, coordinate mappings, boundaries,
minimum-image distances, pair shells, bonds, multiplicity, interaction
identities, weighted channels, dense and sparse coupling inputs, deterministic
canonicalization, and disorder realization provenance. Write tests for tiny
periodic degeneracies, nonsquare row-major mappings, triangular coordinates,
custom geometry, duplicate terms, shell overlap, signed and zero couplings,
symmetric-matrix validation, and scheduling-independent generation before code.

The gate passes when both modern ncli row-major fixtures and canonical SSE
geometry cases agree after explicit conversion, while the legacy x-major path
fails without its named adapter and succeeds with it. Persisted disorder tables
must reproduce without access to the generating RNG.

Execution record (2026-07-19): the initial M3 acceptance tests were added before
production implementation. The focused Rust 1.85 run failed at compilation with
unresolved public geometry and interaction types, which is the intended red
state for this test-first transition. The tests independently specify row-major
rectangular and triangular coordinates, mixed boundaries, periodic-image
multiplicity, canonical weighted terms, duplicate rejection, dense/sparse
coupling parity, and scheduling-independent disorder realization.

Implementation record (2026-07-19): M3 production behavior is now present in
`qslib-core`. Geometry retains explicit periodic-image identity while simple
bonds are endpoint-only; extent-one periodic steps are skipped as self-bonds;
triangular minimum-image selection derives a cell-local closest-vector search
from the actual cell dimensions; shells use a typed absolute or relative
tolerance; and `XMajorAdapter` is the only legacy-order conversion. Interaction
tables preserve named terms and zero-coefficient provenance while exposing an
active numerical view. Dense and sparse coupling constructors validate finite
values, shape, symmetry, diagonals, duplicates, and canonical ordering.
Disorder uses keyed BLAKE3 `qslib-seed-v1` framing and `ChaCha20`, with logical
realization index and structured provenance. Focused geometry, interaction, and
neutral modern-row-major/SSE parity tests are green on Rust 1.85 and stable;
workspace Clippy, tests, rustdoc, and cargo-deny license/source/advisory checks
also pass offline or against the approved local advisory cache.

The final M3 corrections are now implemented: triangular minimum-image search
uses a norm-derived complete candidate interval, ADR-0003 seed framing is exact
(`qslib-seed-v1` plus the required zero byte, u32 domain and index framing), and
all fields participating in periodic-image identity feed a canonical identity
fingerprint. Pinned seed and anisotropic closest-vector tests are green, and
the current architect re-audit is the remaining closure check. Milestone 4
will begin with operator channels and Hamiltonian term assembly only after that
gate closes.

### Milestone 4: operators, Hamiltonians, and models

Implement local operator channels and a deterministic Hamiltonian term list with
constants and per-term complex or real coefficients. Implement validated
builders for inhomogeneous TFIM, isotropic Heisenberg exchange, homogeneous and
disordered J1-J2 geometry, and Rydberg systems with per-site drive and detuning
and per-pair interaction.

Write explicit matrices and connected-state tests before each builder. Verify
Hermiticity, signs, double-counting policy, basis-aware action, local-energy
ratios, homogeneous shorthand expansion, heterogeneous coefficients, and
invalid configurations. The gate passes when matrices and local energies match
independent Python fixtures on small systems and every model can be inspected
as a resolved weighted-term table.

Execution record (2026-07-19): operator and model acceptance tests were added
before implementation. The focused Rust 1.85 run failed at compilation with
unresolved `PauliString`, `ModelError`, and TFIM/Heisenberg/Rydberg builder APIs,
which is the intended M4 red state.

The red contract was then expanded before production work: typed Ising,
Heisenberg, and Rydberg channels; explicit Hermitian matrix checks; TFIM x-basis
rotation; signed heterogeneous three-site exchange; J1-J2 shell shorthand;
ordered Pauli-product reduction; non-Hermitian rejection; and deterministic
duplicate-term handling are now required by the tests.

Implementation record (2026-07-19): M4 now supplies canonical Pauli action,
ordered-product phase reduction, constant folding, deterministic duplicate
combination, Hermitian coefficient validation, explicit matrix application,
and local-energy evaluation with the correct Hermitian matrix-element
orientation. Missing connected amplitudes are errors. TFIM, isotropic
Heisenberg, J1-J2, and driven Rydberg builders validate typed interaction
channels and finite inputs, preserve signed pair coefficients, support the
documented simulation bases, and return `ResolvedModel` provenance wrappers.
The wrapper retains model family, basis, weighted interaction identities, and
site parameters while dereferencing to `Hamiltonian` for numerical consumers.
Independent tests cover heterogeneous pair terms, Rydberg occupation energies,
J1/J2 exact shell counts, triangular rejection, Hermiticity, complex Pauli
local energy, identity reduction, and basis-spectrum invariants. Focused and
workspace Rust 1.85 tests, Clippy with warnings denied, and rustdoc with
warnings denied are green; architect closure review approved M4.

The closure corrections are also implemented: local energy now preserves a
general complex coefficient while conjugating only the Hermitian Pauli matrix
element; duplicate floating-point contributions use a canonical coefficient
ordering; `j1j2_disordered` validates exact shell lengths and preserves signed,
zero, and overlapping shell identities; `ResolvedModel` retains a typed
`ModelSpecification`; and qslib-core tests compare complete model matrices and
local-energy rows against the independent neutral fixtures in
`fixtures/conformance/v1/`. The final architect audit approved M4.
The final qslib-architect audit approved M4 after rerunning 55 qslib-core
integration tests, Clippy with warnings denied, and rustdoc with warnings
denied.

### Milestone 5: symmetry and sectors

Implement site permutations, composition, inverse, finite groups, translations,
rectangle and square point groups, spin inversion, characters, orbit
representatives, canonicalization, sublattice gauges, and group projection.
Separate lattice symmetries from Hamiltonian symmetries by validating that a
candidate action preserves the resolved term list.

Write group-law, permutation, projector-idempotence, character, orbit,
Hamiltonian-commutator, and legacy-adapter tests first. The gate passes when
small exact spectra split into symmetry sectors whose combined multiplicities
match the unrestricted spectrum and rejected actions explain which term breaks
the symmetry. The exact-spectrum multiplicity check depends on the M6 exact
matrix/eigensolver layer; M5 must provide the group, orbit, projection, and
sector-index contracts that M6 consumes, while M5 tests cover group-law and
projector invariants independently.

Execution record (2026-07-19): the M5 acceptance tests were added before
production implementation. The focused Rust 1.85 run failed at compilation
with unresolved permutation, translation, finite-group, spin-inversion,
character, and interaction-symmetry APIs, which is the intended red state.

Implementation record (2026-07-19): M5 symmetry primitives are now present in
`qslib-core`. Gather-direction permutations validate bijections, composition,
inverse, state action, Pauli support mapping, and simple-bond mapping. Geometry
builders cover translations, finite translation groups, rectangle and square
point groups, with explicit open-boundary rejection and periodic wrapping.
Spin inversion, deterministic orbit representatives, trivial-character
projection, and resolved weighted-interaction symmetry checks are covered by
independent tests. Characters are bound to the canonical ordered permutation
set used for validation, and projection rejects reuse with another group even
when the group order matches. Architect closure review approved M5 after the
focused symmetry target passed all 11 tests and the workspace quality suite
passed with formatting, warnings-denied Clippy, and warnings-denied rustdoc.

### Milestone 6: exact ground states and dynamics

Implement full and fixed-sector bases, dense and CSR Hamiltonians, matrix-vector
action, dense Hermitian diagonalization, sparse extremal eigensolving, residual
diagnostics, degenerate-subspace comparison, exact thermal sums, matrix
exponentiation or Krylov evolution, and normalized imaginary-time evolution.

Write tests from analytic spectra and independently assembled matrices first.
Port neutral fixtures from `tests/test_exact_ed_vmc.py`,
`tests/test_exact_quench.py`, and `tests/test_dynamics_exact_gate.py`. Verify
dense-sparse parity, full-sector projection parity, residuals, unitary norm,
energy conservation, second-order or declared integrator convergence, and
Hadamard basis parity. The user-visible gate is a Rust example and CLI-internal
test that obtains the correct ground state for a heterogeneous four-site model.

Execution record (2026-07-19): the M6 acceptance tests were added before
production implementation. The focused Rust 1.85 run failed at compilation
with unresolved exact-basis, matrix, eigensolver, thermal, and evolution APIs,
which is the intended red state.

Implementation record (2026-07-19): `qslib-quantum-exact` now provides ordered
full and fixed-weight bases, dense and CSR matrix construction from the core
Hamiltonian action, matrix-vector products, Hermitian validation, a
deterministic complex-Hermitian reference eigensolver, residual-reporting
ground states, exact thermal summaries, and spectral real- and imaginary-time
evolution. Independent tests cover canonical ordering, dense-CSR parity,
analytic one-site spectra, complex Pauli-Y Hermiticity, non-Hermitian rejection,
thermal sums, normalized evolution, and a heterogeneous four-site ground
state. The sparse extremal path uses deterministic fully reorthogonalized
Lanczos with basis-vector restarts, direct CSR assembly, fixed-sector rejection,
degenerate-projector checks, stable thermal log-sums, signed real-time
evolution, Hadamard spectrum parity, and norm/energy invariants are covered by
independent tests. The heterogeneous four-site Rust example is also green.
The CLI dependency graph gate is also green, and the architect closure review
approved M6 after all 15 exact integration tests passed.

### Milestone 7: observables and statistics

Implement observable definitions and exact or sample estimators for the 1.0
scope. Implement stable weighted online moments, within-chain and between-chain
aggregation, autocorrelation estimates, effective sample size, R-hat with a
named algorithm, and disorder-realization aggregation.

Write product-state, singlet, analytically correlated, exact-enumeration, and
synthetic time-series tests first. Test axes and normalization explicitly. The
gate passes when exact observables match direct matrix evaluation, online and
batch accumulation agree, correlated-chain diagnostics detect a constructed
failure, and disorder uncertainty remains separate from sampling uncertainty.

Execution record (2026-07-19): the M7 statistics and exact-observable tests were
added before production implementation. The focused Rust 1.85 run first failed
with unresolved weighted-moment, complex-statistics, autocorrelation, R-hat,
disorder-aggregation, expectation, and variance APIs, which is the intended red
state.

Implementation record (2026-07-19): `qslib-quantum-variational` provides stable
weighted real and complex online moments with unavailable empty-state errors,
transactional overflow handling, merge invariance, named Geyer
initial-positive-sequence autocorrelation with common-`N` normalization and
effective sample size, structured classic R-hat retaining within/between
variances, and disorder-realization records that reject duplicate IDs, retain
fixed-realization weights, and propagate sampling uncertainty with squared
normalized weights. `qslib-quantum-exact` provides normalized pure-state and
direct Pauli expectations, centered-norm Hermitian variance, Shannon and
bipartite entropy, Pauli and `S=sigma/2` magnetization, exact two-site
correlations and `<S_tot^2>`, raw/connected correlations, labelled
position/q-bound structure factors, weighted sublattice estimators, QFI, and
thermal observables. Product, singlet, exact-enumeration, direct-matrix,
analytic, boundary, overflow, and synthetic-chain tests are green. Full Rust
1.85 workspace tests, Clippy with warnings denied, rustdoc with warnings
denied, formatting, and diff checks pass. Architect closure review approved
M7; M8 is now active.

### Milestone 8: variational and TDVP numerical core

Define interfaces that accept samples, normalized or unnormalized weights,
local energies, amplitude ratios, and caller-computed complex log derivatives.
Implement local-energy aggregation, energy and variance, QGT and force
statistics, dense QGT storage, matrix-vector products, conjugate gradient,
fixed Tikhonov, GCV selection, documented spectral or SNR filtering, direction
clipping, residual diagnostics, and deterministic parameter-layout
fingerprints.

The Rust library does not compute neural-network derivatives in 1.0. Python
bindings adapt Torch or NumPy-produced arrays after validation. Port the
independent formulas and parity cases from `tests/test_variational_tdvp.py` as
neutral fixtures. The gate passes when dense, matrix-vector, and streamed
solutions agree on manufactured positive-semidefinite problems and real- and
imaginary-time signs match the convention document.

Implementation record (2026-07-19): weighted and unweighted caller-supplied
statistics, complex force components, row-oriented local-energy ratio
aggregation, streamed QGT products, dense QGT validation and eigenspectra,
relative-tolerance conjugate gradients, fixed Tikhonov/GCV/spectral-cutoff
regularization, update clipping, projected residual and roundoff-floor
diagnostics, typed mode/representation/solver provenance, and BLAKE3-v1
parameter-layout fingerprints are implemented in `qslib-variational`. The
reference tests cover analytic QGT/force values, complex-coefficient local
energy, streamed/dense parity, scale-invariant spectral cutoffs, degenerate
maximum modes, invalid numerical inputs, and solver diagnostics. The workspace
test, Clippy, rustdoc, formatting, and diff gates pass; the qslib architect
approved closure.

### Milestone 9: integration and checkpointable evolution

Implement a generic velocity callback, Euler for reference, Heun, adaptive
accepted-state transactions, Euclidean and QGT metrics, deterministic stage
seeds, rejection without state mutation, trajectory observation boundaries,
and serializable evolution metadata. Keep model parameter storage outside the
integrator through a checked flat-state abstraction.

Write constant and linear ODE convergence tests, rejection rollback tests,
resume equivalence tests, RNG-stream tests, and exact small-quench comparisons
first. The gate passes when an interrupted and resumed evolution produces the
same accepted trajectory and diagnostics as an uninterrupted run under the
declared reproducibility policy.

Implementation record (2026-07-19): `qslib-quantum-variational` now exposes
checked flat-state Euler and Heun drivers with deterministic stage seeds,
adaptive growth and shrinkage, scale-relative QGT metric validation, complete
JSON checkpoint metadata, accepted-boundary observations, and transactional
rollback. Independent tests cover analytic one-step values, second-order
convergence, QGT norms and indefinite spectra, actual rejection seed invariance,
overflow-safe failure behavior, and adaptive interrupted/resumed trajectories.
Architect closure review approved; workspace tests, Clippy, rustdoc, formatting,
and diff checks pass.

### Milestone 10: SSE migration

Port the standalone SSE implementation into `qslib-sse` through canonical core
types. Do not copy the ambiguous `Spin` meaning. Introduce model-aware legacy
adapters and compare resolved decompositions, propagated operator strings,
energy shifts, sweep statistics, thermodynamic estimates, and deterministic
chain seeding.

Preserve TFIM and Rydberg parameter restrictions unless new sign-safe
decompositions are independently proven. Add exact thermal comparisons at tiny
sizes, trace validation, detailed-balance or update-balance checks where
tractable, operator-string growth tests, thread-count determinism, and
statistical regression tests before retiring old paths. The gate passes when
the qslib SSE examples reproduce the standalone crate's validated physics after
explicit conversion and unsupported sign structures fail clearly.

Implementation record (2026-07-19): `qslib-quantum-sse` now uses canonical
`BasisBit` states, explicit per-bond and per-site coefficients, checked TFIM and
Rydberg `WeightedInteraction` adapters, trace-closed propagation, diagonal and
off-diagonal pair updates, and boundary-state updates. Operator cutoffs grow
only during thermalization and fail explicitly if exhausted during measurement.
The shared ADR-0003 BLAKE3 framing derives 32-byte `sse_chain` seeds for
ChaCha20 streams; the historical SplitMix seed remains an explicitly named
legacy adapter. Independent logical chains are reproducible across worker
counts. A provenance-labelled parity fixture, exact one-site Rydberg and
two-site Ising thermal gates, and 23 SSE conformance tests pass with full
workspace tests, Clippy, rustdoc, formatting, and diff checks.

### Milestone 11: configuration, artifacts, and checkpoints

Implement versioned canonical schemas for geometries, models, interactions,
solver settings, SSE runs, exact runs, manifests, summaries, trajectory rows,
and checkpoints. Support human-readable YAML and JSON input, JSON summaries,
and one documented columnar trajectory format selected by ADR. Store convention
schema, software revision, resolved coefficients, dtype, RNG algorithm, seed
derivation, tolerances, backend, checksums, and parameter-layout fingerprints.

Write unknown-field, schema-version, malformed-array, atomic-write, checksum,
round-trip, partial-write recovery, and legacy-load tests first. The gate passes
when a run can be reconstructed from its resolved configuration and artifacts
without relying on a seed to recreate disorder or on native memory layout to
interpret arrays.

### Milestone 12: Python bindings and ncli adoption

Expose only stable, coarse-grained scientific kernels through `qslib-python`.
Use NumPy-compatible arrays with explicit shape, dtype, order, and ownership.
Start with geometry, coupling resolution, Hamiltonian inspection, exact basis
and matrix construction, exact ground states, observables, and TDVP solves.
Avoid per-sample Python callbacks inside Rust hot loops.

Write Python tests before each binding and compare against existing ncli
behavior. Add a pure-Python fallback or preserve the existing ncli path until
parity and packaging are proven. Migrate ncli consumers incrementally behind a
backend protocol. The gate passes when the exact four-site demonstration and
selected TDVP manufactured problems produce equivalent Rust and Python results,
wheel installation works in a clean environment, and ncli can select either
backend explicitly.

### Milestone 13: command line

Implement `qslib` commands for convention and environment inspection, model
validation and resolution, exact ground-state calculation, exact evolution,
supported SSE runs, artifact inspection, and conformance self-tests. Keep the
CLI thin over public APIs. Error messages must name the invalid physical field
and expected convention.

Write CLI integration and snapshot-free semantic output tests first. Support a
machine-readable JSON mode and a readable default. The gate passes when a new
user can run the documented four-site ground-state and tiny SSE examples from
configuration files and inspect complete provenance without writing Rust.

### Milestone 14: documentation and examples

Complete crate-level and public-item rustdoc, a physics-first book, model and
algorithm guides, installation, CLI reference, Python guide, migration guides
for ncli and standalone SSE, reproducibility guidance, limitations, and worked
examples. Generate API documentation and the book into one local site layout
that can later be deployed to GitHub Pages.

All examples must execute in tests or CI. Check internal links, math rendering,
missing docs, and code snippets. The gate passes when a physicist unfamiliar
with Rust can install a local artifact, run the ground-state and SSE tutorials,
interpret every reported quantity, and follow links to detailed APIs.

### Milestone 15: hardening and stability

Benchmark geometry, interaction resolution, matrix construction, matrix-vector
action, observables, TDVP solves, and SSE sweeps. Establish recorded baselines
without making wall-clock timing assertions in unit tests. Add property tests,
targeted fuzzing for parsers and state conversions, Miri-compatible core tests,
large-size overflow tests, dependency license and vulnerability review, feature
combination checks, MSRV checks, and public API semver checks.

Review every public item and feature. Remove accidental exports, unsupported
prototypes, stale compatibility paths whose retirement criteria are satisfied,
and undocumented panics. The gate passes when all supported feature
combinations build, core tests pass under Miri where applicable, fuzz targets
complete their bounded CI run, benchmarks have no unexplained major regression,
and dependency policy is clean.

### Milestone 16: qslib 1.0 release candidate

Run the complete acceptance matrix from a clean checkout. Build optimized Rust
libraries, CLI binaries for locally available targets, Python source and wheel
artifacts, the documentation site, checksums, license bundle, changelog, and
release notes. Install artifacts into clean temporary environments and rerun the
user-visible demonstrations.

Set the version to 1.0.0 only after API review and every required gate passes.
Do not tag, push, publish, sign, or deploy. Present the owner with the exact
artifact paths, checksums, test matrix, known limitations, migration status,
and optional publication commands. The project is complete under this plan
when nothing remains except explicitly authorized external release actions and
post-1.0 roadmap items.

## Concrete commands

Run commands from the qslib project root after Milestone 1 establishes the
workspace. Update this section if accepted ADRs change tooling.

    cargo fmt --check
    cargo clippy --workspace --all-targets --all-features -- -D warnings
    cargo test --workspace --all-features
    RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps --all-features
    cargo test -p qslib-test-support --test conventions

Run Python binding and ncli parity tests after Milestone 12:

    maturin develop -m crates/qslib-python/Cargo.toml
    python -m pytest python/tests
    uv run pytest tests/test_hamiltonians.py tests/test_lattices.py \
      tests/test_exact_ed_vmc.py tests/test_exact_quench.py \
      tests/test_variational_tdvp.py tests/test_row_major_lattice_pairs.py

Run hardening tools after their installation is recorded by Milestone 15:

    cargo deny check
    cargo audit
    cargo semver-checks check-release
    cargo llvm-cov --workspace --all-features --html

The release-candidate demonstration must include commands equivalent to:

    cargo run -p qslib-quantum-cli -- model validate examples/heisenberg_disordered_4.yaml
    cargo run -p qslib-quantum-cli -- exact ground-state examples/heisenberg_disordered_4.yaml
    cargo run -p qslib-quantum-cli -- exact evolve examples/tfim_4.yaml --t-max 0.1
    cargo run -p qslib-quantum-cli -- sse run examples/tfim_thermal_4.yaml
    python examples/python_exact_ground_state.py

Record concise outputs and numerical tolerances in this plan when the examples
are created. Never invent expected energies. Derive and commit them through the
independent conformance fixtures.

## Validation and acceptance

The 1.0 acceptance matrix requires all of the following:

- every mandatory vector in `docs/conventions.md` passes through the public
  facade and the relevant backend;
- all workspace feature combinations selected by the feature policy compile,
  and the documented core-only build excludes heavy optional dependencies;
- every public item is documented and all doctests pass;
- TFIM, Heisenberg and J1-J2, disordered exchange, and Rydberg exact matrices
  match independent tiny-system references;
- full and sector eigensolvers report residuals within quantity-specific
  tolerances and handle degeneracy through invariant subspaces;
- observables state their axes and normalization and match direct evaluation;
- dense and matrix-vector TDVP solvers agree on manufactured problems and use
  the canonical real- and imaginary-time signs;
- resumed accepted-boundary integration is equivalent to uninterrupted
  integration under the declared RNG policy;
- supported SSE models agree with exact thermal results at tiny sizes within
  stated statistical confidence and preserve per-chain determinism across
  thread counts;
- schema round trips retain realized disorder, conventions, dtype, tolerances,
  algorithm, and checksums;
- Rust and Python exact demonstrations agree within declared tolerances;
- CLI errors and JSON output are tested, stable, and physically interpretable;
- CI passes on Linux, macOS, Windows, stable Rust, and MSRV;
- dependency, vulnerability, formatting, lint, documentation, property, fuzz,
  and semver gates pass according to Milestone 15;
- local release artifacts install and run in clean temporary environments;
- no ignored required test, unexplained tolerance relaxation, undocumented
  public API, unresolved accepted ADR, or required `TODO` remains.

Test coverage is evidence, not the objective. Record workspace line and branch
coverage, require full execution of the convention vectors, and investigate
uncovered public scientific branches. Do not add meaningless tests solely to
reach a numeric percentage.

## Idempotence and recovery

Every migration remains additive until parity tests pass. Keep neutral fixtures
independent from both legacy and new implementations. Write generated outputs
to ignored build or temporary directories and use atomic replacement for
durable artifacts. Never overwrite a user data directory during tests.

If dependency selection fails, return to the relevant ADR and evaluate the next
recorded alternative without changing the public scientific contract. If a
sparse solver is unreliable, keep the dense implementation and size guard
working while the tested qslib Lanczos fallback is developed. If Python binding
packaging fails on one platform, keep the Rust release candidate buildable and
record the platform-specific blocker without weakening memory-safety checks.

If a migration test exposes a legacy convention conflict, preserve the legacy
fixture, write the explicit adapter, and keep canonical qslib behavior
unchanged. Never make qslib silently emulate legacy behavior to obtain parity.

At every stopping point, update `Progress`, `Surprises and discoveries`, and the
`Decision log`. A fresh agent must be able to resume by reading only
`AGENTS.md`, `PLANS.md`, this plan, accepted ADRs, and the working tree.

## External authority and owner gates

### Current gate status, updated 2026-07-20 02:05Z

Second resume record: on 2026-07-19 the owner resolved every reduced Milestone 0
gate. Toolchain and dependency downloads plus local commits are authorized.
Remote push, pull, publication, signing, deployment, and destructive migration
remain prohibited.

- Dedicated repository ownership is accepted and implemented. The repository
  uses branch `main` and contains the project agent assets.
- Remote `origin` is configured as `https://github.com/lere01/qslib.git`.
  Network contact, push, and pull are not authorized.
- Registry identifiers are accepted as project brand `qslib`, Rust and Python
  distribution `qslib-quantum`, Rust library target `qslib`, and Python import
  `qslib_quantum`. Publication remains disabled and names are not reserved.
- Repository-scoped Rust 1.85.0, Cargo, Python, documentation, test, audit, and
  packaging dependency downloads are authorized.
- Unsolicited external contributions remain closed. No public code of conduct
  or private conduct-reporting address is adopted at this stage.
- Local initial and milestone commits are authorized. Push and pull are not.
- External publication, push, tag, signing, deployment, and destructive
  migration are explicitly unauthorized by the execution prompt.
- Remote cross-platform CI execution will require later branch/push authority
  or another owner-approved runner. Workflows can be authored and validated
  locally before that gate.
- Milestone 12 changes to the separate ncli ownership unit require later
  coordinated authority now that qslib is a dedicated repository.
- Local semver evidence is available: `cargo-semver-checks 0.42.0` compared
  `HEAD` with baseline commit `2584261` as an assumed patch release and passed
  165 checks with 12 skips. This is an intra-repository API review, not a
  registry release check.
- Local workspace coverage evidence is available from
  `cargo llvm-cov --locked --workspace --all-features --summary-only`: 78.28%
  line coverage, 71.63% region coverage, and no branch counters emitted by the
  installed tool. The report identified the CLI entry point and Python FFI as
  expected uninstrumented paths in this Rust-only run.

Before the autonomous implementation goal begins, the owner should decide or
authorize these items:

1. Resolved 2026-07-19: qslib is a dedicated Git repository because it has
   independent versioning, releases, CI, and consumers.
2. Resolved 2026-07-19: remote is
   `https://github.com/lere01/qslib.git`; push and pull are prohibited.
3. Resolved 2026-07-19: repository-scoped dependency downloads are authorized.
4. Resolved 2026-07-19: local initial and milestone commits are authorized and
   must never include unrelated ncli changes.
5. Resolved: external release actions remain unauthorized.

The agent must stop and request direction for a normative convention change,
license conflict, secret or signing-key requirement, destructive migration of
user data, remote publication, or a scope choice that changes the 1.0 contract.
It should not stop for ordinary implementation choices already bounded by this
plan.

## Recommended autonomous goal prompt

After the owner gates above are resolved, start Codex from the qslib project
root and use this prompt in Goal mode:

    Execute docs/plans/qslib-v1.md under the rules in AGENTS.md and PLANS.md.
    Use the qslib-architect agent as the primary architectural reviewer and the
    qslib skills at every applicable milestone. Work specification-first and
    test-first. Continue milestone by milestone without asking me for routine
    next steps. Keep the ExecPlan current at every stopping point. One primary
    agent owns overlapping edits; delegate only independent exploration,
    audits, and test runs. Preserve unrelated work. Do not publish, tag, push,
    sign, deploy, or perform destructive migration without explicit authority.
    Stop only when every qslib 1.0 acceptance criterion passes or when an owner
    gate in the plan genuinely requires my decision. Return with local release
    artifacts, checksums, test evidence, documentation paths, known limitations,
    and exact optional publication commands.

If the goal pauses because of connectivity, rate limits, or a host restart,
resume the same goal and instruct it to reread this plan and continue from the
`Progress` section. Do not start a second writer on the same checkout.

## Revision note

2026-07-19: Created the initial qslib 1.0 execution plan from the current qslib
conventions, ncli capability inventory, standalone SSE implementation, and the
agreed specification-driven TDD policy.
