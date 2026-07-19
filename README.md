# qslib

`qslib` is a quantum simulation library written in Rust. It is intended to
provide shared scientific definitions and numerical building blocks for exact
diagonalization, variational Monte Carlo, TDVP, stochastic series expansion,
symmetry analysis, homogeneous and disordered lattice models, observables, and
future simulation programs.

The project is at the convention-first stage. Its first public contract is the
[scientific convention specification](docs/conventions.md). Implementations
must be tested against that specification rather than inheriting implicit
conventions from an earlier program.

The current dependency direction and subsystem boundaries are documented in
the [architecture overview](docs/architecture/README.md). Consequential design
choices are recorded as [architectural decision records](docs/decisions/README.md).
The living [qslib 1.0 execution plan](docs/plans/qslib-v1.md), governed by
[`PLANS.md`](PLANS.md), defines the implementation milestones and completion
criteria for autonomous development.

The current Cargo workspace exposes the `qslib` facade, while capability crates
isolate core, exact, variational, SSE, IO, CLI, Python, and test-support
responsibilities. The core crate now implements the checked basis, geometry,
weighted-interaction, and disorder foundations described by the completed
Milestones 2 and 3. Model assembly and numerical algorithms remain later
milestones.

## Workspace and features

The default facade contains only `qslib-core`. Optional additive features are
`exact`, `variational`, `sse`, and `io`; `full` enables all four Rust library
capabilities. CLI and Python are separate packages, so core users do not compile
interface or heavy-backend dependencies. Cargo package names use the
`qslib-quantum-` prefix, while Rust library targets use concise names such as
`qslib_core` and `qslib_exact`.

The language-neutral fixtures under
[`fixtures/conformance/v1/`](fixtures/conformance/v1/README.md) record the
small analytic systems that later implementations must reproduce. They include
explicit basis order, resolved coefficients, matrix layout, oracle provenance,
comparison policy, and a checksummed manifest. Passing the current harness
means the evidence is internally valid. It does not claim that scientific
algorithms have already been implemented.

The local CI contract is defined in
[`.github/workflows/ci.yml`](.github/workflows/ci.yml). Remote workflow
execution remains pending owner-authorized push activity.

For a physicist-first introduction to site order, boundaries, pair-dependent
couplings, and disorder provenance, see the
[geometry and interactions guide](docs/geometry-interactions.md).

## Project policies

- [qslib 1.0 execution plan](docs/plans/qslib-v1.md)
- [Scientific conventions](docs/conventions.md)
- [Architecture](docs/architecture/README.md)
- [Architectural decisions](docs/decisions/README.md)
- [Rust toolchain and MSRV](docs/toolchain-policy.md)
- [Contributing](CONTRIBUTING.md)
- [Code-of-conduct decision](docs/governance/code-of-conduct-decision.md)
- [Security](SECURITY.md)
- [Changelog](CHANGELOG.md)
- [Apache-2.0 license](LICENSE)

External publication is disabled. The `qslib` crates.io and PyPI distribution
names belong to an unrelated project. The approved Rust and Python distribution
name is `qslib-quantum`; Rust imports use `qslib`, while Python imports use the
collision-safe `qslib_quantum`. Naming does not authorize publication. See
[ADR-0009](docs/decisions/0009-registry-package-names.md).

The project does not currently accept unsolicited external contributions. See
the [contribution policy](CONTRIBUTING.md) before proposing changes.
