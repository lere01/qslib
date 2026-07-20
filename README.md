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

The current Cargo package is an architecture scaffold. It has an explicit,
documented `qslib` library target and still retains the original `Hello, world!`
binary until Milestone 1 replaces the template with the approved workspace.
No scientific API is implemented yet.

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
