# qslib command line

The `qslib` command is a thin interface over the validated Rust kernels. It
lets a physicist inspect conventions, validate a Hamiltonian, calculate a tiny
exact ground state or evolution, run a small TFIM SSE demonstration, and
inspect an artifact directory without writing Rust.

Build it locally from the repository root:

```text
cargo build --release -p qslib-quantum-cli
target/release/qslib inspect conventions
```

Use `--json` for scripts and notebooks. Without it, the command prints one
labelled field per line. Errors are written to stderr and return exit status 2.
They identify the physical field or convention that failed validation.

## Configuration conventions

Configuration files are YAML or JSON. Dense coupling matrices are square,
symmetric, zero-diagonal arrays in canonical row-major site order. Each upper
triangle entry is one resolved pair coefficient, including zero coefficients.
`basis` names the simulation basis (`z`, `x`, or `y` where the model supports
it), while physical observable axes remain separate concepts.

For TFIM,

\[
H=-\sum_{i<j}J_{ij}\sigma_i^z\sigma_j^z-sum_i h_i\sigma_i^x
\]

when `basis: z`. `fields` therefore contains one transverse field per site.
Heisenberg configurations use the same pair matrix and represent
`J_ij S_i . S_j` with `S = sigma/2`. Rydberg configurations use `omega` and
`detuning` arrays in addition to the pair matrix. J1-J2 configurations use
`lx`, `ly`, `j1`, and `j2`, with optional `open` or `periodic` boundaries.

## Commands

```text
qslib inspect conventions [--json]
qslib inspect environment [--json]
qslib model validate CONFIG [--json]
qslib exact ground-state CONFIG [--json]
qslib exact evolve CONFIG --t-max TIME [--imaginary] [--json]
qslib sse run CONFIG [--json]
qslib artifacts inspect PATH [--json]
qslib conformance self-test [--json]
```

`exact ground-state` reports the total energy, exact-basis dimension, and the
residual norm \(\|H|\psi\rangle-E|\psi\rangle\|\). Degenerate ground spaces
should be interpreted through the energy and residual, not a phase-sensitive
individual eigenvector. `exact evolve` starts from canonical basis state zero;
real time uses \(e^{-iHt}\), and imaginary time uses the normalized
\(e^{-H\tau}\) convention.

The SSE command is a small-system TFIM demonstration. It reports the mean
energy per site across independent logical chains. It does not claim an
autocorrelation-corrected confidence interval; production studies must choose
thermalization, sweep, chain, and uncertainty controls for the physical
question.

## Worked examples

```text
cargo run -p qslib-quantum-cli -- model validate examples/heisenberg_disordered_4.yaml
cargo run -p qslib-quantum-cli -- exact ground-state examples/tfim_4.yaml
cargo run -p qslib-quantum-cli -- exact evolve examples/tfim_4.yaml --t-max 0.1
cargo run -p qslib-quantum-cli -- sse run examples/tfim_thermal_4.yaml
```

The examples are executed by the CLI integration tests. The Python equivalent
uses the `qslib_quantum` import documented in [`python.md`](python.md).
