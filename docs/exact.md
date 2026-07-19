# Exact calculations

`qslib-quantum-exact` is the reference small-system backend. It consumes the
canonical `BasisState`, `Hamiltonian`, site order, and bit convention from
`qslib-core`; it does not introduce a second basis encoding.

## Basis and matrix contracts

`ExactBasis::full` enumerates packed little-endian states in increasing integer
order. `ExactBasis::fixed_weight` uses the same order while retaining exactly
the requested Hamming weight. `DenseMatrix::from_hamiltonian` stores rows and
columns in that basis order and uses the conventional action
`y[row] = sum_column H[row,column] x[column]`. `CsrMatrix` is a deterministic
compressed representation of the same matrix and is checked against dense
matrix-vector action in the test suite. CSR assembly resolves Hamiltonian
actions directly into row storage and rejects a non-conserving Hamiltonian when
the supplied basis is a fixed sector.

The builder preserves physical constants and all pair-dependent coefficients.
It reports missing connected states, dimension mismatches, non-Hermitian input,
and non-convergent numerical work as typed errors.

## Spectra, thermal sums, and evolution

`diagonalize_hermitian` is a deterministic reference Hermitian eigensolver. It
returns ascending eigenvalues, normalized complex eigenvectors, and the norm of
`H|v>-lambda|v>` for every pair. `GroundState` selects the lowest eigenvalue;
degenerate vectors must be compared through their invariant subspace rather
than by vector identity. `ground_state_sparse` uses deterministic full
reorthogonalization with basis-vector restarts and the same residual contract
for CSR matrices. `Eigensystem::projector` compares degenerate invariant
subspaces without comparing arbitrary eigenvector representatives.

`ThermalSummary` evaluates `Z`, stable `log Z`, the canonical mean energy, and
`beta^2 Var(H)` from the exact spectrum using an energy-shifted log-sum-exp
calculation. `evolve` uses the same spectral data
for unitary `exp(-i H t)` evolution and normalized imaginary-time
`exp(-H tau)` evolution. Imaginary time and inverse temperature are
non-negative; real time may be signed. All parameters are finite and use the
natural units defined by the project conventions.
If `Z` exceeds the finite `f64` range, `partition_function_overflowed` is true
and `log_partition_function` remains the authoritative finite-scale result.

The backend is intended for validation and small systems. It is not a claim of
production-scale exact diagonalization; later milestones may add a maintained
sparse solver or a qslib-owned Lanczos implementation behind a separate,
diagnosable API.

The exact observable surface also includes normalized Pauli-matrix expectation
and variance, Shannon entropy, pure-state bipartite entropy, axis-labelled
magnetization, raw and connected correlations, position- and momentum-bound
structure factors, weighted sublattice moments, QFI normalization, and thermal
heat-capacity density. These APIs bind site count, physical axis, normalization,
and connected-correlation choice explicitly rather than accepting pre-aggregated
scalars without provenance.
