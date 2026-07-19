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
matrix-vector action in the test suite.

The builder preserves physical constants and all pair-dependent coefficients.
It reports missing connected states, dimension mismatches, non-Hermitian input,
and non-convergent numerical work as typed errors.

## Spectra, thermal sums, and evolution

`diagonalize_hermitian` is a deterministic reference Hermitian eigensolver. It
returns ascending eigenvalues, normalized complex eigenvectors, and the norm of
`H|v>-lambda|v>` for every pair. `GroundState` selects the lowest eigenvalue;
degenerate vectors must be compared through their invariant subspace rather
than by vector identity. `ground_state_sparse` uses deterministic full
reorthogonalization and the same residual contract for CSR matrices.

`ThermalSummary` evaluates `Z`, the canonical mean energy, and
`beta^2 Var(H)` from the exact spectrum. `evolve` uses the same spectral data
for unitary `exp(-i H t)` evolution and normalized imaginary-time
`exp(-H tau)` evolution. Time and inverse temperature are non-negative and
finite, in the natural units defined by the project conventions.

The backend is intended for validation and small systems. It is not a claim of
production-scale exact diagonalization; later milestones may add a maintained
sparse solver or a qslib-owned Lanczos implementation behind a separate,
diagnosable API.
