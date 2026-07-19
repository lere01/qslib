# Operators and model Hamiltonians

qslib represents a physical Hamiltonian as

\[
H=cI+\sum_a c_a O_a,
\]

where `c` is a physical energy constant and every Pauli-string term has its
own resolved coefficient. Algorithmic shifts are not folded into `c`.

## Basis action

Bits use the simulation axis as their diagonal Pauli operator. In a z basis,
`Z` multiplies a state by \(1-2b_i\), while `X` flips site `i` and `Y` flips it
with phase `i` on bit zero and `-i` on bit one. `PauliString::new` accepts
canonical distinct support. `PauliString::product` is the explicit
ordered-product API and reduces repeated factors, including `XY=iZ` and
`YX=-iZ`.

## TFIM

The physical convention is \(H=-\sum_{ij}J_{ij}Z_iZ_j-\sum_i h_iX_i\).
`tfim` requires `InteractionChannel::IsingZZ`. In simulation basis z it emits
`ZZ` bond terms and `X` field terms. In basis x it emits `XX` bond terms and
diagonal `Z` field terms. Positive `J_ij` is ferromagnetic under this sign
convention; signed heterogeneous values remain resolved term by term. The y
basis is rejected explicitly because this first conversion surface does not
yet provide its complex phase convention.

Each builder returns a `ResolvedModel`. It dereferences to the operator
Hamiltonian for matrix and connected-state APIs, and also exposes `family`,
`basis`, the complete weighted interaction identities (including names and
zero coefficients), and a typed `ModelSpecification` containing named site
fields or shell vectors. This prevents a numerical operator from losing the
physical inputs that produced it.

## Heisenberg and J1-J2

The isotropic exchange is
\(H=\sum_{ij}J_{ij}(X_iX_j+Y_iY_j+Z_iZ_j)/4\). The builder requires
`HeisenbergExchange` and supports x, y, and z simulation bases because the
isotropic sum is invariant under a common axis rotation. Positive `J_ij` is
antiferromagnetic and negative values are allowed. The `j1j2` convenience
constructor selects axial nearest-neighbour and diagonal next-nearest-neighbour
shells, names them `j1` and `j2`, and delegates to the same resolved term path.
`j1j2_disordered` accepts one signed or zero coefficient per resolved bond in
each shell, validates both lengths, and retains the shell names even when
periodic geometry makes two physical contributions share an endpoint pair.

## Rydberg

The z-basis builder uses
\(H=-\sum_i\Omega_iX_i/2-\sum_i\Delta_i n_i+\sum_{i<j}V_{ij}n_in_j\),
with `n_i=b_i`. It requires `RydbergDensityDensity`, validates per-site drive
and detuning, and accepts a symmetric finite zero-diagonal coupling matrix.
Its density expansion retains the physical constant and Z terms. Other bases
return an explicit unsupported-basis error rather than silently changing the
occupation convention.

## Application and local energy

`Hamiltonian::apply` returns unique connected states with algebraically combined
coefficients in canonical packed-mask order. `Hamiltonian::local_energy`
computes \(\sum_b H_{ab}\psi_b/\psi_a\) from an explicit amplitude table. Since
`apply` is column-oriented and returns \(H_{ba}\), local energy conjugates the
Pauli matrix element when evaluating the row element while preserving the
physical complex coefficient. It combines and prunes exactly cancelling row
transitions before requesting connected amplitudes, then rejects a missing
uncancelled amplitude or zero reference amplitude. Physical terms and constants
remain separate from future solver shifts.
