# Symmetry actions and sectors

qslib stores a site permutation as `source_for_destination`. Its action is

\[
(g\cdot b)_j=b_{p_g(j)},
\]

so it has the same gather direction as a row-major array read. `Permutation`
validates bijectivity, composes and inverts actions, and can act on dense
binary states or map Pauli supports. A separately named adapter is required
when converting to a destination-for-source convention.

`translation` and `translation_group` derive coordinate actions from the
declared geometry and boundaries. Rectangle and square point-group builders
return validated finite groups. Open-boundary translations that leave the
finite domain are rejected; periodic translations wrap explicitly.

`SpinInversion` applies \(b_i\mapsto1-b_i\). It is an available action, not an
automatic claim that a particular Hamiltonian is invariant. Use
`is_interaction_symmetry` to verify that every resolved weighted interaction,
including its channel, name, bond, and coefficient, is preserved.
`validate_model_symmetry` checks the complete resolved Pauli Hamiltonian,
including onsite fields, and reports the first missing mapped term;
`validate_spin_inversion` performs the corresponding commutator check.

`DiagonalGauge` implements \(\phi(b)=\pi\sum_i a_i b_i\) and applies the
resulting complex phase to amplitudes. Integer coefficients are classified as
sign-only. `sublattice_gauge` supplies the checkerboard gauge for square
geometries and rejects unsupported triangular embeddings.

For a finite group and a one-dimensional character, `FiniteGroup::project_amplitudes`
implements

\[
(P_\chi\psi)(b)=|G|^{-1}\sum_g\chi(g)^*\psi(g^{-1}\cdot b).
\]

Orbit methods return deterministic packed-state representatives. Projection
and orbit utilities operate on explicit dense amplitudes and are intended as
small-system reference primitives before sector-specific exact solvers.
