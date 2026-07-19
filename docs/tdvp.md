# Variational and TDVP kernels

`qslib-quantum-variational` accepts the numerical products of a wavefunction
frontend without depending on an autodiff framework. A caller supplies one
complex local energy and one complex log derivative per sample and parameter,
in row-major sample order. Non-negative normalized or unnormalized sample
weights are accepted; the estimator normalizes them internally.

For real parameters, the estimator forms

```text
S_kl = Re E[conj(delta O_k) delta O_l]
F_k  = E[conj(delta O_k) delta E_loc]
```

and uses `Im(F)` for real time or `-Re(F)` for imaginary time. Energy mean and
variance are retained alongside the dense real QGT. `DenseQgt::matvec` and
`solve_cg` cover dense and caller-supplied matrix-free solves.
`qgt_vector_product_stream` consumes derivative chunks directly, so a frontend
can stream rows without constructing a sample-space identity.

`Regularization::FixedTikhonov` adds a documented diagonal shift,
`Regularization::Gcv` selects a positive shift from a supplied grid, and
`Regularization::SpectralCutoff` removes modes below a relative eigenvalue
threshold. Solve results report clipping, QGT metric norm, projected residual,
and normalized residual. `ParameterLayout` records deterministic names,
shapes, offsets, and a stable fingerprint for checkpoint compatibility.

The Rust layer does not compute neural-network derivatives. Python and other
frontends may adapt their own derivative engines to these checked arrays at a
coarse boundary.

For local energy, `local_energy_from_ratios` expects row-oriented matrix
elements `(H_{b b'}, psi(b')/psi(b))`, including the diagonal separately. A
`Hamiltonian::apply` is column-oriented and returns `c * P_ba`; do not
conjugate that whole coefficient for complex `c`. Prefer a row/local-energy
API, or conjugate only the Pauli matrix element when constructing the pair.
