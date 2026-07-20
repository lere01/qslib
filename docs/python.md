# Python bindings

`qslib-quantum-python` builds the `qslib_quantum` import package with Maturin
and the PyO3 `abi3-py312` interface. The binding accepts NumPy arrays only for
the duration of one call and returns newly allocated C-contiguous arrays; it
does not retain Python buffers or expose mutable Rust handles.

The initial stable kernel surface is intentionally coarse grained:

- `row_major_site_ids(lx, ly)` returns canonical `uint32` site identifiers in
  `(Ly, Lx)` row-major order.
- `resolve_ising_interactions(couplings)` validates a symmetric dense matrix
  and returns integer endpoint pairs plus their resolved `float64` coefficients.
- `basis_states(site_count, weight=None)` returns the canonical exact basis.
- `tfim_matrix` and `tfim_ground_state` expose exact TFIM matrices and residual-
  checked ground states.
- `heisenberg_matrix` and `rydberg_matrix` expose the corresponding canonical
  exact matrices for pair-dependent couplings.
- `tfim_observables` reports energy, energy density, energy variance, and an
  axis-labelled magnetization total and density, with an optional two-site
  correlation. Energy is a total and all densities divide by the number of
  physical sites.
- `tdvp_estimate` and `tdvp_solve` expose owned sufficient statistics and a
  checked dense QGT solve with explicit regularization shift and diagnostics.

Inputs are validated for shape, finite values, symmetry, site order, and model
basis. C-order, Fortran-order, sliced, and read-only NumPy views are copied or
indexed safely. Wrong dtypes and malformed ranks are rejected by the Python
boundary. Errors are exposed as `InputError` or `NumericalError`, both derived
from `QslibError`.

The 1.0 binding keeps the Python GIL while exact dense kernels run. This is a
deliberate small-system policy: calls copy their inputs before entering Rust,
retain no Python-owned buffers, and return promptly bounded results. A future
release may detach the GIL for larger kernels after adding an explicit
threading contract and benchmark evidence.

Build and test a local wheel without publishing it:

```text
maturin build --manifest-path crates/qslib-python/Cargo.toml --out dist
python -m venv /tmp/qslib-venv
/tmp/qslib-venv/bin/pip install dist/qslib_quantum-*.whl pytest
/tmp/qslib-venv/bin/pytest crates/qslib-python/tests/python_contract.py
```

For an authorized PyPI release, install the same binding with:

```text
python -m pip install qslib-quantum==1.0.0
```

The package name and import name are intentionally different: distribution
`qslib-quantum` installs the `qslib_quantum` module. Pin the version for
reproducible research and use the release checksum manifest when installing
from a downloaded wheel instead of PyPI.

The wheel is ABI-compatible with Python 3.12 and later on the platform where
it is built. Publishing, signing, and registry uploads remain owner-gated.

This is the qslib-local backend surface. ncli backend selection and parity
adapters are intentionally deferred until the separately owned parent project
grants coordinated migration authority.
