# Migration from ncli

qslib is the canonical scientific layer for future ncli integration. The
migration boundary is intentionally explicit because ncli and qslib may have
different site ordering, basis labels, array ownership, and optimizer state.

The qslib side is ready with:

- canonical row-major geometry and little-endian packed states;
- typed weighted interactions that preserve pair identity and realized
  disorder;
- exact, TDVP, SSE, and versioned artifact kernels;
- an ABI3 `qslib_quantum` Python surface for coarse matrices, observables,
  ground states, and TDVP statistics.

The parent ncli repository must add and test the backend protocol before a
consumer migration is complete. That protocol should select Rust or Python
explicitly, preserve the existing ncli path until parity is proven, and record
the adapter's site-order, basis, dtype, and normalization conversions. qslib
does not modify the separately owned ncli tree or claim parity that has not run.
