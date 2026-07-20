# Limitations and numerical scope

qslib 1.0 is deliberately a validated finite-system foundation. The exact
backend enumerates full or fixed-weight bases and therefore scales
exponentially with site count. Its dense reference eigensolver is intended for
small matrices and reports residuals so callers can reject an under-resolved
result.

The variational layer supplies statistics, dense QGT solves, regularization,
and transactional integration. It does not compute neural-network derivatives
or prescribe an optimizer. The Python binding exposes coarse kernels and
copies foreign arrays for each call. It retains the GIL during the bounded
small-system kernels in the 1.0 contract.

The SSE backend supports the sign-safe TFIM and Rydberg decompositions that are
covered by conformance tests. Statistical estimates require user-selected
thermalization, chain count, autocorrelation analysis, and confidence
criteria. The CLI's tiny SSE example is a smoke test, not a production error
bar.

The ncli backend-selection adapter remains outside this repository's ownership
boundary. qslib preserves explicit adapters and does not silently reinterpret
legacy site order, spin labels, coupling signs, or basis conventions.
