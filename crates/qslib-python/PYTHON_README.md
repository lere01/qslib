# qslib-quantum Python bindings

`qslib-quantum` provides the `qslib_quantum` Python module for small-system
quantum simulation kernels. It exposes exact TFIM, Heisenberg, and Rydberg
matrices, ground states, convention-labelled observables, and checked TDVP
statistics through NumPy arrays.

Install an authorized release with:

```text
python -m pip install qslib-quantum
```

The binding requires Python 3.12 or newer. Inputs use qslib's canonical
row-major site order and little-endian basis convention. See the main project
documentation at <https://lere01.github.io/qslib/> for physical definitions,
sign conventions, limits, and examples.

This package is distributed under the Apache-2.0 license.
