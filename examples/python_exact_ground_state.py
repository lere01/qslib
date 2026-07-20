"""Compare the Python binding's exact four-site ground-state contract."""

import numpy as np
import qslib_quantum as qslib


couplings = np.array(
    [
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 0.0, 1.0, 0.0],
    ],
    dtype=np.float64,
)
energy, vector, residual = qslib.tfim_ground_state(
    couplings, np.full(4, 0.5), basis="z"
)
print(f"energy = {energy:.12g}")
print(f"residual = {residual:.3e}")
print(f"norm = {np.vdot(vector, vector).real:.12g}")
