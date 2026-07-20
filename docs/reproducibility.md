# Reproducibility

qslib treats reproducibility as a scientific result, not as a promise that a
seed alone recreates every run. A durable record should retain:

- the resolved model family, every pair coefficient and onsite parameter;
- canonical row-major site order, boundaries, physical axes, and simulation
  basis;
- algorithm, integrator, regularizer, cutoff, tolerance, and accepted-step
  controls;
- random algorithm, seed derivation version, logical chain identifiers, and
  checkpoint RNG position;
- software revision, schema version, dtype, parameter-layout fingerprint, and
  checksums for arrays and columnar parts.

For disordered models, persist the realized interaction table. Regenerating it
from a seed can hide changes in ordering, periodic image multiplicity, or a
future random-stream implementation. For stochastic SSE estimates, separate
finite-sample uncertainty, autocorrelation, and finite-disorder ensemble error.

The IO layer provides strict JSON/YAML configuration, typed checkpoints, named
NPY arrays, and append-only Parquet trajectories. The independent verifier in
`tools/verify_io_artifacts.py` is useful in a clean Python environment. The
CLI and Python examples in this repository are small deterministic contract
fixtures; they are not substitutes for production thermalization or error
analysis.
