# Kernel benchmarks

The benchmark target is a small regression smoke suite, not a portable
performance claim. Run it with the pinned minimum toolchain and retain the
machine, compiler, feature set, and workload alongside any comparison:

```text
cargo +1.85.0 bench --locked --all-features --bench kernels
```

The baseline below was recorded on the local development host after the M15
kernel expansion. Times are total wall-clock times for the stated number of
short repetitions and are intended only for same-host comparisons.

```text
periodic geometry and bonds (128 runs): 131.375us
interaction resolution (128 runs): 21.083us
dense matrix construction (32 runs): 632.375us
matrix-vector action (16 runs): 4.583us
dense diagonalization (2 runs): 652.541us
exact expectation (32 runs): 7.125us
TDVP statistics (128 runs): 264.75us
TDVP solve (32 runs): 46.166us
SSE sweeps (4 short runs): 454.083us
```

The target covers geometry and interaction resolution, dense matrix
construction and action, exact diagonalization and expectation, TDVP
statistics and solve, and a short SSE sweep. A future performance gate should
compare like-for-like workloads and report compiler, CPU, feature, and thread
settings rather than treating these values as universal thresholds.
