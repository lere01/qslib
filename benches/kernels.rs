//! Lightweight benchmark baselines for the kernels most likely to regress.
//!
//! This is a reporting benchmark, not a wall-clock unit test. Run with
//! `cargo bench --all-features --bench kernels` and retain the environment,
//! compiler, and system details alongside the printed measurements.

use qslib::exact::{DenseMatrix, ExactBasis, expectation};
use qslib::sse::{BasisSseState, LocalSseModel, Operator, SimulationConfig, run_parallel_chains};
use qslib::variational::{TDVPMode, TDVPSolveOptions, estimate_tdvp};
use qslib::{
    Boundary, Complex64, DenseCouplings, InteractionChannel, InteractionTable, RectangularGeometry,
    SimulationBasis, SiteCount,
};
use std::hint::black_box;
use std::time::Instant;

fn tfim_matrix(sites: usize) -> DenseMatrix {
    let mut values = vec![0.0; sites * sites];
    for site in 0..sites {
        let next = (site + 1) % sites;
        values[site * sites + next] = 1.0;
        values[next * sites + site] = 1.0;
    }
    let dense = DenseCouplings::new(SiteCount::new(sites).unwrap(), values).unwrap();
    let interactions = dense.to_interactions(InteractionChannel::IsingZZ).unwrap();
    let table = InteractionTable::new(
        dense.site_count(),
        interactions
            .into_iter()
            .map(|term| (term.bond(), term.channel().clone(), term.coefficient()))
            .collect(),
    )
    .unwrap();
    let model = qslib::tfim(&table, &vec![0.5; sites], SimulationBasis::Z).unwrap();
    let basis = ExactBasis::full(SiteCount::new(sites).unwrap()).unwrap();
    DenseMatrix::from_hamiltonian(model.hamiltonian(), &basis).unwrap()
}

fn main() {
    let geometry_start = Instant::now();
    for _ in 0..128 {
        black_box(
            RectangularGeometry::new(4, 4, Boundary::Periodic, Boundary::Periodic)
                .unwrap()
                .bonds(qslib::BondMultiplicity::Simple)
                .unwrap(),
        );
    }
    println!(
        "periodic geometry and bonds (128 runs): {:?}",
        geometry_start.elapsed()
    );

    let interaction_start = Instant::now();
    for _ in 0..128 {
        let dense = DenseCouplings::new(SiteCount::new(4).unwrap(), vec![0.0; 16]).unwrap();
        black_box(dense.to_interactions(InteractionChannel::IsingZZ).unwrap());
    }
    println!(
        "interaction resolution (128 runs): {:?}",
        interaction_start.elapsed()
    );

    let construction_start = Instant::now();
    for _ in 0..32 {
        black_box(tfim_matrix(4));
    }
    println!(
        "dense matrix construction (32 runs): {:?}",
        construction_start.elapsed()
    );

    let matrix = tfim_matrix(4);
    let vector = vec![Complex64::new(1.0, 0.0); matrix.dimension()];
    let start = Instant::now();
    for _ in 0..16 {
        black_box(matrix.apply(black_box(&vector)).unwrap());
    }
    println!("matrix-vector action (16 runs): {:?}", start.elapsed());

    let start = Instant::now();
    for _ in 0..2 {
        black_box(qslib::exact::diagonalize_hermitian(black_box(&matrix)).unwrap());
    }
    println!("dense diagonalization (2 runs): {:?}", start.elapsed());

    let start = Instant::now();
    for _ in 0..32 {
        black_box(expectation(black_box(&matrix), black_box(&vector)).unwrap());
    }
    println!("exact expectation (32 runs): {:?}", start.elapsed());

    let weights = vec![1.0; 32];
    let energies = vec![Complex64::new(1.0, 0.5); 32];
    let derivatives = vec![Complex64::new(0.25, 0.1); 32 * 8];
    let stats = estimate_tdvp(&weights, &energies, &derivatives, 8, TDVPMode::RealTime).unwrap();
    let start = Instant::now();
    for _ in 0..128 {
        black_box(
            estimate_tdvp(
                black_box(&weights),
                black_box(&energies),
                black_box(&derivatives),
                8,
                TDVPMode::RealTime,
            )
            .unwrap(),
        );
    }
    println!("TDVP statistics (128 runs): {:?}", start.elapsed());

    let start = Instant::now();
    for _ in 0..32 {
        black_box(stats.solve(TDVPSolveOptions::default()).unwrap());
    }
    println!("TDVP solve (32 runs): {:?}", start.elapsed());

    let start = Instant::now();
    let bonds = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0)];
    let model = LocalSseModel::tfim_weighted(4, &bonds, &[0.5; 4]).unwrap();
    let state = BasisSseState::new(
        vec![qslib::BasisBit::Zero; 4],
        vec![Operator::identity(); 32],
    )
    .unwrap();
    for _ in 0..4 {
        black_box(
            run_parallel_chains(
                model.clone(),
                state.clone(),
                1.0,
                SimulationConfig {
                    thermalization_sweeps: 2,
                    measurement_sweeps: 4,
                    sweeps_per_measurement: 1,
                },
                17,
                1,
                1,
            )
            .unwrap(),
        );
    }
    println!("SSE sweeps (4 short runs): {:?}", start.elapsed());
}
