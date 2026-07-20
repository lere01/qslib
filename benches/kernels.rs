//! Lightweight benchmark baselines for the kernels most likely to regress.
//!
//! This is a reporting benchmark, not a wall-clock unit test. Run with
//! `cargo bench --all-features --bench kernels` and retain the environment,
//! compiler, and system details alongside the printed measurements.

use qslib::exact::{DenseMatrix, ExactBasis, expectation};
use qslib::variational::{TDVPMode, estimate_tdvp};
use qslib::{
    Complex64, DenseCouplings, InteractionChannel, InteractionTable, SimulationBasis, SiteCount,
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
}
