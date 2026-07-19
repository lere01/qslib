use qslib_core::{Complex64, FullBasis, Hamiltonian, Pauli, PauliString, SiteCount, SiteId};
use qslib_exact::{
    ExactBasis, GroundState, ThermalSummary, diagonalize_hermitian, evolve, ground_state_sparse,
};

fn two_level() -> Hamiltonian {
    Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(-2.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap(),
        )],
    )
    .unwrap()
}

#[test]
fn exact_basis_preserves_core_order_and_sector_dimension() {
    let full = ExactBasis::full(SiteCount::new(2).unwrap()).unwrap();
    let packed: Vec<_> = full
        .states()
        .iter()
        .map(|state| state.pack().unwrap().words_le()[0])
        .collect();
    assert_eq!(packed, vec![0, 1, 2, 3]);
    let sector = ExactBasis::fixed_weight(SiteCount::new(4).unwrap(), 2).unwrap();
    assert_eq!(sector.dimension(), 6);
    assert!(
        sector
            .states()
            .iter()
            .all(|state| state.hamming_weight() == 2)
    );
}

#[test]
fn dense_and_csr_matrices_match_direct_hamiltonian_action() {
    let h = two_level();
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let dense = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let csr = qslib_exact::CsrMatrix::from_hamiltonian(&h, &basis).unwrap();
    assert_eq!(dense.dimension(), 2);
    assert_eq!(
        dense.as_slice(),
        &[
            Complex64::new(-2.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(2.0, 0.0)
        ]
    );
    let vector = vec![Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)];
    assert_eq!(dense.apply(&vector).unwrap(), csr.apply(&vector).unwrap());
    assert_eq!(
        dense.apply(&vector).unwrap(),
        vec![Complex64::new(-2.0, 0.0), Complex64::new(4.0, 0.0)]
    );
}

#[test]
fn hermitian_diagonalization_reports_spectrum_and_ground_state_residual() {
    let h = two_level();
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    assert_eq!(spectrum.values(), &[-2.0, 2.0]);
    let ground = GroundState::from_spectrum(&spectrum).unwrap();
    assert_eq!(ground.energy(), -2.0);
    assert!(ground.residual() < 1.0e-12);
}

#[test]
fn complex_hermitian_spectrum_and_nonhermitian_diagnostic_are_explicit() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::Y)]).unwrap(),
        )],
    )
    .unwrap();
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    assert_eq!(spectrum.values(), &[-1.0, 1.0]);
    assert!(
        spectrum
            .residuals()
            .iter()
            .all(|residual| *residual < 1.0e-12)
    );
    let nonhermitian = qslib_exact::DenseMatrix::new(
        2,
        vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(0.0, 0.0),
        ],
    )
    .unwrap();
    assert!(diagonalize_hermitian(&nonhermitian).is_err());
}

#[test]
fn off_diagonal_hermitian_eigenvectors_have_small_residuals() {
    let matrix = qslib_exact::DenseMatrix::new(
        2,
        vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
        ],
    )
    .unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    assert!(
        spectrum
            .residuals()
            .iter()
            .all(|residual| *residual < 1.0e-12)
    );
}

#[test]
fn thermal_sum_and_unitary_or_imaginary_evolution_match_analytic_values() {
    let h = two_level();
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    let thermal = ThermalSummary::from_spectrum(&spectrum, 1.0).unwrap();
    assert!((thermal.partition_function() - 2.0 * 2.0_f64.cosh()).abs() < 1.0e-12);
    assert!((thermal.energy() + 2.0_f64.tanh() * 2.0).abs() < 1.0e-12);
    let initial = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
    let evolved = evolve(&matrix, &initial, std::f64::consts::PI / 4.0, false).unwrap();
    assert!((evolved[0] - Complex64::new(0.0, 1.0)).norm() < 1.0e-12);
    let imaginary = evolve(&matrix, &initial, 0.5, true).unwrap();
    assert!((imaginary.iter().map(|value| value.norm_sqr()).sum::<f64>() - 1.0).abs() < 1.0e-12);
}

#[test]
fn full_basis_iterator_remains_the_independent_order_oracle() {
    let states: Vec<_> = FullBasis::new(SiteCount::new(2).unwrap())
        .unwrap()
        .collect();
    let exact = ExactBasis::full(SiteCount::new(2).unwrap()).unwrap();
    let exact_packed: Vec<_> = exact
        .states()
        .iter()
        .map(|state| state.pack().unwrap())
        .collect();
    assert_eq!(states, exact_packed);
}

#[test]
fn heterogeneous_four_site_ground_state_matches_independent_classical_bound() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(4).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![
            (
                Complex64::new(-1.0, 0.0),
                PauliString::new(vec![(SiteId::new(0), Pauli::Z), (SiteId::new(1), Pauli::Z)])
                    .unwrap(),
            ),
            (
                Complex64::new(-2.0, 0.0),
                PauliString::new(vec![(SiteId::new(1), Pauli::Z), (SiteId::new(2), Pauli::Z)])
                    .unwrap(),
            ),
            (
                Complex64::new(0.5, 0.0),
                PauliString::new(vec![(SiteId::new(2), Pauli::Z), (SiteId::new(3), Pauli::Z)])
                    .unwrap(),
            ),
        ],
    )
    .unwrap();
    let basis = ExactBasis::full(SiteCount::new(4).unwrap()).unwrap();
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let ground = GroundState::from_spectrum(&diagonalize_hermitian(&matrix).unwrap()).unwrap();
    assert!((ground.energy() + 3.5).abs() < 1.0e-12);
    assert!(ground.residual() < 1.0e-11);
    let sparse = qslib_exact::CsrMatrix::from_hamiltonian(&h, &basis).unwrap();
    let sparse_ground = ground_state_sparse(&sparse, 1.0e-10, 16).unwrap();
    assert!((sparse_ground.energy() - ground.energy()).abs() < 1.0e-10);
    assert!(sparse_ground.residual() < 1.0e-10);
}
