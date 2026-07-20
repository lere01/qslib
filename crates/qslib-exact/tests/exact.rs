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

fn assert_real_slices_close(actual: &[f64], expected: &[f64], absolute: f64, relative: f64) {
    assert_eq!(actual.len(), expected.len());
    for (index, (actual_value, expected_value)) in actual.iter().zip(expected).enumerate() {
        let scale = actual_value.abs().max(expected_value.abs());
        let tolerance = absolute + relative * scale;
        assert!(
            (actual_value - expected_value).abs() <= tolerance,
            "spectrum mismatch at index {index}: actual={actual_value:.17e}, \
             expected={expected_value:.17e}, tolerance={tolerance:.3e}"
        );
    }
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
    assert_real_slices_close(spectrum.values(), &[-1.0, 1.0], 1.0e-12, 1.0e-12);
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
fn degenerate_zero_spectrum_is_complete_orthonormal_and_has_identity_projector() {
    let matrix = qslib_exact::DenseMatrix::new(2, vec![Complex64::new(0.0, 0.0); 4]).unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    assert_eq!(spectrum.values(), &[0.0, 0.0]);
    let projector = spectrum.projector(&[0, 1]).unwrap();
    assert_eq!(
        projector.as_slice(),
        &[
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ]
    );
    for left in 0..2 {
        for right in 0..2 {
            let overlap: Complex64 = spectrum
                .vector(left)
                .unwrap()
                .iter()
                .zip(spectrum.vector(right).unwrap())
                .map(|(a, b)| a.conj() * *b)
                .sum();
            assert!((overlap - Complex64::new((left == right) as u8 as f64, 0.0)).norm() < 1.0e-12);
        }
    }
    assert_eq!(
        evolve(
            &matrix,
            &[Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
            -1.0,
            false
        )
        .unwrap()[0],
        Complex64::new(1.0, 0.0)
    );
}

#[test]
fn nearby_distinct_eigenvalues_are_not_collapsed_as_a_degeneracy() {
    let matrix = qslib_exact::DenseMatrix::new(
        2,
        vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(5.0e-14, 0.0),
        ],
    )
    .unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    assert_eq!(spectrum.values(), &[0.0, 5.0e-14]);
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
fn stable_thermal_sums_and_signed_real_time_are_supported() {
    let matrix = qslib_exact::DenseMatrix::new(
        2,
        vec![
            Complex64::new(-1000.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(-999.0, 0.0),
        ],
    )
    .unwrap();
    let spectrum = diagonalize_hermitian(&matrix).unwrap();
    let thermal = ThermalSummary::from_spectrum(&spectrum, 1.0).unwrap();
    assert!((thermal.energy() + 1000.0 - 1.0 / (1.0 + 1.0_f64.exp())).abs() < 1.0e-12);
    assert!(thermal.log_partition_function() > 999.0);
    assert!(thermal.partition_function_overflowed());
    let state = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
    let backward = evolve(&matrix, &state, -0.25, false).unwrap();
    let forward = evolve(&matrix, &state, 0.25, false).unwrap();
    assert!((backward[0] - forward[0].conj()).norm() < 1.0e-12);
    let imaginary = evolve(
        &matrix,
        &[Complex64::new(1.0, 0.0), Complex64::new(1.0, 0.0)],
        1000.0,
        true,
    )
    .unwrap();
    assert!(
        imaginary
            .iter()
            .all(|value| value.re.is_finite() && value.im.is_finite())
    );
}

#[test]
fn hadamard_related_one_site_spectra_and_real_time_invariants_match() {
    let z = Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap(),
        )],
    )
    .unwrap();
    let x = Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::X)]).unwrap(),
        )],
    )
    .unwrap();
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let z_spectrum =
        diagonalize_hermitian(&qslib_exact::DenseMatrix::from_hamiltonian(&z, &basis).unwrap())
            .unwrap();
    let x_matrix = qslib_exact::DenseMatrix::from_hamiltonian(&x, &basis).unwrap();
    let x_spectrum = diagonalize_hermitian(&x_matrix).unwrap();
    assert_real_slices_close(z_spectrum.values(), x_spectrum.values(), 1.0e-12, 1.0e-12);
    let evolved = evolve(
        &x_matrix,
        &[Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        0.37,
        false,
    )
    .unwrap();
    assert!((evolved.iter().map(|value| value.norm_sqr()).sum::<f64>() - 1.0).abs() < 1.0e-12);
    let energy: Complex64 = evolved
        .iter()
        .enumerate()
        .map(|(row, amplitude)| amplitude.conj() * x_matrix.apply(&evolved).unwrap()[row])
        .sum();
    assert!(energy.norm() < 1.0e-12);
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

#[test]
fn sparse_solver_uses_restarts_and_fixed_sectors_reject_nonconserving_actions() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::X)]).unwrap(),
        )],
    )
    .unwrap();
    let full = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let sparse = qslib_exact::CsrMatrix::from_hamiltonian(&h, &full).unwrap();
    let ground = ground_state_sparse(&sparse, 1.0e-10, 4).unwrap();
    assert!((ground.energy() + 1.0).abs() < 1.0e-10);
    let sector = ExactBasis::fixed_weight(SiteCount::new(1).unwrap(), 0).unwrap();
    assert!(matches!(
        qslib_exact::DenseMatrix::from_hamiltonian(&h, &sector),
        Err(qslib_exact::ExactError::StateOutsideBasis { .. })
    ));
    assert!(qslib_exact::DenseMatrix::new(1, vec![Complex64::new(f64::NAN, 0.0)]).is_err());
}

#[test]
fn conserving_fixed_sector_spectrum_matches_the_full_matrix_subset() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(2).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::Z), (SiteId::new(1), Pauli::Z)]).unwrap(),
        )],
    )
    .unwrap();
    let full = ExactBasis::full(SiteCount::new(2).unwrap()).unwrap();
    let sector = ExactBasis::fixed_weight(SiteCount::new(2).unwrap(), 1).unwrap();
    let full_matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &full).unwrap();
    let sector_matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &sector).unwrap();
    let full_values = diagonalize_hermitian(&full_matrix)
        .unwrap()
        .values()
        .to_vec();
    let sector_values = diagonalize_hermitian(&sector_matrix)
        .unwrap()
        .values()
        .to_vec();
    assert_eq!(sector_values, vec![-1.0, -1.0]);
    assert_eq!(
        full_values.iter().filter(|value| **value == -1.0).count(),
        2
    );
}

#[test]
fn tfim_neutral_fixture_matrix_entries_match_the_independent_reference() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(2).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![
            (
                Complex64::new(-2.0, 0.0),
                PauliString::new(vec![(SiteId::new(0), Pauli::Z), (SiteId::new(1), Pauli::Z)])
                    .unwrap(),
            ),
            (
                Complex64::new(-0.5, 0.0),
                PauliString::new(vec![(SiteId::new(0), Pauli::X)]).unwrap(),
            ),
            (
                Complex64::new(-0.5, 0.0),
                PauliString::new(vec![(SiteId::new(1), Pauli::X)]).unwrap(),
            ),
        ],
    )
    .unwrap();
    let basis = ExactBasis::full(SiteCount::new(2).unwrap()).unwrap();
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let expected = [
        -2.0, -0.5, -0.5, 0.0, -0.5, 2.0, 0.0, -0.5, -0.5, 0.0, 2.0, -0.5, 0.0, -0.5, -0.5, -2.0,
    ];
    for (actual, expected) in matrix.as_slice().iter().zip(expected) {
        assert!((actual.re - expected).abs() < 1.0e-12);
        assert_eq!(actual.im, 0.0);
    }
}

#[test]
fn exact_expectation_and_variance_match_direct_matrix_evaluation() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 0.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap(),
        )],
    )
    .unwrap();
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let state = vec![
        Complex64::new(2.0_f64.sqrt().recip(), 0.0),
        Complex64::new(2.0_f64.sqrt().recip(), 0.0),
    ];
    assert!(qslib_exact::expectation(&matrix, &state).unwrap().norm() < 1.0e-12);
    assert!((qslib_exact::variance(&matrix, &state).unwrap().re - 1.0).abs() < 1.0e-12);
}

#[test]
fn variance_rejects_nonhermitian_matrix_instead_of_returning_cancellation_noise() {
    let matrix = qslib_exact::DenseMatrix::new(1, vec![Complex64::new(0.0, 1.0)]).unwrap();
    assert!(qslib_exact::variance(&matrix, &[Complex64::new(1.0, 0.0)]).is_err());
}

#[test]
fn convention_explicit_observable_definitions_keep_totals_densities_and_axes_distinct() {
    let energy = qslib_exact::energy_total_density(-6.0, 3).unwrap();
    assert_eq!((energy.total(), energy.density()), (-6.0, -2.0));
    let (axis, magnetization) =
        qslib_exact::magnetization_total_density(qslib_core::PhysicalAxis::Z, &[1.0, -1.0, 1.0])
            .unwrap();
    assert_eq!(axis, qslib_core::PhysicalAxis::Z);
    assert_eq!(magnetization.total(), 1.0);
    let (raw, connected) = qslib_exact::raw_connected_correlation(0.5, 0.5, -0.5).unwrap();
    assert_eq!((raw, connected), (0.5, 0.75));
    assert_eq!(
        qslib_exact::quantum_fisher_information(0.25, 2)
            .unwrap()
            .total(),
        1.0
    );
    assert_eq!(
        qslib_exact::thermal_heat_capacity(2.0, 0.5, 2)
            .unwrap()
            .density(),
        1.0
    );
    let sublattice =
        qslib_exact::sublattice_moment(&[1.0, -1.0, 1.0, -1.0], &[1.0, -1.0, 1.0, -1.0]).unwrap();
    assert_eq!(
        (
            sublattice.signed(),
            sublattice.absolute(),
            sublattice.squared()
        ),
        (1.0, 1.0, 1.0)
    );
    let structure = qslib_exact::structure_factor(
        qslib_core::PhysicalAxis::Z,
        &[(0.0, 0.0), (1.0, 0.0)],
        (0.0, 0.0),
        &[1.0, 0.5, 0.5, 1.0],
        None,
    )
    .unwrap();
    assert_eq!(structure.axis(), qslib_core::PhysicalAxis::Z);
    assert!(!structure.connected());
    assert_eq!(structure.value(), Complex64::new(1.5, 0.0));
    assert!(qslib_exact::quantum_fisher_information(-0.25, 2).is_err());
}

#[test]
fn pauli_magnetization_and_total_spin_match_direct_matrix_values() {
    let basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let state = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
    let z = PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap();
    assert_eq!(
        qslib_exact::pauli_expectation(&z, &basis, &state).unwrap(),
        Complex64::new(1.0, 0.0)
    );
    let magnetization =
        qslib_exact::magnetization_expectation(qslib_core::PhysicalAxis::Z, &basis, &state)
            .unwrap();
    assert_eq!((magnetization.total(), magnetization.density()), (1.0, 1.0));
    assert_eq!(
        qslib_exact::magnetization_component_magnitude([1.0, 0.0, 0.0]).unwrap(),
        1.0
    );
    assert!(qslib_exact::magnetization_component_magnitude([f64::MAX, f64::MAX, 0.0]).is_err());
}

#[test]
fn product_singlet_and_enumerated_entropy_observables_match_analytic_values() {
    let product = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.0, 0.0),
    ];
    assert!(qslib_exact::shannon_entropy(&product).unwrap().abs() < 1.0e-12);
    assert!(
        qslib_exact::bipartite_entropy(&product, &[0])
            .unwrap()
            .abs()
            < 1.0e-12
    );
    let singlet = [
        Complex64::new(0.0, 0.0),
        Complex64::new(2.0_f64.sqrt().recip(), 0.0),
        Complex64::new(-2.0_f64.sqrt().recip(), 0.0),
        Complex64::new(0.0, 0.0),
    ];
    assert!(
        (qslib_exact::bipartite_entropy(&singlet, &[0]).unwrap() - 2.0_f64.ln()).abs() < 1.0e-12
    );
    let basis = ExactBasis::full(SiteCount::new(2).unwrap()).unwrap();
    let zz =
        qslib_exact::correlation_expectation(qslib_core::PhysicalAxis::Z, 0, 1, &basis, &singlet)
            .unwrap();
    assert!((zz.re + 1.0).abs() < 1.0e-12);
    assert_eq!(
        qslib_exact::correlation_expectation(qslib_core::PhysicalAxis::Z, 0, 0, &basis, &singlet,)
            .unwrap(),
        Complex64::new(1.0, 0.0)
    );
    assert!(
        qslib_exact::total_spin_squared(&basis, &singlet)
            .unwrap()
            .abs()
            < 1.0e-12
    );
    let one_basis = ExactBasis::full(SiteCount::new(1).unwrap()).unwrap();
    let one_state = [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
    assert_eq!(
        qslib_exact::spin_magnetization_expectation(
            qslib_core::PhysicalAxis::Z,
            &one_basis,
            &one_state,
        )
        .unwrap()
        .total(),
        0.5
    );
    assert!(
        (qslib_exact::total_spin_squared(&one_basis, &one_state).unwrap() - 0.75).abs() < 1.0e-12
    );
    assert!(qslib_exact::total_spin_squared(&one_basis, &[Complex64::new(1.0, 0.0)]).is_err());
}

#[test]
fn bipartite_entropy_rejects_empty_non_power_of_two_and_duplicate_subsystems() {
    assert!(qslib_exact::bipartite_entropy(&[], &[0]).is_err());
    assert!(qslib_exact::bipartite_entropy(&[Complex64::new(1.0, 0.0); 3], &[0]).is_err());
    let state = vec![Complex64::new(1.0, 0.0); 4];
    assert!(qslib_exact::bipartite_entropy(&state, &[0, 0]).is_err());
}
