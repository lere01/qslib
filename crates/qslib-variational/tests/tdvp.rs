use qslib_core::{
    BasisBit, BasisState, Complex64, Hamiltonian, Pauli, PauliString, SiteCount, SiteId,
};
use qslib_variational::{
    DenseQgt, ParameterLayout, Regularization, TDVPMode, TDVPSolveOptions, estimate_tdvp,
    gcv_tikhonov_scores, projected_residual_squared, solve_cg,
};

fn batch() -> (Vec<f64>, Vec<Complex64>, Vec<Complex64>) {
    let weights = vec![1.0, 2.0, 1.0, 2.0];
    let energies = vec![
        Complex64::new(0.0, 1.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(4.0, -1.0),
        Complex64::new(6.0, 2.0),
    ];
    let derivatives = vec![
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(1.0, 0.0),
        Complex64::new(1.0, 1.0),
        Complex64::new(1.0, -1.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(0.0, 2.0),
    ];
    (weights, energies, derivatives)
}

#[test]
fn two_sample_reference_pins_qgt_and_force_entries() {
    let weights = [1.0, 1.0];
    let energies = [Complex64::new(0.0, 1.0), Complex64::new(2.0, 3.0)];
    let derivatives = [
        Complex64::new(1.0, 1.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(3.0, -1.0),
        Complex64::new(2.0, 0.0),
    ];
    let statistics =
        estimate_tdvp(&weights, &energies, &derivatives, 2, TDVPMode::RealTime).unwrap();
    assert_eq!(statistics.qgt().as_slice(), &[2.0, 1.0, 1.0, 1.0]);
    assert_eq!(statistics.force().real, vec![0.0, 1.0]);
    assert_eq!(statistics.force().imag, vec![2.0, 1.0]);
    assert_eq!(statistics.rhs(), &[2.0, 1.0]);
}

#[test]
fn tdvp_statistics_match_weighted_centered_definitions_and_mode_signs() {
    let (weights, energies, derivatives) = batch();
    let statistics =
        estimate_tdvp(&weights, &energies, &derivatives, 2, TDVPMode::RealTime).unwrap();
    assert_eq!(statistics.samples(), 4);
    assert_eq!(statistics.parameters(), 2);
    assert!((statistics.energy_mean().re - 10.0 / 3.0).abs() < 1.0e-12);
    assert!(statistics.energy_variance() > 0.0);
    let force = statistics.force();
    let imaginary_rhs = estimate_tdvp(
        &weights,
        &energies,
        &derivatives,
        2,
        TDVPMode::ImaginaryTime,
    )
    .unwrap();
    for index in 0..2 {
        assert!((imaginary_rhs.rhs()[index] + force.real[index]).abs() < 1.0e-12);
        assert!((statistics.rhs()[index] - force.imag[index]).abs() < 1.0e-12);
    }
    assert!(statistics.qgt().is_positive_semidefinite(1.0e-12));
}

#[test]
fn dense_qgt_matvec_and_cg_solve_match_manufactured_system() {
    let qgt = DenseQgt::new(2, vec![4.0, 1.0, 1.0, 3.0]).unwrap();
    assert_eq!(qgt.matvec(&[2.0, -1.0]).unwrap(), vec![7.0, -1.0]);
    let rhs = vec![1.0, 2.0];
    let dense = qgt.solve(&rhs, 0.0, 1.0e-12, 20).unwrap();
    let matrix_free =
        solve_cg(|vector| Ok(qgt.matvec(vector).unwrap()), &rhs, 1.0e-12, 20).unwrap();
    for (left, right) in dense.direction().iter().zip(matrix_free.direction()) {
        assert!((left - right).abs() < 1.0e-10);
    }
    assert!(dense.converged());
}

#[test]
fn regularization_gcv_and_clipping_are_explicit_and_diagnostic() {
    let (weights, energies, derivatives) = batch();
    let statistics =
        estimate_tdvp(&weights, &energies, &derivatives, 2, TDVPMode::RealTime).unwrap();
    let options = TDVPSolveOptions {
        regularization: Regularization::Gcv {
            lambda_grid: vec![1.0e-4, 0.1, 1.0],
        },
        max_update_norm: Some(0.01),
        ..TDVPSolveOptions::default()
    };
    let result = statistics.solve(options).unwrap();
    assert!(result.clipped());
    assert!(result.residual_squared().is_finite());
    assert!(result.normalized_residual().is_finite());
    assert!(result.selected_lambda().is_some());
    assert!(result.linear_residual_squared().is_finite());
    assert!(result.zero_variance_policy().contains("Var(H)=0"));
    let scores = gcv_tikhonov_scores(&[1.0, 2.0], &[1.0, 0.5], &[0.1, 1.0]).unwrap();
    assert_eq!(scores.len(), 2);
}

#[test]
fn spectral_cutoff_keeps_endpoint_and_degenerate_maximum_modes() {
    let statistics = estimate_tdvp(
        &[1.0, 1.0],
        &[Complex64::new(0.0, 1.0), Complex64::new(2.0, 3.0)],
        &[
            Complex64::new(1.0, 1.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(3.0, -1.0),
            Complex64::new(2.0, 0.0),
        ],
        2,
        TDVPMode::RealTime,
    )
    .unwrap();
    let endpoint = statistics
        .solve(TDVPSolveOptions {
            regularization: Regularization::SpectralCutoff {
                relative_cutoff: 1.0,
            },
            ..TDVPSolveOptions::default()
        })
        .unwrap();
    assert!(endpoint.direction().iter().any(|value| value.abs() > 0.0));
    let (floor, applied) = projected_residual_squared(1.0, 1.0, 1.0 + 1.0e-13).unwrap();
    assert_eq!(floor, 0.0);
    assert!(applied);
    assert_eq!(
        projected_residual_squared(1.0, 1.0, 0.0).unwrap(),
        (2.0, false)
    );
    assert!(projected_residual_squared(1.0, 0.0, 2.0).is_err());
    for scale in [1.0e-20_f64, 1.0e20_f64] {
        let amplitude = scale.sqrt();
        let scaled = estimate_tdvp(
            &[1.0, 1.0, 1.0, 1.0],
            &[
                Complex64::new(0.0, 1.0),
                Complex64::new(0.0, -1.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(-1.0, 0.0),
            ],
            &[
                Complex64::new(amplitude, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(-amplitude, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, amplitude),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, -amplitude),
            ],
            2,
            TDVPMode::RealTime,
        )
        .unwrap();
        let result = scaled
            .solve(TDVPSolveOptions {
                regularization: Regularization::SpectralCutoff {
                    relative_cutoff: 1.0,
                },
                ..TDVPSolveOptions::default()
            })
            .unwrap();
        assert!(result.direction().iter().all(|value| value.is_finite()));
        assert!((result.direction()[0] * amplitude - 1.0).abs() < 1.0e-8);
        assert!((result.direction()[1] * amplitude + 1.0).abs() < 1.0e-8);
    }
    let degenerate = DenseQgt::new(2, vec![3.0, 0.0, 0.0, 3.0]).unwrap();
    assert_eq!(degenerate.spectrum().unwrap().eigenvalues, vec![3.0, 3.0]);
    let degenerate_statistics = estimate_tdvp(
        &[1.0, 1.0, 1.0, 1.0],
        &[
            Complex64::new(0.0, 1.0),
            Complex64::new(0.0, -1.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(-1.0, 0.0),
        ],
        &[
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(-1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, -1.0),
        ],
        2,
        TDVPMode::RealTime,
    )
    .unwrap();
    let degenerate_result = degenerate_statistics
        .solve(TDVPSolveOptions {
            regularization: Regularization::SpectralCutoff {
                relative_cutoff: 1.0,
            },
            ..TDVPSolveOptions::default()
        })
        .unwrap();
    assert!(
        degenerate_result
            .direction()
            .iter()
            .all(|value| value.is_finite())
    );
    assert!((degenerate_result.direction()[0] - 1.0).abs() < 1.0e-12);
    assert!((degenerate_result.direction()[1] + 1.0).abs() < 1.0e-12);
}

#[test]
fn streamed_qgt_and_ratio_local_energy_match_the_reference_contract() {
    let rows_a = [Complex64::new(1.0, 1.0), Complex64::new(0.0, 0.0)];
    let rows_b = [Complex64::new(3.0, -1.0), Complex64::new(2.0, 0.0)];
    let product = qslib_variational::qgt_vector_product_stream(
        [
            qslib_variational::QgtSampleChunk {
                weights: &[1.0],
                derivatives: &rows_a,
            },
            qslib_variational::QgtSampleChunk {
                weights: &[1.0],
                derivatives: &rows_b,
            },
        ],
        &[Complex64::new(2.0, 0.0), Complex64::new(1.0, 0.0)],
        2,
        &[0.5, -1.0],
    )
    .unwrap();
    assert_eq!(product, vec![0.0, -0.5]);
    let dense = DenseQgt::new(2, vec![2.0, 1.0, 1.0, 1.0]).unwrap();
    assert_eq!(dense.matvec(&[0.5, -1.0]).unwrap(), product);
    let streamed_solve = solve_cg(
        |vector| {
            qslib_variational::qgt_vector_product_stream(
                [
                    qslib_variational::QgtSampleChunk {
                        weights: &[1.0],
                        derivatives: &rows_a,
                    },
                    qslib_variational::QgtSampleChunk {
                        weights: &[1.0],
                        derivatives: &rows_b,
                    },
                ],
                &[Complex64::new(2.0, 0.0), Complex64::new(1.0, 0.0)],
                2,
                vector,
            )
        },
        &[2.0, 1.0],
        1.0e-12,
        10,
    )
    .unwrap();
    let dense_solve = dense.solve(&[2.0, 1.0], 0.0, 1.0e-12, 10).unwrap();
    assert_eq!(streamed_solve.direction(), dense_solve.direction());
    let energy = qslib_variational::local_energy_from_ratios(
        Complex64::new(1.0, 0.0),
        &[
            (Complex64::new(2.0, 0.0), Complex64::new(3.0, 0.0)),
            (Complex64::new(0.0, 1.0), Complex64::new(0.0, 2.0)),
        ],
    )
    .unwrap();
    assert_eq!(energy, Complex64::new(5.0, 0.0));
    let hamiltonian = Hamiltonian::new(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(
            Complex64::new(1.0, 1.0),
            PauliString::new(vec![(SiteId::new(0), Pauli::Y)]).unwrap(),
        )],
    )
    .unwrap();
    let zero = BasisState::from_bits(&[BasisBit::Zero]).unwrap();
    let one = BasisState::from_bits(&[BasisBit::One]).unwrap();
    let direct = hamiltonian
        .local_energy(
            &zero,
            &[
                (zero.clone(), Complex64::new(1.0, 0.0)),
                (one, Complex64::new(3.0, 0.0)),
            ],
        )
        .unwrap();
    assert_eq!(
        direct,
        qslib_variational::local_energy_from_ratios(
            Complex64::new(0.0, 0.0),
            &[(Complex64::new(1.0, -1.0), Complex64::new(3.0, 0.0))],
        )
        .unwrap()
    );
}

#[test]
fn projected_residual_has_analytic_zero_and_invalid_cases_are_rejected() {
    let energies = [Complex64::new(0.0, 1.0), Complex64::new(2.0, 3.0)];
    let derivatives = [
        Complex64::new(1.0, 1.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(3.0, -1.0),
        Complex64::new(2.0, 0.0),
    ];
    let statistics =
        estimate_tdvp(&[1.0, 1.0], &energies, &derivatives, 2, TDVPMode::RealTime).unwrap();
    let result = statistics
        .solve(TDVPSolveOptions {
            regularization: Regularization::SpectralCutoff {
                relative_cutoff: 0.0,
            },
            ..TDVPSolveOptions::default()
        })
        .unwrap();
    assert!(result.residual_squared() < 1.0e-12);
    assert!(!result.clipped());
    assert!(DenseQgt::new(0, Vec::new()).is_err());
    assert!(DenseQgt::new(2, vec![1.0, 0.0, 1.0, 1.0]).is_err());
    assert!(
        DenseQgt::new(2, vec![1.0, 0.0, 0.0, 1.0])
            .unwrap()
            .matvec(&[f64::NAN, 0.0])
            .is_err()
    );
    let matrix = DenseQgt::new(2, vec![4.0, 1.0, 1.0, 3.0]).unwrap();
    let failed = solve_cg(|vector| matrix.matvec(vector), &[1.0, 2.0], 1.0e-12, 1).unwrap();
    assert!(!failed.converged());
    let breakdown = solve_cg(
        |vector| Ok(vec![0.0; vector.len()]),
        &[1.0, 0.0],
        1.0e-12,
        4,
    )
    .unwrap();
    assert!(!breakdown.converged());
}

#[test]
fn parameter_layout_fingerprint_is_order_and_shape_deterministic() {
    let first = ParameterLayout::new(vec![("a", vec![2]), ("b", vec![1, 1])]).unwrap();
    let same = ParameterLayout::new(vec![("a", vec![2]), ("b", vec![1, 1])]).unwrap();
    let different = ParameterLayout::new(vec![("b", vec![1, 1]), ("a", vec![2])]).unwrap();
    assert_eq!(first.fingerprint(), same.fingerprint());
    assert_eq!(
        first.fingerprint(),
        "blake3-v1-2af2bbc3bcfc15356f1150bf0cb41e2d6bfab9a1e4a83188dfdc23cf5696ad85"
    );
    assert_ne!(first.fingerprint(), different.fingerprint());
    assert_eq!(first.parameters(), 3);
}
