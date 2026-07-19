use qslib_core::BasisBit;
use qslib_sse::{
    BasisSseState, LocalSseModel, Operator, OperatorKind, SimulationConfig, SseModel,
    SseModelError, SseSampler, SseTerm, ThermodynamicAccumulator,
};
use rand_chacha::ChaCha8Rng;
use rand_core::SeedableRng;

#[test]
fn tfim_decomposition_reconstructs_canonical_z_energy() {
    let model = LocalSseModel::tfim(2, &[(0, 1)], 1.25, 0.75).unwrap();
    for bits in [
        [BasisBit::Zero, BasisBit::Zero],
        [BasisBit::Zero, BasisBit::One],
        [BasisBit::One, BasisBit::Zero],
        [BasisBit::One, BasisBit::One],
    ] {
        let diagonal_sum: f64 = model
            .diagonal_term_indices()
            .iter()
            .map(|&index| model.matrix_element(index as usize, &bits).unwrap())
            .sum();
        let z0 = bits[0].pauli_eigenvalue() as f64;
        let z1 = bits[1].pauli_eigenvalue() as f64;
        assert!((model.energy_shift() - diagonal_sum + 1.25 * z0 * z1).abs() < 1.0e-12);
    }
}

#[test]
fn tfim_weighted_constructor_preserves_pair_and_site_coefficients() {
    let model =
        LocalSseModel::tfim_weighted(3, &[(0, 1, 0.5), (1, 2, 1.75)], &[0.1, 0.2, 0.3]).unwrap();
    let bonds = model
        .terms()
        .iter()
        .filter_map(|term| match term {
            SseTerm::TfimBond { coupling, .. } => Some(*coupling),
            _ => None,
        })
        .collect::<Vec<_>>();
    let fields = model
        .terms()
        .iter()
        .filter_map(|term| match term {
            SseTerm::SpinFlip { amplitude, .. } => Some(*amplitude),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(bonds, vec![0.5, 1.75]);
    assert_eq!(fields, vec![0.1, 0.2, 0.3]);
}

#[test]
fn rydberg_decomposition_uses_bit_one_as_occupation() {
    let model = LocalSseModel::rydberg(2, &[0.0, 1.0], &[(0, 1, 2.0)], 1.2).unwrap();
    let bits = [BasisBit::One, BasisBit::Zero];
    let diagonal_sum: f64 = model
        .diagonal_term_indices()
        .iter()
        .map(|&index| model.matrix_element(index as usize, &bits).unwrap())
        .sum();
    assert!((model.energy_shift() - diagonal_sum).abs() < 1.0e-12);
    assert_eq!(model.occupation(&bits, 0).unwrap(), 1.0);
}

#[test]
fn rydberg_pair_coefficients_remain_distinct_and_signed() {
    let model =
        LocalSseModel::rydberg(3, &[0.0, 0.0, 0.0], &[(0, 1, 2.0), (0, 2, -3.0)], 0.0).unwrap();
    let bits = [BasisBit::One, BasisBit::One, BasisBit::One];
    let pair_values = model
        .terms()
        .iter()
        .enumerate()
        .filter(|(_, term)| matches!(term, SseTerm::RydbergInteraction { .. }))
        .map(|(index, _)| model.matrix_element(index, &bits).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(pair_values, vec![0.0, 3.0]);
}

#[test]
fn propagation_requires_trace_closure_and_preserves_operator_kinds() {
    let model = LocalSseModel::tfim(1, &[], 0.0, 0.75).unwrap();
    let state = BasisSseState::new(
        vec![BasisBit::Zero],
        vec![Operator::off_diagonal(1), Operator::off_diagonal(1)],
    )
    .unwrap();
    let result = state.propagate(&model).unwrap();
    assert!(result.trace_closed);
    assert_eq!(result.final_state, vec![BasisBit::Zero]);
    assert_eq!(result.operator_count, 2);
    assert_eq!(model.operator_kind(1).unwrap(), OperatorKind::OffDiagonal);
    assert!(state.validate_trace(&model).is_ok());
}

#[test]
fn thermodynamic_accumulator_matches_expansion_order_moments() {
    let mut accumulator = ThermodynamicAccumulator::default();
    accumulator.record(2);
    accumulator.record(4);
    let result = accumulator.results(2.0, 4.0, 2).unwrap();
    assert_eq!(result.samples, 2);
    assert!((result.energy - 2.5).abs() < 1.0e-12);
    assert!(result.heat_capacity.is_finite());
}

#[test]
fn unsupported_signs_and_invalid_terms_fail_explicitly() {
    assert!(matches!(
        LocalSseModel::tfim(2, &[(0, 1)], -1.0, 0.5),
        Err(SseModelError::UnsupportedSign { .. })
    ));
    assert!(matches!(
        LocalSseModel::tfim(2, &[(0, 2)], 1.0, 0.5),
        Err(SseModelError::InvalidSite { .. })
    ));
    let mut short = vec![BasisBit::Zero];
    assert!(matches!(
        model_for_test().apply_off_diagonal(1, &mut short),
        Err(SseModelError::InvalidLength { .. })
    ));
}

fn model_for_test() -> LocalSseModel {
    LocalSseModel::tfim(2, &[(0, 1)], 1.0, 0.5).unwrap()
}

#[test]
fn seeded_diagonal_sampler_preserves_trace_and_reports_statistics() {
    let model = LocalSseModel::tfim(2, &[(0, 1)], 1.0, 0.5).unwrap();
    let state = BasisSseState::new(
        vec![BasisBit::Zero, BasisBit::Zero],
        vec![Operator::identity(); 32],
    )
    .unwrap();
    let mut sampler = SseSampler::new(model, state, 0.8, ChaCha8Rng::seed_from_u64(12)).unwrap();
    let results = sampler
        .run(SimulationConfig {
            thermalization_sweeps: 2,
            measurement_sweeps: 4,
            sweeps_per_measurement: 1,
        })
        .unwrap();
    assert_eq!(results.thermodynamics.samples, 4);
    assert!(results.thermodynamics.energy.is_finite());
    assert!(sampler.state().validate_trace(sampler.model()).is_ok());
}

#[test]
fn one_site_tfim_thermal_energy_matches_exact_limit() {
    let field = 0.7;
    let beta = 1.1;
    let model = LocalSseModel::tfim(1, &[], 0.0, field).unwrap();
    let state = BasisSseState::new(vec![BasisBit::Zero], vec![Operator::identity(); 96]).unwrap();
    let mut sampler = SseSampler::new(model, state, beta, ChaCha8Rng::seed_from_u64(99)).unwrap();
    let result = sampler
        .run(SimulationConfig {
            thermalization_sweeps: 2_000,
            measurement_sweeps: 8_000,
            sweeps_per_measurement: 1,
        })
        .unwrap();
    let exact = -field * (beta * field).tanh();
    assert!((result.thermodynamics.energy - exact).abs() < 0.12);
    assert!(result.off_diagonal.accepted > 0);
}

#[test]
fn diagonal_updates_use_the_propagated_time_slice_state() {
    let model = LocalSseModel::rydberg(1, &[1.0], &[], 0.5).unwrap();
    let state = BasisSseState::new(
        vec![BasisBit::Zero],
        vec![
            Operator::off_diagonal(1),
            Operator::diagonal(2),
            Operator::off_diagonal(1),
        ],
    )
    .unwrap();
    let mut sampler = SseSampler::new(model, state, 0.8, ChaCha8Rng::seed_from_u64(7)).unwrap();
    assert!(sampler.diagonal_sweep().is_ok());
    assert!(sampler.state().validate_trace(sampler.model()).is_ok());
}
