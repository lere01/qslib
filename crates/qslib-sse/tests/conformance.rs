use qslib_core::{
    BasisBit, Bond, InteractionChannel, QSLIB_SEED_SCHEME, SiteId, WeightedInteraction,
    derive_seed, expand_master_seed,
};
use qslib_sse::{
    BasisSseState, LegacyModelKind, LegacySpin, LocalSseModel, Operator, OperatorKind,
    SimulationConfig, SimulationResults, SseModel, SseModelError, SseSampler, SseTerm,
    ThermodynamicAccumulator, convert_legacy_bits, derive_chain_seed, derive_legacy_chain_seed,
    logical_chain_seeds, run_parallel_chains,
};
use rand_chacha::ChaCha20Rng;
use rand_core::SeedableRng;
use serde::Deserialize;

fn independent_energy_estimate(results: &[SimulationResults]) -> (f64, f64) {
    assert!(results.len() >= 2);
    let count = results.len() as f64;
    let mean = results
        .iter()
        .map(|result| result.thermodynamics.energy)
        .sum::<f64>()
        / count;
    let variance = results
        .iter()
        .map(|result| (result.thermodynamics.energy - mean).powi(2))
        .sum::<f64>()
        / (count - 1.0);
    (mean, (variance / count).sqrt())
}

#[derive(Deserialize)]
struct ParityFixture {
    schema: String,
    provenance: Provenance,
    legacy_chain_seed: u64,
    sweep_comparison_policy: String,
    standalone_capture: StandaloneCapture,
    tfim: TfimParity,
    rydberg: RydbergParity,
}
#[derive(Deserialize)]
struct Provenance {
    source: String,
    revision: String,
    capture: String,
    seed_scheme: String,
    command: String,
    config: String,
}
#[derive(Deserialize)]
struct StandaloneCapture {
    chains: usize,
    threads: usize,
    master_seed: u64,
    thermalization_sweeps: usize,
    measurement_sweeps: usize,
    energy_per_site: f64,
    chain_standard_error: f64,
    chain_seeds: Vec<u64>,
    energy_per_chain: Vec<f64>,
    insertions_accepted: Vec<usize>,
    removals_accepted: Vec<usize>,
    off_diagonal_proposals_accepted: Vec<usize>,
}
#[derive(Deserialize)]
struct TfimParity {
    num_sites: usize,
    bond: [u32; 2],
    coupling: f64,
    field: f64,
    energy_shift: f64,
    term_count: usize,
    propagation: PropagationParity,
}
#[derive(Deserialize)]
struct PropagationParity {
    initial_bits: Vec<u8>,
    operator: (String, u32),
    final_bits: Vec<u8>,
    log_weight: f64,
}
#[derive(Deserialize)]
struct RydbergParity {
    num_sites: usize,
    detunings: Vec<f64>,
    interaction: [f64; 3],
    omega: f64,
    energy_shift: f64,
    term_count: usize,
}

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
fn standalone_parity_fixture_matches_resolved_models_and_propagation() {
    let fixture: ParityFixture =
        serde_json::from_str(include_str!("fixtures/sse_parity_v1.json")).unwrap();
    assert_eq!(fixture.schema, "qslib-sse-parity-v1");
    assert!(fixture.provenance.source.contains("standalone-sse"));
    assert_eq!(
        fixture.provenance.revision,
        "f6b7650e7f2af230a77ee48a0613d388448bea29"
    );
    assert_eq!(fixture.provenance.capture, "2026-07-19");
    assert_eq!(fixture.provenance.seed_scheme, "legacy-splitmix-u64");
    assert!(
        fixture
            .provenance
            .command
            .contains("configs/tfim-chain.yaml")
    );
    assert_eq!(fixture.provenance.config, "sse/configs/tfim-chain.yaml");
    assert_eq!(fixture.legacy_chain_seed, derive_legacy_chain_seed(42, 3));
    assert!(
        fixture
            .sweep_comparison_policy
            .contains("different proposal kernels")
    );
    assert_eq!(fixture.standalone_capture.chains, 4);
    assert_eq!(fixture.standalone_capture.threads, 4);
    assert_eq!(fixture.standalone_capture.master_seed, 24301);
    assert_eq!(fixture.standalone_capture.thermalization_sweeps, 1000);
    assert_eq!(fixture.standalone_capture.measurement_sweeps, 10000);
    assert!((fixture.standalone_capture.energy_per_site + 1.0630828125).abs() < 1.0e-12);
    assert!(
        (fixture.standalone_capture.chain_standard_error - 0.0022881635749160936).abs() < 1.0e-15
    );
    assert_eq!(fixture.standalone_capture.chain_seeds.len(), 4);
    for (index, seed) in fixture.standalone_capture.chain_seeds.iter().enumerate() {
        assert_eq!(*seed, derive_legacy_chain_seed(24301, index as u64));
    }
    assert_eq!(fixture.standalone_capture.energy_per_chain.len(), 4);
    assert_eq!(fixture.standalone_capture.insertions_accepted.len(), 4);
    assert_eq!(fixture.standalone_capture.removals_accepted.len(), 4);
    assert_eq!(
        fixture.standalone_capture.off_diagonal_proposals_accepted,
        vec![10000; 4]
    );
    let tfim = LocalSseModel::tfim(
        fixture.tfim.num_sites,
        &[tuple_pair(fixture.tfim.bond)],
        fixture.tfim.coupling,
        fixture.tfim.field,
    )
    .unwrap();
    assert_eq!(tfim.terms().len(), fixture.tfim.term_count);
    assert!((tfim.energy_shift() - fixture.tfim.energy_shift).abs() < 1.0e-12);
    let tfim_state = BasisSseState::new(
        fixture
            .tfim
            .propagation
            .initial_bits
            .iter()
            .map(|bit| {
                if *bit == 0 {
                    BasisBit::Zero
                } else {
                    BasisBit::One
                }
            })
            .collect(),
        vec![Operator::diagonal(fixture.tfim.propagation.operator.1)],
    )
    .unwrap();
    let propagation = tfim_state.propagate(&tfim).unwrap();
    assert_eq!(fixture.tfim.propagation.operator.0, "diagonal");
    assert_eq!(propagation.final_state, tfim_state.bits());
    assert_eq!(
        propagation
            .final_state
            .iter()
            .map(|bit| if *bit == BasisBit::Zero { 0 } else { 1 })
            .collect::<Vec<_>>(),
        fixture.tfim.propagation.final_bits
    );
    assert!((propagation.log_weight - fixture.tfim.propagation.log_weight).abs() < 1.0e-12);

    let rydberg = LocalSseModel::rydberg(
        fixture.rydberg.num_sites,
        &fixture.rydberg.detunings,
        &[(
            fixture.rydberg.interaction[0] as u32,
            fixture.rydberg.interaction[1] as u32,
            fixture.rydberg.interaction[2],
        )],
        fixture.rydberg.omega,
    )
    .unwrap();
    assert_eq!(rydberg.terms().len(), fixture.rydberg.term_count);
    assert!((rydberg.energy_shift() - fixture.rydberg.energy_shift).abs() < 1.0e-12);
}

#[test]
fn converted_four_site_standalone_capture_matches_qslib_thermal_energy() {
    let fixture: ParityFixture =
        serde_json::from_str(include_str!("fixtures/sse_parity_v1.json")).unwrap();
    let capture = &fixture.standalone_capture;
    let model = LocalSseModel::tfim(4, &[(0, 1), (1, 2), (2, 3), (3, 0)], 1.0, 0.5).unwrap();
    let state = BasisSseState::new(
        vec![BasisBit::Zero, BasisBit::One, BasisBit::Zero, BasisBit::One],
        vec![Operator::identity(); 128],
    )
    .unwrap();
    let results = run_parallel_chains(
        model,
        state,
        4.0,
        SimulationConfig {
            thermalization_sweeps: capture.thermalization_sweeps,
            measurement_sweeps: capture.measurement_sweeps,
            sweeps_per_measurement: 1,
        },
        capture.master_seed,
        capture.chains,
        capture.threads,
    )
    .unwrap();
    let chain_means = results
        .iter()
        .map(|result| result.thermodynamics.energy_per_site)
        .collect::<Vec<_>>();
    assert!(results.iter().all(|result| {
        result.thermodynamics.samples == capture.measurement_sweeps as u64
            && result.diagonal.insertions_accepted > 0
            && result.diagonal.removals_accepted > 0
            && result.basis.accepted > 0
    }));
    let qslib_mean = chain_means.iter().sum::<f64>() / chain_means.len() as f64;
    let qslib_variance = chain_means
        .iter()
        .map(|value| (value - qslib_mean).powi(2))
        .sum::<f64>()
        / (chain_means.len() as f64 - 1.0);
    let qslib_se = (qslib_variance / chain_means.len() as f64).sqrt();
    let combined_se = (qslib_se.powi(2) + capture.chain_standard_error.powi(2)).sqrt();
    assert!(
        (qslib_mean - capture.energy_per_site).abs() < 3.0 * combined_se,
        "qslib {qslib_mean} +/- {qslib_se} versus standalone {} +/- {}",
        capture.energy_per_site,
        capture.chain_standard_error
    );
}

fn tuple_pair(pair: [u32; 2]) -> (u32, u32) {
    (pair[0], pair[1])
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
fn tfim_resolved_constructor_consumes_canonical_weighted_interactions() {
    let interaction = WeightedInteraction::new(
        Bond::new(SiteId::new(0), SiteId::new(1)).unwrap(),
        InteractionChannel::IsingZZ,
        1.75,
    )
    .unwrap();
    let model = LocalSseModel::tfim_resolved(2, &[interaction], &[0.2, 0.3]).unwrap();
    assert!(model.terms().iter().any(|term| matches!(
        term,
        SseTerm::TfimBond { coupling, .. } if (*coupling - 1.75).abs() < 1.0e-12
    )));
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
fn rydberg_resolved_constructor_requires_density_density_channels() {
    let interaction = WeightedInteraction::new(
        Bond::new(SiteId::new(0), SiteId::new(1)).unwrap(),
        InteractionChannel::RydbergDensityDensity,
        -2.0,
    )
    .unwrap();
    let model = LocalSseModel::rydberg_resolved(2, &[0.1, 0.2], &[interaction], 0.5).unwrap();
    assert!(model.terms().iter().any(|term| matches!(
        term,
        SseTerm::RydbergInteraction { interaction, .. } if (*interaction + 2.0).abs() < 1.0e-12
    )));
}

#[test]
fn canonical_chain_seed_uses_versioned_domain_separated_32_byte_streams() {
    let master = expand_master_seed(42);
    let chain = derive_seed(&master, "sse_chain", &[3]);
    assert_eq!(QSLIB_SEED_SCHEME, "qslib-seed-v1");
    assert_eq!(chain.len(), 32);
    assert_eq!(
        chain,
        [
            220, 238, 78, 156, 26, 9, 222, 111, 111, 9, 7, 120, 219, 151, 73, 193, 227, 58, 173,
            168, 25, 200, 84, 59, 226, 82, 251, 236, 153, 61, 77, 89,
        ]
    );
    assert_ne!(chain, derive_seed(&master, "sse_chain", &[4]));
    assert_ne!(chain, derive_seed(&master, "disorder", &[3]));
    assert_eq!(chain, derive_chain_seed(42, 3));
}

#[test]
fn legacy_chain_seed_is_explicitly_separate_from_canonical_seed() {
    assert_eq!(derive_legacy_chain_seed(42, 3), 6_904_877_152_625_194_467);
    assert_ne!(
        derive_legacy_chain_seed(42, 3).to_le_bytes().len(),
        derive_chain_seed(42, 3).len()
    );
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
    let mut sampler = SseSampler::new(model, state, 0.8, ChaCha20Rng::seed_from_u64(12)).unwrap();
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
    let results = run_parallel_chains(
        model,
        state,
        beta,
        SimulationConfig {
            thermalization_sweeps: 500,
            measurement_sweeps: 2_000,
            sweeps_per_measurement: 1,
        },
        99,
        8,
        3,
    )
    .unwrap();
    let (energy, standard_error) = independent_energy_estimate(&results);
    let exact = -field * (beta * field).tanh();
    assert!((energy - exact).abs() < 3.0 * standard_error);
    assert!(
        results
            .iter()
            .any(|result| result.off_diagonal.accepted > 0)
    );
}

#[test]
fn one_site_rydberg_boundary_updates_match_exact_thermal_energy() {
    let detuning = 0.9;
    let beta: f64 = 1.1;
    let model = LocalSseModel::rydberg(1, &[detuning], &[], 0.0).unwrap();
    let state = BasisSseState::new(vec![BasisBit::Zero], vec![Operator::identity(); 64]).unwrap();
    let results = run_parallel_chains(
        model,
        state,
        beta,
        SimulationConfig {
            thermalization_sweeps: 500,
            measurement_sweeps: 2_000,
            sweeps_per_measurement: 1,
        },
        123,
        8,
        3,
    )
    .unwrap();
    let (energy, standard_error) = independent_energy_estimate(&results);
    let exact = -detuning * (beta * detuning).exp() / (1.0 + (beta * detuning).exp());
    assert!((energy - exact).abs() < 3.0 * standard_error);
    assert!(results.iter().any(|result| result.basis.accepted > 0));
}

#[test]
fn two_site_tfim_diagonal_energy_matches_exact_ising_limit() {
    let beta: f64 = 0.8;
    let model = LocalSseModel::tfim(2, &[(0, 1)], 1.0, 0.0).unwrap();
    let state = BasisSseState::new(
        vec![BasisBit::Zero, BasisBit::One],
        vec![Operator::identity(); 64],
    )
    .unwrap();
    let results = run_parallel_chains(
        model,
        state,
        beta,
        SimulationConfig {
            thermalization_sweeps: 500,
            measurement_sweeps: 2_000,
            sweeps_per_measurement: 1,
        },
        321,
        8,
        3,
    )
    .unwrap();
    let (energy, standard_error) = independent_energy_estimate(&results);
    let exact = -beta.tanh();
    assert!((energy - exact).abs() < 3.0 * standard_error);
}

#[test]
fn legacy_spin_adapters_are_model_explicit() {
    let spins = [LegacySpin::Up, LegacySpin::Down];
    assert_eq!(
        convert_legacy_bits(LegacyModelKind::Tfim, &spins),
        vec![BasisBit::Zero, BasisBit::One]
    );
    assert_eq!(
        convert_legacy_bits(LegacyModelKind::Rydberg, &spins),
        vec![BasisBit::One, BasisBit::Zero]
    );
}

#[test]
fn logical_chain_seeds_are_stable_and_partition_independent() {
    let all = logical_chain_seeds(77, 8);
    assert_eq!(all[3], derive_chain_seed(77, 3));
    let first = logical_chain_seeds(77, 4);
    let second = (4..8)
        .map(|index| derive_chain_seed(77, index))
        .collect::<Vec<_>>();
    assert_eq!([first, second].concat(), all);
    assert_ne!(all[0], all[1]);
}

#[test]
fn logical_chain_results_are_independent_of_worker_count() {
    let model = LocalSseModel::tfim(1, &[], 0.0, 0.3).unwrap();
    let state = BasisSseState::new(vec![BasisBit::Zero], vec![Operator::identity(); 32]).unwrap();
    let config = SimulationConfig {
        thermalization_sweeps: 10,
        measurement_sweeps: 20,
        sweeps_per_measurement: 1,
    };
    let serial = run_parallel_chains(model.clone(), state.clone(), 0.7, config, 42, 4, 1).unwrap();
    let parallel = run_parallel_chains(model, state, 0.7, config, 42, 4, 3).unwrap();
    assert_eq!(serial, parallel);
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
    let mut sampler = SseSampler::new(model, state, 0.8, ChaCha20Rng::seed_from_u64(7)).unwrap();
    assert!(sampler.diagonal_sweep().is_ok());
    assert!(sampler.state().validate_trace(sampler.model()).is_ok());
}

#[test]
fn operator_string_growth_preserves_existing_vertices() {
    let mut state = BasisSseState::new(
        vec![BasisBit::Zero],
        vec![Operator::off_diagonal(1), Operator::off_diagonal(1)],
    )
    .unwrap();
    state.grow_operator_string(8);
    assert_eq!(
        state.operator_string()[..2],
        [Operator::off_diagonal(1), Operator::off_diagonal(1)]
    );
    assert_eq!(state.operator_string().len(), 8);
}

#[test]
fn sampler_grows_cutoff_before_identity_headroom_is_exhausted() {
    let model = LocalSseModel::tfim(1, &[], 0.0, 0.5).unwrap();
    let mut operators = vec![Operator::identity(); 16];
    operators[0] = Operator::diagonal(0);
    let state = BasisSseState::new(vec![BasisBit::Zero], operators).unwrap();
    let mut sampler = SseSampler::new(model, state, 1.0, ChaCha20Rng::seed_from_u64(8)).unwrap();
    assert!(sampler.ensure_operator_headroom());
    assert_eq!(sampler.state().operator_string().len(), 32);
}

#[test]
fn measurement_phase_rejects_an_exhausted_cutoff_instead_of_mixing_ensembles() {
    let model = LocalSseModel::tfim(1, &[], 0.0, 0.5).unwrap();
    let mut operators = vec![Operator::identity(); 16];
    operators[0] = Operator::diagonal(0);
    let state = BasisSseState::new(vec![BasisBit::Zero], operators).unwrap();
    let mut sampler = SseSampler::new(model, state, 1.0, ChaCha20Rng::seed_from_u64(9)).unwrap();
    let error = sampler
        .run(SimulationConfig {
            thermalization_sweeps: 0,
            measurement_sweeps: 1,
            sweeps_per_measurement: 1,
        })
        .unwrap_err();
    assert!(matches!(
        error,
        qslib_sse::SamplerError::InvalidConfig("operator_string_cutoff")
    ));
}
