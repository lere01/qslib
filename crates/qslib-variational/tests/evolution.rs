use qslib_io::{
    AcceptedStateMetadata, Checkpoint, CheckpointArray, ConventionMetadata, EvolutionControls,
    EvolutionErrorMetric, EvolutionMethod, ModelMetadata, RngStateMetadata, ScientificConfig,
    SolverMetadata, read_checkpoint_bundle, write_checkpoint_bundle,
};
use qslib_variational::{
    DenseQgt, ErrorMetric, EvolutionConfig, EvolutionDriver, EvolutionError, EvolutionMetadata,
    FlatState, IntegrationMethod, Velocity,
};
use std::fs;

fn state(values: &[f64]) -> FlatState {
    FlatState::new("layout-test", values.to_vec()).unwrap()
}

#[test]
fn io_checkpoint_restores_evolution_driver_without_seed_regeneration() {
    let controls = EvolutionControls {
        method: qslib_io::EvolutionMethod::Euler,
        error_metric: qslib_io::EvolutionErrorMetric::Euclidean,
        adaptive: false,
        step_tolerance: 1.0e-6,
        dt_min: 1.0e-8,
        dt_max: 1.0,
        safety_factor: 0.9,
        seed: 17,
        seed_algorithm_version: 1,
    };
    let accepted =
        AcceptedStateMetadata::new(0.2, 0.1, 2, 0, "layout-test", vec![1], controls.clone())
            .unwrap();
    let payload = accepted.to_bytes().unwrap();
    let rng = RngStateMetadata::new([4; 32], 3);
    let rng_payload = rng.to_bytes().unwrap();
    let config = ScientificConfig::new(
        "evolution-restore",
        ConventionMetadata::new("qslib-conventions-v1", "row_major", "little_endian", "z"),
        ModelMetadata::new("empty", 1, Vec::new(), "z"),
        SolverMetadata::new(
            "fixture",
            "f64",
            "chacha20",
            "qslib-seed-v1",
            [3; 32],
            vec![("tolerance".into(), 1.0e-12)],
        ),
    );
    let (array, _) = CheckpointArray::from_values("parameters", vec![1], &[1.2]).unwrap();
    let checkpoint = Checkpoint::new(
        "evolution-restore",
        2,
        "layout-test",
        qslib_io::checksum(&rng_payload),
        &payload,
    )
    .bind_config(&config)
    .unwrap()
    .with_state_metadata(accepted, rng)
    .unwrap()
    .with_arrays(vec![array.clone()])
    .unwrap();
    let root = std::env::temp_dir().join(format!("qslib-evolution-restore-{}", std::process::id()));
    write_checkpoint_bundle(
        &root,
        &checkpoint,
        &payload,
        &rng_payload,
        &[(&array, &[1.2])],
    )
    .unwrap();
    let bundle = read_checkpoint_bundle(&root).unwrap();
    let accepted = AcceptedStateMetadata::from_bytes(&bundle.payload).unwrap();
    assert_eq!(bundle.checkpoint.accepted_step, accepted.accepted_steps);
    assert_eq!(bundle.checkpoint.rng_state.unwrap().word_position(), 48);

    let method = match accepted.evolution.method {
        EvolutionMethod::Euler => IntegrationMethod::Euler,
        EvolutionMethod::Heun => IntegrationMethod::Heun,
    };
    let metric = match accepted.evolution.error_metric {
        EvolutionErrorMetric::Euclidean => ErrorMetric::Euclidean,
        EvolutionErrorMetric::Qgt => ErrorMetric::Qgt,
    };
    let metadata = EvolutionMetadata::from_parts(
        accepted.physical_time,
        accepted.next_step,
        accepted.accepted_steps,
        accepted.rejected_steps,
        accepted.evolution.seed,
        accepted.parameter_layout_fingerprint.clone(),
        method,
        metric,
        accepted.evolution.adaptive,
        accepted.evolution.step_tolerance,
        accepted.evolution.dt_min,
        accepted.evolution.dt_max,
        accepted.evolution.safety_factor,
        accepted.evolution.seed_algorithm_version,
    );
    let evolution_config = EvolutionConfig {
        method,
        error_metric: metric,
        adaptive: accepted.evolution.adaptive,
        dt: accepted.next_step,
        step_tolerance: accepted.evolution.step_tolerance,
        dt_min: accepted.evolution.dt_min,
        dt_max: accepted.evolution.dt_max,
        safety_factor: accepted.evolution.safety_factor,
        seed: accepted.evolution.seed,
    };
    let mut restored =
        EvolutionDriver::from_parts(state(&[1.2]), metadata, evolution_config).unwrap();
    let mut uninterrupted = EvolutionDriver::new(
        state(&[1.0]),
        EvolutionConfig {
            method: IntegrationMethod::Euler,
            dt: 0.1,
            seed: 17,
            ..EvolutionConfig::default()
        },
    )
    .unwrap();
    uninterrupted
        .run(3, |_| {}, |_, _, _| Ok(Velocity::new(vec![1.0])))
        .unwrap();
    restored
        .advance(|_, _, _| Ok(Velocity::new(vec![1.0])))
        .unwrap();
    assert!(
        (restored.state().flat_state().values()[0]
            - uninterrupted.state().flat_state().values()[0])
            .abs()
            < 1.0e-15
    );
    assert!(
        (restored.state().metadata().time() - uninterrupted.state().metadata().time()).abs()
            < 1.0e-15
    );
    let _ = fs::remove_dir_all(root);
}

#[test]
fn euler_and_heun_match_constant_and_linear_odes() {
    let mut euler = EvolutionDriver::new(
        state(&[1.0]),
        EvolutionConfig {
            method: IntegrationMethod::Euler,
            adaptive: false,
            dt: 0.1,
            ..EvolutionConfig::default()
        },
    )
    .unwrap();
    euler
        .advance(|_, _, _| Ok(Velocity::new(vec![2.0])))
        .unwrap();
    assert_eq!(euler.state().flat_state().values(), &[1.2]);

    let mut heun = EvolutionDriver::new(
        state(&[1.0]),
        EvolutionConfig {
            method: IntegrationMethod::Heun,
            adaptive: false,
            dt: 0.1,
            ..EvolutionConfig::default()
        },
    )
    .unwrap();
    heun.advance(|parameters, _, _| Ok(Velocity::new(vec![parameters[0]])))
        .unwrap();
    assert!((heun.state().flat_state().values()[0] - 1.105).abs() < 1.0e-12);
}

#[test]
fn adaptive_rejection_is_transactional_and_observes_only_accepted_states() {
    let mut driver = EvolutionDriver::new(
        state(&[0.0]),
        EvolutionConfig {
            method: IntegrationMethod::Heun,
            adaptive: true,
            dt: 1.0,
            step_tolerance: 0.01,
            dt_min: 0.01,
            dt_max: 1.0,
            safety_factor: 0.5,
            ..EvolutionConfig::default()
        },
    )
    .unwrap();
    let mut observed = Vec::new();
    let mut stage_seeds = Vec::new();
    driver
        .run(
            2,
            |accepted| observed.push(accepted.flat_state().values()[0]),
            |parameters, _, seed| {
                stage_seeds.push(seed);
                Ok(Velocity::new(vec![1.0 + parameters[0]]))
            },
        )
        .unwrap();
    assert_eq!(observed.len(), 2);
    assert_eq!(driver.state().metadata().accepted_steps(), 2);
    assert!(driver.state().metadata().rejected_steps() > 0);
    assert!(driver.state().metadata().proposed_dt() < 1.0);
    for stage in 0..=1 {
        let values = stage_seeds
            .iter()
            .filter(|seed| seed.accepted_step() == 0 && seed.stage() == stage)
            .map(|seed| seed.value())
            .collect::<Vec<_>>();
        assert!(values.len() > 1);
        assert!(values.windows(2).all(|window| window[0] == window[1]));
    }
}

#[test]
fn qgt_error_metric_and_stage_seeds_are_explicit_and_deterministic() {
    let mut driver = EvolutionDriver::new(
        state(&[1.0]),
        EvolutionConfig {
            adaptive: true,
            error_metric: ErrorMetric::Qgt,
            seed: 42,
            dt: 0.1,
            step_tolerance: 0.05,
            ..EvolutionConfig::default()
        },
    )
    .unwrap();
    let qgt = DenseQgt::new(1, vec![4.0]).unwrap();
    let mut seeds = Vec::new();
    let outcome = driver
        .advance(|parameters, _, seed| {
            seeds.push(seed);
            Ok(Velocity::with_qgt(vec![2.0 + parameters[0]], qgt.clone()))
        })
        .unwrap();
    assert!(!seeds.is_empty());
    assert_eq!(seeds[0].master_seed(), 42);
    assert_eq!(seeds[0], seeds[0]);
    assert!((outcome.error_norm() - 0.03).abs() < 1.0e-12);
}

#[test]
fn qgt_metric_rejects_tiny_and_large_indefinite_spectra_scale_relatively() {
    for eigenvalue in [-1.0e-20, -1.0e20] {
        let mut driver = EvolutionDriver::new(
            state(&[1.0]),
            EvolutionConfig {
                adaptive: true,
                error_metric: ErrorMetric::Qgt,
                step_tolerance: 0.1,
                ..EvolutionConfig::default()
            },
        )
        .unwrap();
        let qgt = DenseQgt::new(1, vec![eigenvalue]).unwrap();
        let error = driver
            .advance(|_, _, _| Ok(Velocity::with_qgt(vec![1.0], qgt.clone())))
            .unwrap_err();
        assert!(matches!(
            error,
            qslib_variational::EvolutionError::InvalidParameter(_)
        ));
    }
}

#[test]
fn metadata_round_trip_and_layout_mismatch_are_checked() {
    let driver = EvolutionDriver::new(state(&[1.0, 2.0]), EvolutionConfig::default()).unwrap();
    let encoded = driver.state().metadata().to_json().unwrap();
    let decoded = EvolutionMetadata::from_json(&encoded).unwrap();
    assert_eq!(decoded, *driver.state().metadata());
    let mismatch = FlatState::new("other-layout", vec![1.0, 2.0]).unwrap();
    assert!(EvolutionDriver::from_parts(mismatch, decoded, EvolutionConfig::default()).is_err());
}

#[test]
fn interrupted_and_resumed_trajectory_matches_uninterrupted_run() {
    let config = EvolutionConfig {
        method: IntegrationMethod::Heun,
        adaptive: false,
        dt: 0.05,
        seed: 91,
        ..EvolutionConfig::default()
    };
    let mut uninterrupted = EvolutionDriver::new(state(&[0.25]), config).unwrap();
    uninterrupted
        .run(
            4,
            |_| {},
            |parameters, _, _| Ok(Velocity::new(vec![1.0 + parameters[0]])),
        )
        .unwrap();
    let mut prefix = EvolutionDriver::new(state(&[0.25]), config).unwrap();
    prefix
        .run(
            2,
            |_| {},
            |parameters, _, _| Ok(Velocity::new(vec![1.0 + parameters[0]])),
        )
        .unwrap();
    let metadata =
        EvolutionMetadata::from_json(&prefix.state().metadata().to_json().unwrap()).unwrap();
    let checkpoint_state = FlatState::new(
        prefix.state().flat_state().fingerprint(),
        prefix.state().flat_state().values().to_vec(),
    )
    .unwrap();
    let mut resumed = EvolutionDriver::from_parts(checkpoint_state, metadata, config).unwrap();
    resumed
        .run(
            2,
            |_| {},
            |parameters, _, _| Ok(Velocity::new(vec![1.0 + parameters[0]])),
        )
        .unwrap();
    assert_eq!(resumed.state(), uninterrupted.state());
}

#[test]
fn heun_has_second_order_convergence_on_a_linear_quench() {
    fn integrate(dt: f64) -> f64 {
        let steps = (0.4 / dt).round() as u64;
        let mut driver = EvolutionDriver::new(
            state(&[1.0, 0.0]),
            EvolutionConfig {
                method: IntegrationMethod::Heun,
                adaptive: false,
                dt,
                ..EvolutionConfig::default()
            },
        )
        .unwrap();
        driver
            .run(
                steps,
                |_| {},
                |parameters, _, _| Ok(Velocity::new(vec![-parameters[1], parameters[0]])),
            )
            .unwrap();
        let exact = 0.4_f64.cos();
        (driver.state().flat_state().values()[0] - exact).abs()
    }
    let coarse = integrate(0.1);
    let fine = integrate(0.05);
    assert!(coarse > 0.0 && fine > 0.0);
    assert!(coarse / fine > 3.5);
}

#[test]
fn adaptive_resume_with_rejections_preserves_seed_dependent_trajectory() {
    let config = EvolutionConfig {
        adaptive: true,
        dt: 0.5,
        step_tolerance: 0.005,
        dt_min: 0.001,
        dt_max: 0.5,
        safety_factor: 0.7,
        seed: 177,
        ..EvolutionConfig::default()
    };
    let velocity = |parameters: &[f64], seed: qslib_variational::StageSeed| {
        let noise = (seed.value() % 7) as f64 * 1.0e-3;
        vec![1.0 + parameters[0] + noise]
    };
    let mut uninterrupted = EvolutionDriver::new(state(&[0.0]), config).unwrap();
    let mut uninterrupted_observations = Vec::new();
    let outcomes = uninterrupted
        .run(
            3,
            |accepted| uninterrupted_observations.push(accepted.flat_state().values()[0]),
            |parameters, _, seed| Ok(Velocity::new(velocity(parameters, seed))),
        )
        .unwrap();
    assert!(
        outcomes
            .iter()
            .any(|outcome| outcome.rejected_attempts() > 0)
    );

    let mut prefix = EvolutionDriver::new(state(&[0.0]), config).unwrap();
    let mut prefix_observations = Vec::new();
    prefix
        .run(
            1,
            |accepted| prefix_observations.push(accepted.flat_state().values()[0]),
            |parameters, _, seed| Ok(Velocity::new(velocity(parameters, seed))),
        )
        .unwrap();
    let metadata =
        EvolutionMetadata::from_json(&prefix.state().metadata().to_json().unwrap()).unwrap();
    let checkpoint_state = FlatState::new(
        prefix.state().flat_state().fingerprint(),
        prefix.state().flat_state().values().to_vec(),
    )
    .unwrap();
    let mut resumed = EvolutionDriver::from_parts(checkpoint_state, metadata, config).unwrap();
    let mut resumed_observations = Vec::new();
    let resumed_outcomes = resumed
        .run(
            2,
            |accepted| resumed_observations.push(accepted.flat_state().values()[0]),
            |parameters, _, seed| Ok(Velocity::new(velocity(parameters, seed))),
        )
        .unwrap();
    assert_eq!(resumed_outcomes, outcomes[1..]);
    assert_eq!(prefix_observations, uninterrupted_observations[..1]);
    assert_eq!(resumed_observations, uninterrupted_observations[1..]);
    assert_eq!(resumed.state(), uninterrupted.state());
}

#[test]
fn invalid_callbacks_and_qgt_metric_fail_without_corrupting_state() {
    let mut driver = EvolutionDriver::new(
        state(&[3.0]),
        EvolutionConfig {
            error_metric: ErrorMetric::Qgt,
            adaptive: false,
            ..EvolutionConfig::default()
        },
    )
    .unwrap();
    let before = driver.state().clone();
    let error = driver
        .advance(|_, _, _| Err(EvolutionError::Callback("failed".into())))
        .unwrap_err();
    assert!(matches!(error, EvolutionError::Callback(_)));
    assert_eq!(driver.state(), &before);
}
