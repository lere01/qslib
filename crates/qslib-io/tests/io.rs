use qslib_io::{
    AcceptedStateMetadata, ArtifactEntry, ArtifactManifest, Checkpoint, CheckpointBundle,
    ColumnarTrajectory, ConventionMetadata, EvolutionControls, ExactRun, ModelMetadata,
    ParquetDatasetManifest, RngStateMetadata, RunSummary, ScientificConfig, SolverMetadata, SseRun,
    SseRunSpec, TrajectoryRow, atomic_write, checksum,
};
use std::fs;

fn config() -> ScientificConfig {
    let mut model = ModelMetadata::new("tfim", 2, vec![1.0, 0.5, 0.5], "z");
    model.geometry.kind = "chain".into();
    model.interactions = vec![qslib_io::InteractionMetadata {
        identity: "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
        first: 0,
        second: 1,
        channel: "ising_zz".into(),
        name: None,
        coefficient: 1.0,
        multiplicity: 1,
        image_x: 0,
        image_y: 0,
        source: 0,
        direction_x: 0,
        direction_y: 0,
    }];
    model.onsite_terms = vec![
        qslib_io::OnsiteMetadata {
            identity: "field-0".into(),
            site: 0,
            role: "x_field".into(),
            coefficient: 0.5,
        },
        qslib_io::OnsiteMetadata {
            identity: "field-1".into(),
            site: 1,
            role: "x_field".into(),
            coefficient: 0.5,
        },
    ];
    model.resolved_term_ids = vec![
        "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
        "field-0".into(),
        "field-1".into(),
    ];
    ScientificConfig::new(
        "run-1",
        ConventionMetadata::new("qslib-conventions-v1", "row_major", "little_endian", "z"),
        model,
        SolverMetadata::new(
            "dense",
            "f64",
            "chacha20",
            "qslib-seed-v1",
            [7; 32],
            vec![("tolerance".to_owned(), 1.0e-8)],
        ),
    )
}

fn checkpoint_payloads(
    step: u64,
    layout: &str,
) -> (Vec<u8>, Vec<u8>, AcceptedStateMetadata, RngStateMetadata) {
    let accepted = AcceptedStateMetadata::new(
        0.25,
        0.01,
        step,
        2,
        layout,
        vec![2],
        EvolutionControls {
            seed: 17,
            ..EvolutionControls::default()
        },
    )
    .unwrap();
    let rng = RngStateMetadata::new([9; 32], step);
    (
        accepted.to_bytes().unwrap(),
        rng.to_bytes().unwrap(),
        accepted,
        rng,
    )
}

#[test]
fn strict_json_and_yaml_round_trip_preserves_provenance() {
    let value = config();
    assert_eq!(
        ScientificConfig::from_json(&value.to_json().unwrap()).unwrap(),
        value
    );
    assert_eq!(
        ScientificConfig::from_yaml(&value.to_yaml().unwrap()).unwrap(),
        value
    );
    assert!(
        ScientificConfig::from_json(
            r#"{"schema_version":"qslib-config-v1","run_id":"x","extra":1}"#
        )
        .is_err()
    );
    let original = value.model.to_resolved_model().unwrap();
    let reloaded = ScientificConfig::from_json(&value.to_json().unwrap()).unwrap();
    let reconstructed = reloaded.model.to_resolved_model().unwrap();
    assert_eq!(original.hamiltonian(), reconstructed.hamiltonian());
}

#[test]
fn resolved_model_reconstruction_preserves_multiplicity_and_constant() {
    let mut value = config();
    value.model.physical_constant = 3.25;
    value.model.interactions[0].multiplicity = 2;
    let model = value.model.to_resolved_model().unwrap();
    assert_eq!(model.hamiltonian().constant().re, 3.25);
    assert_eq!(model.interactions()[0].coefficient(), 2.0);
}

#[test]
fn duplicate_onsite_roles_on_one_site_are_rejected() {
    let mut value = config();
    value.model.onsite_terms.push(qslib_io::OnsiteMetadata {
        identity: "field-0-duplicate".into(),
        site: 0,
        role: "x_field".into(),
        coefficient: 0.25,
    });
    value
        .model
        .resolved_term_ids
        .push("field-0-duplicate".into());
    value.model.resolved_coefficients.push(0.25);
    assert!(value.validate_semantics().is_err());
}

#[test]
fn positional_coefficients_without_a_typed_term_table_are_rejected() {
    let invalid = ScientificConfig::new(
        "invalid",
        ConventionMetadata::new("qslib-conventions-v1", "row_major", "little_endian", "z"),
        ModelMetadata::new("tfim", 2, vec![1.0], "z"),
        SolverMetadata::new(
            "dense",
            "f64",
            "chacha20",
            "qslib-seed-v1",
            [1; 32],
            vec![("tolerance".into(), 1.0e-8)],
        ),
    );
    assert!(invalid.validate_semantics().is_err());
    assert!(invalid.to_json().is_err());
}

#[test]
fn restart_metadata_rejects_out_of_bounds_steps_and_unknown_seed_versions() {
    let (_, _, accepted, _) = checkpoint_payloads(4, "layout");
    let mut invalid_step = accepted.clone();
    invalid_step.next_step = invalid_step.evolution.dt_max * 2.0;
    assert!(invalid_step.to_bytes().is_err());

    let mut invalid_seed_version = accepted;
    invalid_seed_version.evolution.seed_algorithm_version = 2;
    assert!(invalid_seed_version.to_bytes().is_err());
}

#[test]
fn schema_versions_and_column_lengths_are_checked() {
    let mut json = config().to_json().unwrap();
    json = json.replace("qslib-config-v1", "qslib-config-v0");
    let error = ScientificConfig::from_json(&json).unwrap_err().to_string();
    assert!(error.contains("qslib-config-v1"));
    let legacy = ScientificConfig::from_json(r#"{"run_id":"legacy"}"#).unwrap_err();
    assert!(legacy.to_string().contains("compatibility adapter"));
    let table = ColumnarTrajectory::new(vec![
        TrajectoryRow::new(0, 0.0, -1.0),
        TrajectoryRow::new(1, 0.1, -0.9),
    ])
    .unwrap();
    assert_eq!(table.len(), 2);
    let bound_table = table.clone().bind_config(&config()).unwrap();
    assert_eq!(
        ColumnarTrajectory::from_json(&bound_table.to_json().unwrap()).unwrap(),
        bound_table
    );
    assert!(ColumnarTrajectory::from_columns(vec![0, 1], vec![0.0], vec![-1.0, -0.9]).is_err());
}

#[test]
fn role_specific_run_schemas_reject_invalid_controls() {
    let sse = SseRun::new(SseRunSpec {
        num_sites: 2,
        beta: 1.0,
        operator_string_length: 64,
        thermalization_sweeps: 10,
        measurement_sweeps: 20,
        sweeps_per_measurement: 5,
        chains: 2,
        threads: 2,
        initial_bits: vec![0, 1],
    })
    .unwrap()
    .bind_config(&config())
    .unwrap();
    assert_eq!(SseRun::from_yaml(&sse.to_yaml().unwrap()).unwrap(), sse);
    let invalid = |num_sites, beta, sweeps_per_measurement, initial_bits| {
        SseRun::new(SseRunSpec {
            num_sites,
            beta,
            operator_string_length: 64,
            thermalization_sweeps: 10,
            measurement_sweeps: 20,
            sweeps_per_measurement,
            chains: 2,
            threads: 2,
            initial_bits,
        })
    };
    assert!(invalid(2, 0.0, 5, vec![0, 1]).is_err());
    assert!(invalid(3, 1.0, 5, vec![0, 1]).is_err());
    assert!(invalid(2, 1.0, 0, vec![0, 1]).is_err());
    let exact = ExactRun::new("dense", 1.0e-10, None)
        .unwrap()
        .bind_config(&config())
        .unwrap();
    assert_eq!(
        ExactRun::from_json(&exact.to_json().unwrap()).unwrap(),
        exact
    );
    assert!(ExactRun::new("", 1.0e-10, None).is_err());
    let summary = RunSummary::new(-1.0, 0.01, 20, true)
        .unwrap()
        .bind_config(&config())
        .unwrap();
    assert_eq!(
        RunSummary::from_json(&summary.to_json().unwrap()).unwrap(),
        summary
    );
    assert!(RunSummary::new(-1.0, f64::NAN, 20, true).is_err());
}

#[test]
fn manifest_and_checkpoint_bind_checksums_and_layout() {
    let payload = b"parameters";
    let manifest = ArtifactManifest::new(
        config().checksum().unwrap(),
        vec![ArtifactEntry::new("checkpoint.bin", checksum(payload), 10)],
    )
    .unwrap();
    assert!(manifest.validate_config(&config()).is_ok());
    let manifest_json = manifest.to_json().unwrap();
    assert!(ArtifactManifest::from_json(&manifest_json[..manifest_json.len() - 1]).is_err());
    assert!(
        manifest
            .validate_artifact("checkpoint.bin", payload)
            .is_ok()
    );
    assert!(
        manifest
            .validate_artifact("checkpoint.bin", b"changed")
            .is_err()
    );
    let mut invalid_manifest = manifest.clone();
    invalid_manifest.convention_schema = "legacy-unspecified".into();
    assert!(invalid_manifest.to_json().is_err());
    let (first_array, _) =
        qslib_io::CheckpointArray::from_values("parameters", vec![2], &[1.0, 2.0]).unwrap();
    let (state_payload, rng_payload, accepted, rng) = checkpoint_payloads(4, "layout-hash");
    let checkpoint = Checkpoint::new(
        "run-1",
        4,
        "layout-hash",
        checksum(&rng_payload),
        &state_payload,
    )
    .bind_config(&config())
    .unwrap()
    .with_state_metadata(accepted, rng)
    .unwrap()
    .with_arrays(vec![first_array])
    .unwrap();
    assert!(checkpoint.validate_payload(&state_payload).is_ok());
    assert!(checkpoint.validate_payload(b"changed").is_err());
    let checkpoint_json = checkpoint.to_json().unwrap();
    assert!(Checkpoint::from_json(&checkpoint_json[..checkpoint_json.len() - 1]).is_err());
    let npy = std::env::temp_dir().join(format!("qslib-array-{}.npy", std::process::id()));
    let array = qslib_io::CheckpointArray::new("parameters", "f64", vec![2], "C", &[0; 16]);
    array.write_npy(&npy, &[1.0, 2.0]).unwrap();
    assert_eq!(
        qslib_io::CheckpointArray::read_npy(&npy).unwrap(),
        (vec![1.0, 2.0], vec![2])
    );
    let _ = fs::remove_file(npy);
    let (array, _) =
        qslib_io::CheckpointArray::from_values("parameters", vec![2], &[1.0, 2.0]).unwrap();
    let (state_payload, rng_payload, accepted, rng) = checkpoint_payloads(5, "layout-hash");
    let checkpoint = Checkpoint::new(
        "run-1",
        5,
        "layout-hash",
        checksum(&rng_payload),
        &state_payload,
    )
    .bind_config(&config())
    .unwrap()
    .with_state_metadata(accepted, rng)
    .unwrap()
    .with_arrays(vec![array.clone()])
    .unwrap();
    let bundle = std::env::temp_dir().join(format!("qslib-checkpoint-{}", std::process::id()));
    qslib_io::write_checkpoint_bundle(
        &bundle,
        &checkpoint,
        &state_payload,
        &rng_payload,
        &[(&array, &[1.0, 2.0])],
    )
    .unwrap();
    assert!(bundle.join("checkpoint.json").exists());
    assert!(bundle.join("parameters.npy").exists());
    let loaded: CheckpointBundle = qslib_io::read_checkpoint_bundle(&bundle).unwrap();
    assert_eq!(loaded.checkpoint, checkpoint);
    assert_eq!(loaded.payload, state_payload);
    assert_eq!(loaded.rng_state, rng_payload);
    assert_eq!(loaded.arrays[0].1, vec![1.0, 2.0]);
    let _ = fs::remove_dir_all(bundle);
}

#[test]
fn checkpoint_array_sets_are_exact_and_truncated_payloads_are_rejected() {
    let (first, _) =
        qslib_io::CheckpointArray::from_values("parameters", vec![2], &[1.0, 2.0]).unwrap();
    let (second, _) = qslib_io::CheckpointArray::from_values("extra", vec![1], &[3.0]).unwrap();
    let (state_payload, rng_payload, accepted, rng) = checkpoint_payloads(6, "layout");
    let checkpoint = Checkpoint::new("run-1", 6, "layout", checksum(&rng_payload), &state_payload)
        .bind_config(&config())
        .unwrap()
        .with_state_metadata(accepted, rng)
        .unwrap()
        .with_arrays(vec![first.clone(), second.clone()])
        .unwrap();
    let root = std::env::temp_dir().join(format!("qslib-checkpoint-set-{}", std::process::id()));
    assert!(
        qslib_io::write_checkpoint_bundle(
            &root,
            &checkpoint,
            &state_payload,
            &rng_payload,
            &[(&first, &[1.0, 2.0]), (&first, &[1.0, 2.0])],
        )
        .is_err()
    );
    assert!(!root.exists());
    qslib_io::write_checkpoint_bundle(
        &root,
        &checkpoint,
        &state_payload,
        &rng_payload,
        &[(&first, &[1.0, 2.0]), (&second, &[3.0])],
    )
    .unwrap();
    let npy = root.join("parameters.npy");
    let mut bytes = fs::read(&npy).unwrap();
    bytes.pop();
    fs::write(&npy, bytes).unwrap();
    assert!(qslib_io::read_checkpoint_bundle(&root).is_err());
    let _ = fs::remove_dir_all(root);
}

#[test]
fn atomic_write_replaces_only_the_target_and_round_trips_bytes() {
    let root = std::env::temp_dir().join(format!("qslib-io-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let path = root.join("artifact.json");
    atomic_write(&path, b"first").unwrap();
    atomic_write(&path, b"second").unwrap();
    assert_eq!(fs::read(&path).unwrap(), b"second");
    let _ = fs::remove_dir_all(root);
}

#[test]
fn columnar_trajectory_writes_parquet_parts() {
    let root = std::env::temp_dir().join(format!("qslib-parquet-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let path = root.join("part-000.parquet");
    let table = ColumnarTrajectory::new(vec![TrajectoryRow::new(0, 0.0, -1.0)])
        .unwrap()
        .bind_config(&config())
        .unwrap();
    table.write_parquet(&path).unwrap();
    assert_eq!(ColumnarTrajectory::read_parquet(&path).unwrap(), table);
    let bytes = fs::read(&path).unwrap();
    assert_eq!(&bytes[..4], b"PAR1");
    assert_eq!(&bytes[bytes.len() - 4..], b"PAR1");
    fs::write(&path, &bytes[..bytes.len() - 1]).unwrap();
    assert!(ColumnarTrajectory::read_parquet(&path).is_err());
    let _ = fs::remove_dir_all(root);
}

#[test]
fn parquet_dataset_manifest_is_append_only_and_requires_completion_marker() {
    let root = std::env::temp_dir().join(format!("qslib-dataset-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let mut manifest = ParquetDatasetManifest::new(config().checksum().unwrap()).unwrap();
    let table = ColumnarTrajectory::new(vec![TrajectoryRow::new(0, 0.0, -1.0)])
        .unwrap()
        .bind_config(&config())
        .unwrap();
    manifest
        .append_part(&root, "part-000.parquet", &table)
        .unwrap();
    manifest.write_manifest(&root).unwrap();
    assert!(ParquetDatasetManifest::load(&root).is_err());
    manifest.finish(&root).unwrap();
    let loaded = ParquetDatasetManifest::load(&root).unwrap();
    assert!(loaded.complete);
    assert_eq!(loaded.parts.len(), 1);
    assert_eq!(ParquetDatasetManifest::inspect(&root).unwrap(), loaded);
    assert_eq!(
        ColumnarTrajectory::read_parquet(&root.join("part-000.parquet")).unwrap(),
        table
    );
    assert!(
        manifest
            .append_part(&root, "part-000.parquet", &table)
            .is_err()
    );
    fs::remove_file(root.join(qslib_io::DATASET_COMPLETE_MARKER)).unwrap();
    assert!(ParquetDatasetManifest::inspect(&root).is_err());
    assert!(!root.join(qslib_io::DATASET_COMPLETE_MARKER).exists());
    let recovered = ParquetDatasetManifest::load(&root).unwrap();
    assert!(recovered.complete);
    assert_eq!(
        fs::read(root.join(qslib_io::DATASET_COMPLETE_MARKER)).unwrap(),
        b"qslib-dataset-complete-v1\n"
    );
    fs::write(root.join(qslib_io::DATASET_COMPLETE_MARKER), b"wrong\n").unwrap();
    assert!(ParquetDatasetManifest::inspect(&root).is_err());
    assert_eq!(
        fs::read(root.join(qslib_io::DATASET_COMPLETE_MARKER)).unwrap(),
        b"wrong\n"
    );
    assert!(ParquetDatasetManifest::load(&root).unwrap().complete);
    let _ = fs::remove_dir_all(root);
}

#[test]
fn parquet_finish_validates_parts_before_marking_complete() {
    let root = std::env::temp_dir().join(format!("qslib-dataset-recovery-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let mut manifest = ParquetDatasetManifest::new(config().checksum().unwrap()).unwrap();
    let table = ColumnarTrajectory::new(vec![TrajectoryRow::new(0, 0.0, -1.0)])
        .unwrap()
        .bind_config(&config())
        .unwrap();
    manifest
        .append_part(&root, "part-000.parquet", &table)
        .unwrap();
    fs::remove_file(root.join("part-000.parquet")).unwrap();
    assert!(manifest.finish(&root).is_err());
    assert!(!manifest.complete);
    assert!(!root.join(qslib_io::DATASET_COMPLETE_MARKER).exists());
    let _ = fs::remove_dir_all(root);
}

#[test]
fn model_metadata_retains_interaction_identity_and_realized_disorder() {
    let mut model = ModelMetadata::new("j1j2", 4, vec![1.0, 0.5], "z");
    model.geometry.kind = "chain".into();
    model.geometry.dimensions = vec![4];
    model.geometry.boundaries = vec!["periodic".into()];
    model.interactions = vec![
        qslib_io::InteractionMetadata {
            identity: "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
            first: 0,
            second: 1,
            channel: "ising_zz".into(),
            name: None,
            coefficient: 1.0,
            multiplicity: 1,
            image_x: 0,
            image_y: 0,
            source: 0,
            direction_x: 0,
            direction_y: 0,
        },
        qslib_io::InteractionMetadata {
            identity: "bond:1,2;image:0,0;source:1;direction:0,0;channel:8:ising_zz;name:0:".into(),
            first: 1,
            second: 2,
            channel: "ising_zz".into(),
            name: None,
            coefficient: 0.5,
            multiplicity: 1,
            image_x: 0,
            image_y: 0,
            source: 1,
            direction_x: 0,
            direction_y: 0,
        },
    ];
    model.resolved_term_ids = vec![
        "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
        "bond:1,2;image:0,0;source:1;direction:0,0;channel:8:ising_zz;name:0:".into(),
    ];
    model.disorder = Some(qslib_io::DisorderMetadata {
        seed: [4; 32],
        seed_scheme: "qslib-seed-v1".into(),
        rng_algorithm: "chacha20".into(),
        distribution: "normal".into(),
        realization_index: 7,
        identity_mapping: vec![
            "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
            "bond:1,2;image:0,0;source:1;direction:0,0;channel:8:ising_zz;name:0:".into(),
        ],
        coefficients: vec![1.0, 0.5],
    });
    let mut value = config();
    value.model = model;
    assert!(value.validate_semantics().is_ok());
    let roundtrip = ScientificConfig::from_json(&value.to_json().unwrap()).unwrap();
    assert_eq!(
        roundtrip.model.interactions[0].identity,
        "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:"
    );
    assert_eq!(roundtrip.model.disorder.unwrap().identity_mapping.len(), 2);
}

#[test]
fn resolved_interaction_identity_round_trips_through_core_types() {
    let mut model = ModelMetadata::new("heisenberg", 2, vec![2.0], "z");
    model.geometry.kind = "chain".into();
    model.interactions = vec![qslib_io::InteractionMetadata {
        identity:
            "bond:0,1;image:0,0;source:0;direction:1,0;channel:19:heisenberg_exchange;name:2:j1"
                .into(),
        first: 0,
        second: 1,
        channel: "heisenberg_exchange".into(),
        name: Some("j1".into()),
        coefficient: 2.0,
        multiplicity: 1,
        image_x: 0,
        image_y: 0,
        source: 0,
        direction_x: 1,
        direction_y: 0,
    }];
    model.resolved_term_ids = vec![
        "bond:0,1;image:0,0;source:0;direction:1,0;channel:19:heisenberg_exchange;name:2:j1".into(),
    ];
    let interaction = model.interactions[0].to_weighted_interaction().unwrap();
    let encoded =
        qslib_io::InteractionMetadata::from_weighted_interaction(&interaction, 1).unwrap();
    assert_eq!(encoded.first, 0);
    assert_eq!(encoded.second, 1);
    assert_eq!(encoded.name.as_deref(), Some("j1"));
    assert_eq!(encoded.coefficient, 2.0);
    assert_eq!(encoded.direction_x, 1);
    assert_eq!(encoded.channel, "heisenberg_exchange");
}

#[test]
fn contradictory_disorder_coefficients_are_rejected() {
    let mut model = ModelMetadata::new("tfim", 2, vec![1.0], "z");
    model.geometry.kind = "chain".into();
    model.interactions = vec![qslib_io::InteractionMetadata {
        identity: "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
        first: 0,
        second: 1,
        channel: "ising_zz".into(),
        name: None,
        coefficient: 1.0,
        multiplicity: 1,
        image_x: 0,
        image_y: 0,
        source: 0,
        direction_x: 0,
        direction_y: 0,
    }];
    model.resolved_term_ids =
        vec!["bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into()];
    model.disorder = Some(qslib_io::DisorderMetadata {
        seed: [1; 32],
        seed_scheme: "qslib-seed-v1".into(),
        rng_algorithm: "chacha20".into(),
        distribution: "normal".into(),
        realization_index: 0,
        identity_mapping: vec![
            "bond:0,1;image:0,0;source:0;direction:0,0;channel:8:ising_zz;name:0:".into(),
        ],
        coefficients: vec![2.0],
    });
    let mut value = config();
    value.model = model;
    assert!(value.validate_semantics().is_err());
}

#[test]
fn atomic_write_cleans_failed_temporary_file() {
    let root = std::env::temp_dir().join(format!("qslib-atomic-failure-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let result = atomic_write(&root, b"not a file");
    assert!(result.is_err());
    assert_eq!(fs::read_dir(&root).unwrap().count(), 0);
    let _ = fs::remove_dir_all(root);
}
