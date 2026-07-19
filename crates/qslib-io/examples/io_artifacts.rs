//! Generate a small durable qslib IO fixture for independent readers.

use qslib_io::{
    AcceptedStateMetadata, Checkpoint, CheckpointArray, ColumnarTrajectory, ConventionMetadata,
    EvolutionControls, ModelMetadata, ParquetDatasetManifest, RngStateMetadata, ScientificConfig,
    SolverMetadata, TrajectoryRow, write_checkpoint_bundle,
};
use std::env;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = PathBuf::from(
        env::args()
            .nth(1)
            .ok_or("usage: io_artifacts <new-output-directory>")?,
    );
    if root.exists() {
        return Err("output directory already exists".into());
    }
    fs::create_dir_all(&root)?;
    let config = ScientificConfig::new(
        "io-fixture",
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
    let accepted = AcceptedStateMetadata::new(
        0.25,
        0.01,
        1,
        0,
        "layout",
        vec![2],
        EvolutionControls {
            seed: 17,
            ..EvolutionControls::default()
        },
    )?;
    let rng = RngStateMetadata::new([4; 32], 1);
    let payload = accepted.to_bytes()?;
    let rng_payload = rng.to_bytes()?;
    let (array, _) = CheckpointArray::from_values("parameters", vec![2], &[1.0, 2.0])?;
    let checkpoint = Checkpoint::new(
        "io-fixture",
        1,
        "layout",
        qslib_io::checksum(&rng_payload),
        &payload,
    )
    .bind_config(&config)?
    .with_state_metadata(accepted, rng)?
    .with_arrays(vec![array.clone()])?;
    write_checkpoint_bundle(
        &root.join("checkpoint"),
        &checkpoint,
        &payload,
        &rng_payload,
        &[(&array, &[1.0, 2.0])],
    )?;
    let trajectory =
        ColumnarTrajectory::new(vec![TrajectoryRow::new(0, 0.0, -1.0)])?.bind_config(&config)?;
    let mut dataset = ParquetDatasetManifest::new(config.checksum()?)?;
    dataset.append_part(&root.join("dataset"), "part-000.parquet", &trajectory)?;
    dataset.write_manifest(&root.join("dataset"))?;
    dataset.finish(&root.join("dataset"))?;
    println!("wrote {}", root.display());
    Ok(())
}
