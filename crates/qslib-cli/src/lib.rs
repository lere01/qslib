//! Command-line orchestration for physicists using qslib.

use qslib::exact::{self, DenseMatrix, ExactBasis};
use qslib::sse::{BasisSseState, LocalSseModel, Operator, SimulationConfig, run_parallel_chains};
use qslib::{
    BasisBit, Boundary, DenseCouplings, InteractionChannel, InteractionTable, ModelSpecification,
    RectangularGeometry, ResolvedModel, SimulationBasis, SiteCount,
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::fmt::{self, Display, Formatter};
use std::fs;
use std::path::Path;

/// Versioned schema for the compact physicist-facing model input envelope.
pub const MODEL_INPUT_SCHEMA_VERSION: &str = "qslib-model-input-v1";

/// An error reported by the command-line boundary.
#[derive(Debug)]
pub struct CliError(String);

impl CliError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl Display for CliError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for CliError {}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ConventionInput {
    schema: String,
    site_order: String,
    byte_order: String,
    basis: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ModelInput {
    schema_version: String,
    conventions: ConventionInput,
    model: String,
    #[serde(default)]
    couplings: Option<Vec<Vec<f64>>>,
    #[serde(default)]
    fields: Option<Vec<f64>>,
    #[serde(default)]
    omega: Option<Vec<f64>>,
    #[serde(default)]
    detuning: Option<Vec<f64>>,
    #[serde(default)]
    basis: Option<String>,
    #[serde(default)]
    lx: Option<usize>,
    #[serde(default)]
    ly: Option<usize>,
    #[serde(default)]
    boundary_x: Option<String>,
    #[serde(default)]
    boundary_y: Option<String>,
    #[serde(default)]
    j1: Option<f64>,
    #[serde(default)]
    j2: Option<f64>,
    #[serde(default)]
    beta: Option<f64>,
    #[serde(default)]
    thermalization_sweeps: Option<usize>,
    #[serde(default)]
    measurement_sweeps: Option<usize>,
    #[serde(default)]
    sweeps_per_measurement: Option<usize>,
    #[serde(default)]
    chains: Option<usize>,
    #[serde(default)]
    workers: Option<usize>,
    #[serde(default)]
    seed: Option<u64>,
    #[serde(default)]
    operator_cutoff: Option<usize>,
}

/// Execute one qslib command and return either human-readable or JSON output.
///
/// ```
/// let output = qslib_cli::run(&["inspect".into(), "conventions".into()])?;
/// assert!(output.contains("site_order"));
/// # Ok::<(), qslib_cli::CliError>(())
/// ```
pub fn run(args: &[String]) -> Result<String, CliError> {
    if args.is_empty() || args.iter().any(|arg| arg == "--help" || arg == "-h") {
        return Ok(usage());
    }
    let json_output = args.iter().any(|arg| arg == "--json");
    validate_command_args(args)?;
    let value = match args[0].as_str() {
        "inspect" => inspect_command(args)?,
        "model" => model_command(args)?,
        "exact" => exact_command(args)?,
        "sse" => sse_command(args)?,
        "artifacts" => artifact_command(args)?,
        "conformance" => conformance_command(args)?,
        command => {
            return Err(CliError::new(format!(
                "unknown command {command:?}; use --help"
            )));
        }
    };
    if json_output {
        serde_json::to_string_pretty(&value).map_err(|error| CliError::new(error.to_string()))
    } else {
        Ok(human_output(&value))
    }
}

fn usage() -> String {
    "qslib - quantum simulation tools\n\n\
Commands:\n\
  inspect conventions [--json]\n\
  inspect environment [--json]\n\
  model validate CONFIG [--json]\n\
  exact ground-state CONFIG [--json]\n\
  exact evolve CONFIG --t-max TIME [--imaginary] [--json]\n\
  sse run CONFIG [--json]\n\
  artifacts inspect PATH [--json]\n\
  conformance self-test [--json]\n\n\
CONFIG is YAML or JSON and uses canonical row-major sites and explicit physical axes.\n"
        .to_owned()
}

fn validate_command_args(args: &[String]) -> Result<(), CliError> {
    let command = args.first().map(String::as_str).unwrap_or_default();
    let expected_positionals = match (command, args.get(1).map(String::as_str)) {
        ("inspect", Some("conventions" | "environment")) | ("conformance", Some("self-test")) => 2,
        ("model", Some("validate"))
        | ("exact", Some("ground-state"))
        | ("sse", Some("run"))
        | ("artifacts", Some("inspect")) => 3,
        ("exact", Some("evolve")) => 3,
        _ => return Ok(()),
    };
    let mut positional = 2_usize;
    let mut index = 2_usize;
    while index < args.len() {
        match args[index].as_str() {
            "--json" => index += 1,
            "--imaginary"
                if command == "exact" && args.get(1).map(String::as_str) == Some("evolve") =>
            {
                index += 1;
            }
            "--t-max"
                if command == "exact" && args.get(1).map(String::as_str) == Some("evolve") =>
            {
                if args.get(index + 1).is_none() || args[index + 1].starts_with("--") {
                    return Err(CliError::new("--t-max requires one numeric value"));
                }
                index += 2;
            }
            value if value.starts_with('-') => {
                return Err(CliError::new(format!("unknown command option {value:?}")));
            }
            _ => {
                positional += 1;
                index += 1;
            }
        }
    }
    if positional != expected_positionals {
        return Err(CliError::new("unexpected or missing command argument"));
    }
    if command == "exact"
        && args.get(1).map(String::as_str) == Some("evolve")
        && option_value(args, "--t-max").is_none()
    {
        return Err(CliError::new("exact evolve requires --t-max TIME"));
    }
    Ok(())
}

fn inspect_command(args: &[String]) -> Result<Value, CliError> {
    match args.get(1).map(String::as_str) {
        Some("conventions") => Ok(json!({
            "schema": "qslib-conventions-v1",
            "site_order": "row_major",
            "packed_bits": "little_endian_site_zero_is_least_significant",
            "basis_bit_zero": "+1 eigenvalue of the simulation axis",
            "physical_axes": ["x", "y", "z"],
            "energy_normalization": "total unless a field is named density or per_site"
        })),
        Some("environment") => Ok(json!({
            "package": "qslib-quantum",
            "version": env!("CARGO_PKG_VERSION"),
            "rust_version": runtime_rust_version(),
            "features": ["exact", "sse", "io"]
        })),
        Some(command) => Err(CliError::new(format!(
            "unknown inspect command {command:?}"
        ))),
        None => Err(CliError::new("inspect needs conventions or environment")),
    }
}

fn runtime_rust_version() -> String {
    std::process::Command::new("rustc")
        .arg("--version")
        .output()
        .ok()
        .filter(|output| output.status.success())
        .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_owned())
        .filter(|version| !version.is_empty())
        .unwrap_or_else(|| "unavailable".to_owned())
}

fn model_command(args: &[String]) -> Result<Value, CliError> {
    if args.get(1).map(String::as_str) != Some("validate") {
        return Err(CliError::new(
            "model currently supports only validate CONFIG",
        ));
    }
    let input = read_input(args.get(2))?;
    validate_input_contract(&input)?;
    let model = resolve_model(&input)?;
    Ok(model_summary(&input, &model))
}

fn exact_command(args: &[String]) -> Result<Value, CliError> {
    let subcommand = args.get(1).map(String::as_str);
    let input = read_input(args.get(2))?;
    validate_input_contract(&input)?;
    let model = resolve_model(&input)?;
    let matrix = exact_matrix(&model)?;
    match subcommand {
        Some("ground-state") => {
            let spectrum = exact::diagonalize_hermitian(&matrix).map_err(core_error)?;
            let ground = exact::GroundState::from_spectrum(&spectrum).map_err(core_error)?;
            Ok(json!({
                "model": model.family(),
                "basis": model.basis().to_string(),
                "site_count": model.hamiltonian().site_count().get(),
                "dimension": matrix.dimension(),
                "energy": ground.energy(),
                "residual": ground.residual(),
                "normalization": "total energy",
                "provenance": provenance(&input, &model)
            }))
        }
        Some("evolve") => {
            let time = option_value(args, "--t-max")
                .ok_or_else(|| CliError::new("exact evolve requires --t-max TIME"))?
                .parse::<f64>()
                .map_err(|_| CliError::new("--t-max must be finite numeric time"))?;
            let imaginary = args.iter().any(|arg| arg == "--imaginary");
            let mut initial = vec![qslib::Complex64::new(0.0, 0.0); matrix.dimension()];
            initial[0] = qslib::Complex64::new(1.0, 0.0);
            let evolved = exact::evolve(&matrix, &initial, time, imaginary).map_err(core_error)?;
            let norm = evolved.iter().map(|value| value.norm_sqr()).sum::<f64>();
            let energy = exact::expectation(&matrix, &evolved).map_err(core_error)?;
            Ok(json!({
                "model": model.family(),
                "time": time,
                "imaginary_time": imaginary,
                "norm": norm,
                "energy": [energy.re, energy.im],
                "provenance": provenance(&input, &model)
            }))
        }
        Some(command) => Err(CliError::new(format!("unknown exact command {command:?}"))),
        None => Err(CliError::new("exact needs ground-state or evolve")),
    }
}

fn sse_command(args: &[String]) -> Result<Value, CliError> {
    if args.get(1).map(String::as_str) != Some("run") {
        return Err(CliError::new("sse currently supports only run CONFIG"));
    }
    let input = read_input(args.get(2))?;
    validate_input_contract(&input)?;
    if !input.model.eq_ignore_ascii_case("tfim") {
        return Err(CliError::new(
            "SSE CLI currently supports the TFIM decomposition",
        ));
    }
    let couplings = dense_from_input(&input)?;
    let fields = input
        .fields
        .clone()
        .ok_or_else(|| CliError::new("SSE TFIM requires fields"))?;
    let n = couplings.site_count().get();
    if fields.len() != n {
        return Err(CliError::new(format!(
            "field length {}/{} does not match sites",
            fields.len(),
            n
        )));
    }
    let mut bonds = Vec::new();
    let coupling_rows = input
        .couplings
        .as_ref()
        .ok_or_else(|| CliError::new("couplings must be supplied as a square matrix"))?;
    for (first, row) in coupling_rows.iter().enumerate() {
        for (second, value) in row.iter().copied().enumerate().skip(first + 1) {
            if value != 0.0 {
                bonds.push((first as u32, second as u32, value));
            }
        }
    }
    let model = LocalSseModel::tfim_weighted(n, &bonds, &fields).map_err(core_error)?;
    let state = BasisSseState::new(
        vec![BasisBit::Zero; n],
        vec![Operator::identity(); checked_operator_cutoff(&input)?],
    )
    .map_err(core_error)?;
    let results = run_parallel_chains(
        model,
        state,
        input.beta.unwrap_or(1.0),
        SimulationConfig {
            thermalization_sweeps: input.thermalization_sweeps.unwrap_or(100),
            measurement_sweeps: input.measurement_sweeps.unwrap_or(100),
            sweeps_per_measurement: input.sweeps_per_measurement.unwrap_or(1),
        },
        input.seed.unwrap_or(0),
        input.chains.unwrap_or(1),
        input.workers.unwrap_or(1),
    )
    .map_err(core_error)?;
    let energy_per_site = results
        .iter()
        .map(|result| result.thermodynamics.energy_per_site)
        .sum::<f64>()
        / results.len() as f64;
    Ok(json!({
        "model": "tfim",
        "chains": results.len(),
        "energy_per_site": energy_per_site,
        "statistical_note": "independent chain aggregation; inspect autocorrelation before final confidence intervals"
    }))
}

fn artifact_command(args: &[String]) -> Result<Value, CliError> {
    if args.get(1).map(String::as_str) != Some("inspect") {
        return Err(CliError::new(
            "artifacts currently supports only inspect PATH",
        ));
    }
    let path = args
        .get(2)
        .ok_or_else(|| CliError::new("artifact inspect needs PATH"))?;
    let metadata = fs::metadata(path)
        .map_err(|error| CliError::new(format!("cannot inspect {path}: {error}")))?;
    if !metadata.is_dir() {
        return Err(CliError::new(
            "artifact inspection requires a qslib Parquet dataset directory",
        ));
    }
    let manifest =
        qslib::io::ParquetDatasetManifest::inspect(Path::new(path)).map_err(|error| {
            CliError::new(format!("invalid qslib artifact dataset {path}: {error}"))
        })?;
    Ok(json!({
        "path": path,
        "kind": "directory",
        "bytes": metadata.len(),
        "complete_marker": true,
        "complete": manifest.complete,
        "schema_version": manifest.schema_version,
        "convention_schema": manifest.convention_schema,
        "config_checksum": manifest.config_checksum,
        "part_count": manifest.parts.len(),
        "parts": manifest.parts
    }))
}

fn conformance_command(args: &[String]) -> Result<Value, CliError> {
    if args.get(1).map(String::as_str) != Some("self-test") {
        return Err(CliError::new(
            "conformance currently supports only self-test",
        ));
    }
    let model = qslib::tfim(
        &InteractionTable::new(SiteCount::new(1).map_err(core_error)?, Vec::new())
            .map_err(core_error)?,
        &[2.0],
        SimulationBasis::Z,
    )
    .map_err(core_error)?;
    let matrix = exact_matrix(&model)?;
    let expected = [[0.0, -2.0], [-2.0, 0.0]];
    for (row, values) in expected.iter().enumerate() {
        for (column, expected) in values.iter().enumerate() {
            let actual = matrix
                .get(row, column)
                .ok_or_else(|| CliError::new("conformance fixture matrix index overflowed"))?
                .re;
            if (actual - expected).abs() > 1.0e-12 {
                return Err(CliError::new("conformance self-test matrix mismatch"));
            }
        }
    }
    Ok(json!({
        "status": "smoke_pass",
        "scope": "one independent one-site TFIM matrix fixture",
        "fixture_count": 1,
        "fixture": "one-site-tfim-z-basis"
    }))
}

fn read_input(path: Option<&String>) -> Result<ModelInput, CliError> {
    let path =
        path.ok_or_else(|| CliError::new("a YAML or JSON configuration path is required"))?;
    let text = fs::read_to_string(path)
        .map_err(|error| CliError::new(format!("cannot read {path}: {error}")))?;
    serde_yaml_ng::from_str(&text)
        .map_err(|error| CliError::new(format!("invalid configuration {path}: {error}")))
}

fn resolve_model(input: &ModelInput) -> Result<ResolvedModel, CliError> {
    validate_input_contract(input)?;
    let basis = input
        .basis
        .as_deref()
        .unwrap_or("z")
        .parse::<SimulationBasis>()
        .map_err(core_error)?;
    match input.model.to_ascii_lowercase().as_str() {
        "tfim" => {
            let dense = dense_from_input(input)?;
            let fields = input
                .fields
                .as_deref()
                .ok_or_else(|| CliError::new("TFIM requires fields"))?;
            let interactions = dense
                .to_interactions(InteractionChannel::IsingZZ)
                .map_err(core_error)?;
            let table = InteractionTable::new(
                dense.site_count(),
                interactions
                    .into_iter()
                    .map(|term| (term.bond(), term.channel().clone(), term.coefficient()))
                    .collect(),
            )
            .map_err(core_error)?;
            qslib::tfim(&table, fields, basis).map_err(core_error)
        }
        "heisenberg" => {
            let dense = dense_from_input(input)?;
            let interactions = dense
                .to_interactions(InteractionChannel::HeisenbergExchange)
                .map_err(core_error)?;
            let table = InteractionTable::new(
                dense.site_count(),
                interactions
                    .into_iter()
                    .map(|term| (term.bond(), term.channel().clone(), term.coefficient()))
                    .collect(),
            )
            .map_err(core_error)?;
            qslib::heisenberg(&table, basis).map_err(core_error)
        }
        "rydberg" => {
            let dense = dense_from_input(input)?;
            qslib::rydberg(
                &dense,
                input
                    .omega
                    .as_deref()
                    .ok_or_else(|| CliError::new("Rydberg requires omega"))?,
                input
                    .detuning
                    .as_deref()
                    .ok_or_else(|| CliError::new("Rydberg requires detuning"))?,
                basis,
            )
            .map_err(core_error)
        }
        "j1j2" | "j1-j2" => {
            let lx = input.lx.ok_or_else(|| CliError::new("J1-J2 requires lx"))?;
            let ly = input.ly.ok_or_else(|| CliError::new("J1-J2 requires ly"))?;
            let geometry = RectangularGeometry::new(
                lx,
                ly,
                parse_boundary(input.boundary_x.as_deref().unwrap_or("open"))?,
                parse_boundary(input.boundary_y.as_deref().unwrap_or("open"))?,
            )
            .map_err(core_error)?;
            qslib::j1j2(
                &geometry,
                input.j1.ok_or_else(|| CliError::new("J1-J2 requires j1"))?,
                input.j2.ok_or_else(|| CliError::new("J1-J2 requires j2"))?,
                basis,
            )
            .map_err(core_error)
        }
        model => Err(CliError::new(format!("unsupported model family {model:?}"))),
    }
}

fn parse_boundary(value: &str) -> Result<Boundary, CliError> {
    match value {
        "open" => Ok(Boundary::Open),
        "periodic" => Ok(Boundary::Periodic),
        _ => Err(CliError::new(format!(
            "boundary must be open or periodic, got {value:?}"
        ))),
    }
}

fn dense_from_input(input: &ModelInput) -> Result<DenseCouplings, CliError> {
    let rows = input
        .couplings
        .as_ref()
        .ok_or_else(|| CliError::new("couplings must be supplied as a square matrix"))?;
    let n = rows.len();
    if n == 0 || rows.iter().any(|row| row.len() != n) {
        return Err(CliError::new("couplings must be a non-empty square matrix"));
    }
    let values = rows
        .iter()
        .flat_map(|row| row.iter().copied())
        .collect::<Vec<_>>();
    DenseCouplings::new(SiteCount::new(n).map_err(core_error)?, values).map_err(core_error)
}

fn exact_matrix(model: &ResolvedModel) -> Result<DenseMatrix, CliError> {
    let dimension = exact_dimension(model.hamiltonian().site_count().get())?;
    if dimension != 1 {
        let entries = dimension
            .checked_mul(dimension)
            .ok_or_else(|| CliError::new("exact matrix size overflowed"))?;
        let bytes = entries
            .checked_mul(std::mem::size_of::<qslib::Complex64>())
            .ok_or_else(|| CliError::new("exact matrix byte budget overflowed"))?;
        const MAX_EXACT_MATRIX_BYTES: usize = 256 * 1024 * 1024;
        if bytes > MAX_EXACT_MATRIX_BYTES {
            return Err(CliError::new(format!(
                "exact calculation needs {bytes} bytes for a dense matrix; limit is {MAX_EXACT_MATRIX_BYTES} bytes"
            )));
        }
    }
    let basis = ExactBasis::full(model.hamiltonian().site_count()).map_err(core_error)?;
    DenseMatrix::from_hamiltonian(model.hamiltonian(), &basis).map_err(core_error)
}

fn exact_dimension(site_count: usize) -> Result<usize, CliError> {
    let shift = u32::try_from(site_count)
        .map_err(|_| CliError::new("exact calculation site count is too large"))?;
    1_usize
        .checked_shl(shift)
        .ok_or_else(|| CliError::new("exact Hilbert-space dimension overflowed"))
}

fn checked_operator_cutoff(input: &ModelInput) -> Result<usize, CliError> {
    const MAX_OPERATOR_CUTOFF: usize = 1_000_000;
    let cutoff = input.operator_cutoff.unwrap_or(128);
    if cutoff == 0 || cutoff > MAX_OPERATOR_CUTOFF {
        return Err(CliError::new(format!(
            "operator_cutoff must be between 1 and {MAX_OPERATOR_CUTOFF}"
        )));
    }
    Ok(cutoff)
}

fn validate_input_contract(input: &ModelInput) -> Result<(), CliError> {
    if input.schema_version != MODEL_INPUT_SCHEMA_VERSION {
        return Err(CliError::new(format!(
            "unsupported configuration schema {:?}; expected {}",
            input.schema_version, MODEL_INPUT_SCHEMA_VERSION
        )));
    }
    if input.conventions.schema != "qslib-conventions-v1"
        || input.conventions.site_order != "row_major"
        || input.conventions.byte_order != "little_endian"
        || !matches!(input.conventions.basis.as_str(), "x" | "y" | "z")
    {
        return Err(CliError::new(
            "configuration conventions must use qslib-conventions-v1, row_major, little_endian, and basis x, y, or z",
        ));
    }
    let selected_basis = input.basis.as_deref().unwrap_or("z");
    if input.conventions.basis != selected_basis {
        return Err(CliError::new(format!(
            "convention basis {:?} does not match model basis {:?}",
            input.conventions.basis, selected_basis
        )));
    }
    match input.model.to_ascii_lowercase().as_str() {
        "tfim" => reject_fields(
            input,
            &[
                "omega",
                "detuning",
                "j1",
                "j2",
                "lx",
                "ly",
                "boundary_x",
                "boundary_y",
            ],
        )?,
        "heisenberg" => reject_fields(
            input,
            &[
                "fields",
                "omega",
                "detuning",
                "j1",
                "j2",
                "lx",
                "ly",
                "boundary_x",
                "boundary_y",
            ],
        )?,
        "rydberg" => reject_fields(
            input,
            &["fields", "j1", "j2", "lx", "ly", "boundary_x", "boundary_y"],
        )?,
        "j1j2" | "j1-j2" => reject_fields(input, &["fields", "omega", "detuning", "couplings"])?,
        _ => {}
    }
    Ok(())
}

fn reject_fields(input: &ModelInput, names: &[&str]) -> Result<(), CliError> {
    for name in names {
        let present = match *name {
            "omega" => input.omega.is_some(),
            "detuning" => input.detuning.is_some(),
            "fields" => input.fields.is_some(),
            "j1" => input.j1.is_some(),
            "j2" => input.j2.is_some(),
            "lx" => input.lx.is_some(),
            "ly" => input.ly.is_some(),
            "boundary_x" => input.boundary_x.is_some(),
            "boundary_y" => input.boundary_y.is_some(),
            "couplings" => input.couplings.is_some(),
            _ => false,
        };
        if present {
            return Err(CliError::new(format!(
                "configuration field {name:?} is not valid for model {:?}",
                input.model
            )));
        }
    }
    Ok(())
}

fn provenance(input: &ModelInput, model: &ResolvedModel) -> Value {
    json!({
        "schema_version": input.schema_version,
        "convention_schema": input.conventions.schema,
        "site_order": input.conventions.site_order,
        "byte_order": input.conventions.byte_order,
        "simulation_basis": model.basis().to_string(),
        "model_family": model.family(),
        "resolved_interactions": model.interactions().len(),
        "operator_terms": model.hamiltonian().terms().len()
    })
}

fn model_summary(input: &ModelInput, model: &ResolvedModel) -> Value {
    let interactions = model
        .interactions()
        .iter()
        .map(|term| {
            json!({
                "first": term.bond().first().get(),
                "second": term.bond().second().get(),
                "channel": interaction_channel_name(term.channel()),
                "name": term.identity().name(),
                "coefficient": term.coefficient(),
                "image_translation": [term.bond().image_translation().0, term.bond().image_translation().1],
                "source": term.bond().source().get(),
                "direction": [term.bond().direction().0, term.bond().direction().1]
            })
        })
        .collect::<Vec<_>>();
    let specification = match model.specification() {
        ModelSpecification::Tfim { fields } => json!({"fields": fields}),
        ModelSpecification::Heisenberg => json!({}),
        ModelSpecification::Rydberg { omega, detuning } => {
            json!({"omega": omega, "detuning": detuning})
        }
        ModelSpecification::J1J2 { geometry, j1, j2 } => json!({
            "lx": geometry.lx(),
            "ly": geometry.ly(),
            "boundary_x": boundary_name(geometry.boundary_x()),
            "boundary_y": boundary_name(geometry.boundary_y()),
            "j1": j1,
            "j2": j2
        }),
    };
    json!({
        "model": model.family(),
        "basis": model.basis().to_string(),
        "site_count": model.hamiltonian().site_count().get(),
        "interaction_terms": model.interactions().len(),
        "operator_terms": model.hamiltonian().terms().len(),
        "status": "valid",
        "provenance": provenance(input, model),
        "resolved_interactions": interactions,
        "resolved_specification": specification
    })
}

fn interaction_channel_name(channel: &InteractionChannel) -> &'static str {
    match channel {
        InteractionChannel::IsingZZ => "ising_zz",
        InteractionChannel::HeisenbergExchange => "heisenberg_exchange",
        InteractionChannel::RydbergDensityDensity => "rydberg_density_density",
        InteractionChannel::Generic(_) => "generic",
    }
}

fn boundary_name(boundary: Boundary) -> &'static str {
    match boundary {
        Boundary::Open => "open",
        Boundary::Periodic => "periodic",
    }
}

fn option_value(args: &[String], name: &str) -> Option<String> {
    args.windows(2)
        .find(|pair| pair[0] == name)
        .map(|pair| pair[1].clone())
}

fn core_error(error: impl Display) -> CliError {
    CliError::new(error.to_string())
}

fn human_output(value: &Value) -> String {
    match value {
        Value::Object(object) => object
            .iter()
            .map(|(key, value)| format!("{key}: {}", human_value(value)))
            .collect::<Vec<_>>()
            .join("\n"),
        value => human_value(value),
    }
}

fn human_value(value: &Value) -> String {
    match value {
        Value::String(value) => value.clone(),
        _ => value.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn unique_temp_path(prefix: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock must be after the Unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}-{}-{nonce}", std::process::id()))
    }

    #[test]
    fn ground_state_json_is_machine_readable_and_physically_labelled() {
        let path = unique_temp_path("qslib-cli").with_extension("yaml");
        fs::write(
            &path,
            format!(
                "{}model: tfim\ncouplings: [[0.0]]\nfields: [2.0]\nbasis: z\n",
                valid_config_header()
            ),
        )
        .unwrap();
        let output = super::run(&[
            "exact".into(),
            "ground-state".into(),
            path.display().to_string(),
            "--json".into(),
        ])
        .unwrap();
        let value: serde_json::Value = serde_json::from_str(&output).unwrap();
        assert_eq!(value["model"], "tfim");
        assert_eq!(value["site_count"], 1);
        let energy = value["energy"].as_f64().unwrap();
        assert!(
            (energy - (-2.0)).abs() <= 1.0e-12,
            "one-site TFIM ground-state energy drifted: {energy}"
        );
        assert!(value["residual"].as_f64().unwrap() < 1.0e-12);
        let negative_time = super::run(&[
            "exact".into(),
            "evolve".into(),
            path.display().to_string(),
            "--t-max".into(),
            "-0.1".into(),
            "--json".into(),
        ])
        .unwrap();
        assert!(negative_time.contains("imaginary_time"));
        assert!(
            super::run(&["inspect".into(), "conventions".into(), "--imaginary".into()]).is_err()
        );
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn invalid_model_reports_the_physical_field() {
        let path = unique_temp_path("qslib-cli-invalid").with_extension("yaml");
        fs::write(
            &path,
            format!(
                "{}model: tfim\ncouplings: [[0.0]]\nfields: [NaN]\n",
                valid_config_header()
            ),
        )
        .unwrap();
        let error = super::run(&[
            "model".into(),
            "validate".into(),
            path.display().to_string(),
        ])
        .unwrap_err();
        assert!(error.to_string().contains("field"));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn config_requires_versioned_conventions_and_rejects_unknown_fields() {
        let missing_metadata = serde_yaml_ng::from_str::<super::ModelInput>(
            "model: tfim\ncouplings: [[0.0]]\nfields: [2.0]\n",
        )
        .unwrap_err();
        assert!(missing_metadata.to_string().contains("schema_version"));

        let unknown = format!(
            "{}model: tfim\ncouplings: [[0.0]]\nfields: [2.0]\nunexpected: true\n",
            valid_config_header()
        );
        let error = serde_yaml_ng::from_str::<super::ModelInput>(&unknown).unwrap_err();
        assert!(error.to_string().contains("unexpected"));
    }

    #[test]
    fn convention_and_artifact_commands_have_stable_json_fields() {
        let conventions =
            super::run(&["inspect".into(), "conventions".into(), "--json".into()]).unwrap();
        let value: serde_json::Value = serde_json::from_str(&conventions).unwrap();
        assert_eq!(value["site_order"], "row_major");

        let directory = unique_temp_path("qslib-cli-artifact");
        fs::create_dir(&directory).unwrap();
        fs::write(directory.join("COMPLETE"), b"complete\n").unwrap();
        let malformed = super::run(&[
            "artifacts".into(),
            "inspect".into(),
            directory.display().to_string(),
            "--json".into(),
        ]);
        assert!(malformed.is_err(), "a marker is not a complete dataset");

        fs::remove_file(directory.join("COMPLETE")).unwrap();
        let mut manifest = qslib::io::ParquetDatasetManifest::new("a".repeat(64)).unwrap();
        manifest.finish(&directory).unwrap();
        let inspected = super::run(&[
            "artifacts".into(),
            "inspect".into(),
            directory.display().to_string(),
            "--json".into(),
        ])
        .unwrap();
        let value: serde_json::Value = serde_json::from_str(&inspected).unwrap();
        assert_eq!(value["kind"], "directory");
        assert_eq!(value["complete_marker"], true);
        assert_eq!(value["schema_version"], "qslib-parquet-dataset-v1");
        assert_eq!(value["complete"], true);
        assert_eq!(value["part_count"], 0);

        let marker_path = directory.join("COMPLETE");
        let marker = fs::read(&marker_path).unwrap();
        fs::remove_file(&marker_path).unwrap();
        let missing_marker = super::run(&[
            "artifacts".into(),
            "inspect".into(),
            directory.display().to_string(),
            "--json".into(),
        ]);
        assert!(missing_marker.is_err());
        assert!(!marker_path.exists(), "inspection must not repair markers");
        fs::write(&marker_path, marker).unwrap();

        let manifest_path = directory.join("manifest.json");
        let manifest = fs::read_to_string(&manifest_path).unwrap();
        fs::write(
            &manifest_path,
            manifest.replace("qslib-parquet-dataset-v1", "qslib-parquet-dataset-v0"),
        )
        .unwrap();
        let wrong_schema = super::run(&[
            "artifacts".into(),
            "inspect".into(),
            directory.display().to_string(),
            "--json".into(),
        ]);
        assert!(
            wrong_schema.is_err(),
            "unsupported artifact schema must fail"
        );
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn bounded_config_fuzz_never_panics() {
        for seed in 0_u64..1_000 {
            let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15);
            let mut text = String::with_capacity(64);
            for _ in 0..64 {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                let byte = (state >> 56) as u8;
                text.push(char::from(b'a' + (byte % 26)));
            }
            let result =
                std::panic::catch_unwind(|| serde_yaml_ng::from_str::<super::ModelInput>(&text));
            assert!(result.is_ok(), "parser panicked for seed {seed}");
        }
    }

    #[test]
    fn structured_config_fuzz_never_panics_during_resolution() {
        let mut resolved_count = [0_usize; 4];
        for seed in 0_u64..256 {
            let model_index = (seed % 4) as usize;
            let model = match model_index {
                0 => "tfim",
                1 => "heisenberg",
                2 => "rydberg",
                _ => "j1j2",
            };
            let body = match model_index {
                0 => "\"couplings\":[[0.0,1.0],[1.0,0.0]],\"fields\":[0.5,0.5]",
                1 => "\"couplings\":[[0.0,1.0],[1.0,0.0]]",
                2 => {
                    "\"couplings\":[[0.0,1.0],[1.0,0.0]],\"omega\":[1.0,1.0],\"detuning\":[0.2,0.2]"
                }
                _ => "\"lx\":2,\"ly\":2,\"j1\":1.0,\"j2\":0.25",
            };
            let unknown = if (seed / 4) % 2 == 1 {
                format!(",\"unicode\":\"λ🙂{seed}\"")
            } else {
                String::new()
            };
            let payload = format!(
                "{{\"schema_version\":\"qslib-model-input-v1\",\"conventions\":{{\"schema\":\"qslib-conventions-v1\",\"site_order\":\"row_major\",\"byte_order\":\"little_endian\",\"basis\":\"z\"}},\"model\":\"{model}\",{body}{unknown}}}"
            );
            let parsed =
                std::panic::catch_unwind(|| serde_yaml_ng::from_str::<super::ModelInput>(&payload))
                    .expect("structured parser must not panic");
            if let Ok(input) = parsed {
                let resolution = std::panic::catch_unwind(|| super::resolve_model(&input))
                    .expect("model validation must not panic");
                assert!(resolution.is_ok(), "valid structured fixture must resolve");
                resolved_count[model_index] += 1;
            }
        }
        assert_eq!(resolved_count, [32, 32, 32, 32]);
    }

    fn valid_config_header() -> &'static str {
        "schema_version: qslib-model-input-v1\nconventions:\n  schema: qslib-conventions-v1\n  site_order: row_major\n  byte_order: little_endian\n  basis: z\n"
    }
}
