//! Versioned scientific configuration, artifacts, and checkpoints.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

use qslib_core::{
    Bond, InteractionChannel, InteractionIdentity, InteractionTable, SimulationBasis, SiteCount,
    SiteId, WeightedInteraction,
};
use serde::{Deserialize, Serialize, de::DeserializeOwned};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::{self, Display, Formatter};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::sync::Arc;

/// Current schema version for scientific run configurations.
pub const CONFIG_SCHEMA_VERSION: &str = "qslib-config-v1";
/// Current schema version for artifact manifests.
pub const MANIFEST_SCHEMA_VERSION: &str = "qslib-manifest-v1";
/// Current schema version for checkpoints.
pub const CHECKPOINT_SCHEMA_VERSION: &str = "qslib-checkpoint-v1";
/// Current schema version for JSON columnar trajectories.
pub const TRAJECTORY_SCHEMA_VERSION: &str = "qslib-trajectory-v1";
/// Current schema version for SSE run controls.
pub const SSE_RUN_SCHEMA_VERSION: &str = "qslib-sse-run-v1";
/// Current schema version for exact run controls.
pub const EXACT_RUN_SCHEMA_VERSION: &str = "qslib-exact-run-v1";
/// Current schema version for run summaries.
pub const SUMMARY_SCHEMA_VERSION: &str = "qslib-summary-v1";
/// Current schema version for append-only Parquet dataset manifests.
pub const PARQUET_DATASET_SCHEMA_VERSION: &str = "qslib-parquet-dataset-v1";
/// Name of the marker written only after a dataset manifest is complete.
pub const DATASET_COMPLETE_MARKER: &str = "COMPLETE";
const CONVENTION_SCHEMA: &str = "qslib-conventions-v1";
const SEED_SCHEME: &str = "qslib-seed-v1";
const RNG_ALGORITHM: &str = "chacha20";
const ACCEPTED_STATE_SCHEMA: &str = "qslib-accepted-state-v1";
const RNG_STATE_SCHEMA: &str = "qslib-rng-state-v1";

/// Errors raised by versioned qslib IO contracts.
#[derive(Debug)]
pub enum IoError {
    /// Serialization or deserialization failed.
    Serialization(String),
    /// A schema version is unsupported.
    UnsupportedVersion {
        /// Required version.
        expected: &'static str,
        /// Received version.
        actual: String,
    },
    /// A recognized unversioned legacy document that requires an explicit adapter.
    LegacyUnsupported {
        /// Human-readable detected format.
        detected: String,
        /// Migration instruction.
        migration: &'static str,
    },
    /// A column or payload has an invalid shape or checksum.
    InvalidData(String),
    /// Filesystem operation failed.
    Filesystem(std::io::Error),
}
impl Display for IoError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Serialization(error) => write!(f, "serialization error: {error}"),
            Self::UnsupportedVersion { expected, actual } => {
                write!(f, "expected schema {expected}, got {actual}")
            }
            Self::LegacyUnsupported {
                detected,
                migration,
            } => {
                write!(f, "legacy document {detected} is unsupported: {migration}")
            }
            Self::InvalidData(error) => f.write_str(error),
            Self::Filesystem(error) => write!(f, "filesystem error: {error}"),
        }
    }
}
impl std::error::Error for IoError {}
impl From<std::io::Error> for IoError {
    fn from(error: std::io::Error) -> Self {
        Self::Filesystem(error)
    }
}

/// Convention metadata persisted with every reconstructible run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConventionMetadata {
    /// Convention schema identifier.
    pub schema: String,
    /// Site ordering.
    pub site_order: String,
    /// Packed-state byte order.
    pub byte_order: String,
    /// Stored simulation basis.
    pub basis: String,
}
impl ConventionMetadata {
    /// Construct convention metadata.
    pub fn new(schema: &str, site_order: &str, byte_order: &str, basis: &str) -> Self {
        Self {
            schema: schema.into(),
            site_order: site_order.into(),
            byte_order: byte_order.into(),
            basis: basis.into(),
        }
    }
}

/// Resolved physical model metadata and coefficients.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelMetadata {
    /// Model family identifier.
    pub kind: String,
    /// Number of sites.
    pub num_sites: usize,
    /// Resolved coefficients in canonical term order.
    pub resolved_coefficients: Vec<f64>,
    /// Stable identity for each entry in `resolved_coefficients`.
    pub resolved_term_ids: Vec<String>,
    /// Stored basis identifier.
    pub basis: String,
    /// Typed resolved geometry and boundary metadata.
    pub geometry: GeometryMetadata,
    /// Resolved pair interaction identities and coefficients.
    pub interactions: Vec<InteractionMetadata>,
    /// Resolved site-local terms.
    pub onsite_terms: Vec<OnsiteMetadata>,
    /// Physical constant contribution.
    pub physical_constant: f64,
    /// Realized disorder provenance and coefficients, if present.
    pub disorder: Option<DisorderMetadata>,
}
/// Reconstructible geometry metadata.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeometryMetadata {
    /// Geometry family.
    pub kind: String,
    /// Dimensions or explicit site count.
    pub dimensions: Vec<usize>,
    /// Boundary condition per dimension.
    pub boundaries: Vec<String>,
    /// Canonical site order.
    pub site_order: String,
    /// Explicit coordinates when the geometry is not generated from dimensions.
    pub coordinates: Vec<Vec<f64>>,
}
/// One resolved pair interaction.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InteractionMetadata {
    /// Stable identity for this resolved interaction.
    pub identity: String,
    /// First endpoint.
    pub first: u32,
    /// Second endpoint.
    pub second: u32,
    /// Canonical operator channel.
    pub channel: String,
    /// Optional physical interaction name, such as `j1` or `j2`.
    pub name: Option<String>,
    /// Resolved coefficient.
    pub coefficient: f64,
    /// Number of physical images represented by this interaction.
    pub multiplicity: usize,
    /// Periodic image translation along x.
    pub image_x: i32,
    /// Periodic image translation along y.
    pub image_y: i32,
    /// Directed source site used to generate this interaction.
    pub source: u32,
    /// Directed generation direction along x.
    pub direction_x: i8,
    /// Directed generation direction along y.
    pub direction_y: i8,
}

impl InteractionMetadata {
    /// Reconstruct the canonical qslib weighted interaction identity.
    pub fn to_weighted_interaction(&self) -> Result<WeightedInteraction, IoError> {
        let bond = Bond::from_parts(
            SiteId::new(self.first),
            SiteId::new(self.second),
            self.image_x,
            self.image_y,
            SiteId::new(self.source),
            self.direction_x,
            self.direction_y,
        )
        .map_err(|error| IoError::InvalidData(format!("invalid bond identity: {error}")))?;
        let channel = match self.channel.as_str() {
            "ising_zz" => InteractionChannel::IsingZZ,
            "heisenberg_exchange" => InteractionChannel::HeisenbergExchange,
            "rydberg_density_density" => InteractionChannel::RydbergDensityDensity,
            value => {
                let generic = value.strip_prefix("generic:").ok_or_else(|| {
                    IoError::InvalidData(format!("unknown interaction channel {value:?}"))
                })?;
                InteractionChannel::generic(generic)
                    .map_err(|error| IoError::InvalidData(format!("invalid channel: {error}")))?
            }
        };
        let identity = match &self.name {
            Some(name) => InteractionIdentity::named(bond, channel, name).map_err(|error| {
                IoError::InvalidData(format!("invalid interaction name: {error}"))
            })?,
            None => InteractionIdentity::new(bond, channel),
        };
        WeightedInteraction::from_identity(identity, self.coefficient).map_err(|error| {
            IoError::InvalidData(format!("invalid interaction coefficient: {error}"))
        })
    }

    /// Encode a canonical qslib interaction as durable identity metadata.
    pub fn from_weighted_interaction(
        interaction: &WeightedInteraction,
        multiplicity: usize,
    ) -> Result<Self, IoError> {
        if multiplicity == 0 {
            return Err(IoError::InvalidData(
                "interaction multiplicity is zero".into(),
            ));
        }
        let bond = interaction.bond();
        let channel = match interaction.channel() {
            InteractionChannel::IsingZZ => "ising_zz".into(),
            InteractionChannel::HeisenbergExchange => "heisenberg_exchange".into(),
            InteractionChannel::RydbergDensityDensity => "rydberg_density_density".into(),
            InteractionChannel::Generic(value) => format!("generic:{value}"),
        };
        let mut encoded = Self {
            identity: String::new(),
            first: bond.first().get(),
            second: bond.second().get(),
            channel,
            name: interaction.identity().name().map(str::to_owned),
            coefficient: interaction.coefficient(),
            multiplicity,
            image_x: bond.image_translation().0,
            image_y: bond.image_translation().1,
            source: bond.source().get(),
            direction_x: bond.direction().0,
            direction_y: bond.direction().1,
        };
        encoded.identity = canonical_interaction_identity(&encoded);
        Ok(encoded)
    }
}
/// One resolved site-local term.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OnsiteMetadata {
    /// Stable identity for this resolved local term.
    pub identity: String,
    /// Site identifier.
    pub site: u32,
    /// Physical role.
    pub role: String,
    /// Resolved coefficient.
    pub coefficient: f64,
}
/// Realized disorder provenance.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DisorderMetadata {
    /// Master seed bytes.
    pub seed: [u8; 32],
    /// Derivation scheme.
    pub seed_scheme: String,
    /// RNG algorithm used to draw the realization.
    pub rng_algorithm: String,
    /// Distribution identifier and parameterization name.
    pub distribution: String,
    /// Index of this realization in an ensemble.
    pub realization_index: u64,
    /// Interaction identities corresponding positionally to `coefficients`.
    pub identity_mapping: Vec<String>,
    /// Realized coefficients in interaction identity order.
    pub coefficients: Vec<f64>,
}
impl ModelMetadata {
    /// Construct model metadata.
    pub fn new(kind: &str, num_sites: usize, resolved_coefficients: Vec<f64>, basis: &str) -> Self {
        let resolved_term_ids = (0..resolved_coefficients.len())
            .map(|index| format!("resolved-{index}"))
            .collect();
        Self {
            kind: kind.into(),
            num_sites,
            resolved_coefficients,
            resolved_term_ids,
            basis: basis.into(),
            geometry: GeometryMetadata {
                kind: "explicit".into(),
                dimensions: vec![num_sites],
                boundaries: vec!["open".into()],
                site_order: "row_major".into(),
                coordinates: Vec::new(),
            },
            interactions: Vec::new(),
            onsite_terms: Vec::new(),
            physical_constant: 0.0,
            disorder: None,
        }
    }

    /// Validate that resolved model metadata is sufficient to reconstruct the
    /// physical term table without regenerating disorder or guessing layout.
    pub fn validate_semantics(&self, site_order: &str, basis: &str) -> Result<(), IoError> {
        let dimensions_product = self
            .geometry
            .dimensions
            .iter()
            .try_fold(1usize, |product, &dimension| product.checked_mul(dimension));
        let mut identities = BTreeSet::new();
        if self.kind.trim().is_empty()
            || self.num_sites == 0
            || self.basis != basis
            || self.geometry.site_order != site_order
            || !matches!(
                self.geometry.kind.as_str(),
                "chain" | "rectangular" | "triangular" | "custom" | "explicit"
            )
            || self.geometry.dimensions.is_empty()
            || self.geometry.dimensions.contains(&0)
            || dimensions_product != Some(self.num_sites)
            || self.geometry.boundaries.len() != self.geometry.dimensions.len()
            || self
                .geometry
                .boundaries
                .iter()
                .any(|boundary| !matches!(boundary.as_str(), "open" | "periodic"))
            || self.geometry.coordinates.iter().any(|coordinate| {
                coordinate.is_empty() || coordinate.iter().any(|value| !value.is_finite())
            })
            || (!self.geometry.coordinates.is_empty()
                && self.geometry.coordinates.len() != self.num_sites)
            || (self.geometry.kind == "custom" && self.geometry.coordinates.len() != self.num_sites)
            || self
                .resolved_coefficients
                .iter()
                .any(|value| !value.is_finite())
            || self.resolved_term_ids.len() != self.resolved_coefficients.len()
            || self
                .resolved_term_ids
                .iter()
                .any(|id| id.trim().is_empty() || !identities.insert(id.clone()))
        {
            return Err(IoError::InvalidData(
                "invalid resolved model metadata".into(),
            ));
        }
        let mut interaction_ids = BTreeSet::new();
        for interaction in &self.interactions {
            if interaction.identity.trim().is_empty()
                || !interaction_ids.insert(interaction.identity.clone())
                || interaction.first >= interaction.second
                || interaction.second as usize >= self.num_sites
                || interaction.source as usize >= self.num_sites
                || (interaction.source != interaction.first
                    && interaction.source != interaction.second)
                || interaction.identity != canonical_interaction_identity(interaction)
                || interaction.channel.trim().is_empty()
                || !is_canonical_channel(&interaction.channel)
                || interaction.name.as_ref().is_some_and(|name| {
                    name.trim().is_empty() || name.chars().any(char::is_whitespace)
                })
                || !interaction.coefficient.is_finite()
                || interaction.multiplicity == 0
            {
                return Err(IoError::InvalidData(
                    "invalid resolved interaction metadata".into(),
                ));
            }
        }
        let mut onsite_ids = BTreeSet::new();
        let mut onsite_roles = BTreeSet::new();
        for onsite in &self.onsite_terms {
            if onsite.identity.trim().is_empty()
                || !onsite_ids.insert(onsite.identity.clone())
                || !onsite_roles.insert((onsite.site, onsite.role.clone()))
                || onsite.site as usize >= self.num_sites
                || onsite.role.trim().is_empty()
                || !onsite.coefficient.is_finite()
            {
                return Err(IoError::InvalidData(
                    "invalid resolved onsite metadata".into(),
                ));
            }
        }
        if (!self.interactions.is_empty() || !self.onsite_terms.is_empty())
            && self.resolved_coefficients.len() != self.interactions.len() + self.onsite_terms.len()
        {
            return Err(IoError::InvalidData(
                "resolved coefficient table does not match typed terms".into(),
            ));
        }
        if self.resolved_coefficients.is_empty()
            != (self.interactions.is_empty() && self.onsite_terms.is_empty())
        {
            return Err(IoError::InvalidData(
                "resolved coefficients require a complete typed term table".into(),
            ));
        }
        if (!self.interactions.is_empty() || !self.onsite_terms.is_empty())
            && self
                .interactions
                .iter()
                .map(|term| term.identity.as_str())
                .chain(self.onsite_terms.iter().map(|term| term.identity.as_str()))
                .any(|identity| !self.resolved_term_ids.iter().any(|id| id == identity))
        {
            return Err(IoError::InvalidData(
                "typed term identity is absent from resolved coefficient table".into(),
            ));
        }
        if !self.interactions.is_empty() || !self.onsite_terms.is_empty() {
            let typed_terms = self
                .interactions
                .iter()
                .map(|term| (&term.identity, term.coefficient))
                .chain(
                    self.onsite_terms
                        .iter()
                        .map(|term| (&term.identity, term.coefficient)),
                )
                .collect::<Vec<_>>();
            if self
                .resolved_term_ids
                .iter()
                .zip(self.resolved_coefficients.iter())
                .zip(typed_terms.iter())
                .any(
                    |((resolved_id, resolved_coefficient), (typed_id, typed_coefficient))| {
                        resolved_id != *typed_id || resolved_coefficient != typed_coefficient
                    },
                )
            {
                return Err(IoError::InvalidData(
                    "resolved terms disagree with typed interaction coefficients".into(),
                ));
            }
        }
        if !self.physical_constant.is_finite() {
            return Err(IoError::InvalidData("non-finite physical constant".into()));
        }
        if let Some(disorder) = &self.disorder {
            let mut mapped = BTreeSet::new();
            if disorder.seed_scheme != SEED_SCHEME
                || disorder.rng_algorithm != RNG_ALGORITHM
                || disorder.distribution.trim().is_empty()
                || disorder.identity_mapping.is_empty()
                || disorder.identity_mapping.len() != disorder.coefficients.len()
                || disorder
                    .identity_mapping
                    .iter()
                    .any(|identity| identity.trim().is_empty() || !mapped.insert(identity.clone()))
                || disorder.coefficients.iter().any(|value| !value.is_finite())
            {
                return Err(IoError::InvalidData("invalid disorder provenance".into()));
            }
            if disorder
                .identity_mapping
                .iter()
                .any(|identity| !interaction_ids.contains(identity))
            {
                return Err(IoError::InvalidData(
                    "disorder identity is not a resolved interaction".into(),
                ));
            }
            for (identity, coefficient) in disorder
                .identity_mapping
                .iter()
                .zip(disorder.coefficients.iter())
            {
                let resolved = self
                    .interactions
                    .iter()
                    .find(|term| &term.identity == identity)
                    .map(|term| term.coefficient);
                if resolved != Some(*coefficient) {
                    return Err(IoError::InvalidData(
                        "disorder coefficient disagrees with resolved interaction".into(),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Reconstruct a canonical qslib model from the resolved term table.
    ///
    /// This adapter accepts only model families whose complete typed terms are
    /// represented by this schema. It never regenerates a coupling from a
    /// seed or infers a missing onsite term.
    pub fn to_resolved_model(&self) -> Result<qslib_core::ResolvedModel, IoError> {
        self.validate_semantics(&self.geometry.site_order, &self.basis)?;
        let basis = self
            .basis
            .parse::<SimulationBasis>()
            .map_err(|error| IoError::InvalidData(format!("invalid simulation basis: {error}")))?;
        let identities = self
            .interactions
            .iter()
            .map(|term| {
                let mut weighted = term.to_weighted_interaction()?;
                let coefficient = weighted.coefficient() * term.multiplicity as f64;
                if !coefficient.is_finite() {
                    return Err(IoError::InvalidData(
                        "interaction multiplicity overflowed its coefficient".into(),
                    ));
                }
                weighted =
                    WeightedInteraction::from_identity(weighted.identity().clone(), coefficient)
                        .map_err(|error| {
                            IoError::InvalidData(format!("invalid interaction: {error}"))
                        })?;
                Ok((weighted.identity().clone(), weighted.coefficient()))
            })
            .collect::<Result<Vec<_>, IoError>>()?;
        let table = InteractionTable::new_with_identities(
            SiteCount::new(self.num_sites)
                .map_err(|error| IoError::InvalidData(format!("invalid site count: {error}")))?,
            identities,
        )
        .map_err(|error| IoError::InvalidData(format!("invalid interaction table: {error}")))?;
        match self.kind.as_str() {
            "tfim" => {
                let mut fields = vec![0.0; self.num_sites];
                for term in &self.onsite_terms {
                    if term.role != "x_field" {
                        return Err(IoError::InvalidData(
                            "TFIM onsite role must be x_field".into(),
                        ));
                    }
                    fields[term.site as usize] = term.coefficient;
                }
                qslib_core::tfim(&table, &fields, basis)
                    .and_then(|model| model.with_physical_constant(self.physical_constant))
                    .map_err(|error| {
                        IoError::InvalidData(format!("TFIM reconstruction failed: {error}"))
                    })
            }
            "heisenberg" => {
                if !self.onsite_terms.is_empty() {
                    return Err(IoError::InvalidData(
                        "Heisenberg reconstruction does not accept onsite terms".into(),
                    ));
                }
                qslib_core::heisenberg(&table, basis)
                    .and_then(|model| model.with_physical_constant(self.physical_constant))
                    .map_err(|error| {
                        IoError::InvalidData(format!("Heisenberg reconstruction failed: {error}"))
                    })
            }
            value => Err(IoError::InvalidData(format!(
                "model family {value:?} has no qslib-io reconstruction adapter"
            ))),
        }
    }
}

/// Solver and reproducibility metadata.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolverMetadata {
    /// Backend identifier.
    pub backend: String,
    /// Scientific scalar dtype.
    pub dtype: String,
    /// RNG algorithm identifier.
    pub rng_algorithm: String,
    /// Seed derivation scheme identifier.
    pub seed_scheme: String,
    /// Canonical 32-byte master seed.
    pub master_seed: [u8; 32],
    /// Named numerical tolerances.
    pub tolerances: BTreeMap<String, f64>,
    /// qslib semantic version.
    pub qslib_version: String,
    /// Source revision used for the run.
    pub software_revision: String,
}
impl SolverMetadata {
    /// Construct solver metadata from named tolerances.
    pub fn new(
        backend: &str,
        dtype: &str,
        rng_algorithm: &str,
        seed_scheme: &str,
        master_seed: [u8; 32],
        tolerances: Vec<(String, f64)>,
    ) -> Self {
        Self {
            backend: backend.into(),
            dtype: dtype.into(),
            rng_algorithm: rng_algorithm.into(),
            seed_scheme: seed_scheme.into(),
            master_seed,
            tolerances: tolerances.into_iter().collect(),
            qslib_version: env!("CARGO_PKG_VERSION").into(),
            software_revision: "unknown".into(),
        }
    }
}

/// Canonical versioned scientific run configuration.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScientificConfig {
    /// Schema version.
    pub schema_version: String,
    /// Stable run identifier.
    pub run_id: String,
    /// Physical convention metadata.
    pub conventions: ConventionMetadata,
    /// Resolved model metadata.
    pub model: ModelMetadata,
    /// Solver and RNG metadata.
    pub solver: SolverMetadata,
}

/// Typed SSE run controls persisted beside a scientific configuration.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SseRun {
    /// Schema version.
    pub schema_version: String,
    /// Checksum of the bound resolved configuration.
    pub config_checksum: String,
    /// Scientific convention schema used by the run.
    pub convention_schema: String,
    /// Number of sites represented by `initial_bits`.
    pub num_sites: usize,
    /// Inverse temperature.
    pub beta: f64,
    /// Initial operator-string cutoff.
    pub operator_string_length: usize,
    /// Thermalization sweeps.
    pub thermalization_sweeps: usize,
    /// Measurement sweeps.
    pub measurement_sweeps: usize,
    /// Sweeps between recorded measurements.
    pub sweeps_per_measurement: usize,
    /// Logical chain count.
    pub chains: usize,
    /// Worker count.
    pub threads: usize,
    /// Canonical initial basis bits.
    pub initial_bits: Vec<u8>,
}
/// Constructor parameters for an SSE run document.
#[derive(Clone, Debug, PartialEq)]
pub struct SseRunSpec {
    /// Number of sites represented by `initial_bits`.
    pub num_sites: usize,
    /// Inverse temperature.
    pub beta: f64,
    /// Initial operator-string cutoff.
    pub operator_string_length: usize,
    /// Thermalization sweeps.
    pub thermalization_sweeps: usize,
    /// Measurement sweeps.
    pub measurement_sweeps: usize,
    /// Sweeps between recorded measurements.
    pub sweeps_per_measurement: usize,
    /// Logical chain count.
    pub chains: usize,
    /// Worker count.
    pub threads: usize,
    /// Canonical initial basis bits.
    pub initial_bits: Vec<u8>,
}
impl SseRun {
    /// Construct validated SSE controls.
    pub fn new(spec: SseRunSpec) -> Result<Self, IoError> {
        if !spec.beta.is_finite()
            || spec.beta <= 0.0
            || spec.num_sites == 0
            || spec.operator_string_length == 0
            || spec.measurement_sweeps == 0
            || spec.sweeps_per_measurement == 0
            || spec.chains == 0
            || spec.threads == 0
            || spec.initial_bits.len() != spec.num_sites
            || spec.initial_bits.iter().any(|&bit| bit > 1)
        {
            return Err(IoError::InvalidData("invalid SSE run controls".into()));
        }
        Ok(Self {
            schema_version: SSE_RUN_SCHEMA_VERSION.into(),
            config_checksum: String::new(),
            convention_schema: String::new(),
            num_sites: spec.num_sites,
            beta: spec.beta,
            operator_string_length: spec.operator_string_length,
            thermalization_sweeps: spec.thermalization_sweeps,
            measurement_sweeps: spec.measurement_sweeps,
            sweeps_per_measurement: spec.sweeps_per_measurement,
            chains: spec.chains,
            threads: spec.threads,
            initial_bits: spec.initial_bits,
        })
    }
    /// Bind controls to a resolved configuration before durable serialization.
    pub fn bind_config(mut self, config: &ScientificConfig) -> Result<Self, IoError> {
        self.config_checksum = config.checksum()?;
        self.convention_schema = config.conventions.schema.clone();
        self.validate_config(config)?;
        Ok(self)
    }
    /// Validate controls and their configuration binding.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.schema_version != SSE_RUN_SCHEMA_VERSION
            || self.config_checksum.len() != 64
            || !is_hex(&self.config_checksum)
            || self.convention_schema != CONVENTION_SCHEMA
            || self.num_sites == 0
            || !self.beta.is_finite()
            || self.beta <= 0.0
            || self.operator_string_length == 0
            || self.measurement_sweeps == 0
            || self.sweeps_per_measurement == 0
            || self.chains == 0
            || self.threads == 0
            || self.initial_bits.len() != self.num_sites
            || self.initial_bits.iter().any(|&bit| bit > 1)
        {
            return Err(IoError::InvalidData("invalid SSE run controls".into()));
        }
        Ok(())
    }
    /// Verify that this role document is bound to the supplied configuration.
    pub fn validate_config(&self, config: &ScientificConfig) -> Result<(), IoError> {
        self.validate_semantics()?;
        if self.config_checksum != config.checksum()?
            || self.convention_schema != config.conventions.schema
        {
            return Err(IoError::InvalidData("SSE/config binding mismatch".into()));
        }
        if self.num_sites != config.model.num_sites {
            return Err(IoError::InvalidData("SSE/model site-count mismatch".into()));
        }
        Ok(())
    }
    /// Serialize as strict JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        strict_json(self)
    }
    /// Deserialize and validate strict JSON.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self = strict_from_json(value)?;
        if parsed.schema_version != SSE_RUN_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: SSE_RUN_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Serialize as strict YAML.
    pub fn to_yaml(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_yaml_ng::to_string(self).map_err(|error| IoError::Serialization(error.to_string()))
    }
    /// Deserialize and validate strict YAML.
    pub fn from_yaml(value: &str) -> Result<Self, IoError> {
        let parsed: Self = serde_yaml_ng::from_str(value)
            .map_err(|error| IoError::Serialization(error.to_string()))?;
        if parsed.schema_version != SSE_RUN_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: SSE_RUN_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
}
/// Typed exact-solver controls.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExactRun {
    /// Schema version.
    pub schema_version: String,
    /// Checksum of the bound resolved configuration.
    pub config_checksum: String,
    /// Scientific convention schema used by the run.
    pub convention_schema: String,
    /// Solver backend.
    pub backend: String,
    /// Residual tolerance.
    pub tolerance: f64,
    /// Optional fixed-sector weight.
    pub sector_weight: Option<usize>,
}
impl ExactRun {
    /// Construct exact-solver controls.
    pub fn new(
        backend: &str,
        tolerance: f64,
        sector_weight: Option<usize>,
    ) -> Result<Self, IoError> {
        if backend.trim().is_empty() || !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(IoError::InvalidData("invalid exact run controls".into()));
        }
        Ok(Self {
            schema_version: EXACT_RUN_SCHEMA_VERSION.into(),
            config_checksum: String::new(),
            convention_schema: String::new(),
            backend: backend.into(),
            tolerance,
            sector_weight,
        })
    }
    /// Bind controls to a resolved configuration before durable serialization.
    pub fn bind_config(mut self, config: &ScientificConfig) -> Result<Self, IoError> {
        self.config_checksum = config.checksum()?;
        self.convention_schema = config.conventions.schema.clone();
        self.validate_config(config)?;
        Ok(self)
    }
    /// Validate exact-solver controls and their configuration binding.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.schema_version != EXACT_RUN_SCHEMA_VERSION
            || self.config_checksum.len() != 64
            || !is_hex(&self.config_checksum)
            || self.convention_schema != CONVENTION_SCHEMA
            || self.backend.trim().is_empty()
            || !self.tolerance.is_finite()
            || self.tolerance <= 0.0
        {
            return Err(IoError::InvalidData("invalid exact run controls".into()));
        }
        Ok(())
    }
    /// Verify that this role document is bound to the supplied configuration.
    pub fn validate_config(&self, config: &ScientificConfig) -> Result<(), IoError> {
        self.validate_semantics()?;
        if self.config_checksum != config.checksum()?
            || self.convention_schema != config.conventions.schema
        {
            return Err(IoError::InvalidData("exact/config binding mismatch".into()));
        }
        if self
            .sector_weight
            .is_some_and(|weight| weight > config.model.num_sites)
        {
            return Err(IoError::InvalidData(
                "exact sector weight exceeds model site count".into(),
            ));
        }
        Ok(())
    }
    /// Serialize as strict JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        strict_json(self)
    }
    /// Deserialize and validate strict JSON.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self = strict_from_json(value)?;
        if parsed.schema_version != EXACT_RUN_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: EXACT_RUN_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Serialize as strict YAML.
    pub fn to_yaml(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_yaml_ng::to_string(self).map_err(|error| IoError::Serialization(error.to_string()))
    }
    /// Deserialize and validate strict YAML.
    pub fn from_yaml(value: &str) -> Result<Self, IoError> {
        let parsed: Self = serde_yaml_ng::from_str(value)
            .map_err(|error| IoError::Serialization(error.to_string()))?;
        if parsed.schema_version != EXACT_RUN_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: EXACT_RUN_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
}
/// Typed summary of a completed run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunSummary {
    /// Schema version.
    pub schema_version: String,
    /// Checksum of the bound resolved configuration.
    pub config_checksum: String,
    /// Scientific convention schema used by the run.
    pub convention_schema: String,
    /// Total-energy estimate.
    pub energy: f64,
    /// Standard error or residual bound.
    pub uncertainty: f64,
    /// Number of accepted samples or solver states.
    pub samples: u64,
    /// Whether all declared artifacts are complete.
    pub complete: bool,
}
impl RunSummary {
    /// Construct a finite run summary.
    pub fn new(
        energy: f64,
        uncertainty: f64,
        samples: u64,
        complete: bool,
    ) -> Result<Self, IoError> {
        if !energy.is_finite() || !uncertainty.is_finite() || uncertainty < 0.0 {
            return Err(IoError::InvalidData("invalid run summary".into()));
        }
        Ok(Self {
            schema_version: SUMMARY_SCHEMA_VERSION.into(),
            config_checksum: String::new(),
            convention_schema: String::new(),
            energy,
            uncertainty,
            samples,
            complete,
        })
    }
    /// Bind the summary to a resolved configuration before durable serialization.
    pub fn bind_config(mut self, config: &ScientificConfig) -> Result<Self, IoError> {
        self.config_checksum = config.checksum()?;
        self.convention_schema = config.conventions.schema.clone();
        self.validate_config(config)?;
        Ok(self)
    }
    /// Validate summary values and their configuration binding.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.schema_version != SUMMARY_SCHEMA_VERSION
            || self.config_checksum.len() != 64
            || !is_hex(&self.config_checksum)
            || self.convention_schema != CONVENTION_SCHEMA
            || !self.energy.is_finite()
            || !self.uncertainty.is_finite()
            || self.uncertainty < 0.0
        {
            return Err(IoError::InvalidData("invalid run summary".into()));
        }
        Ok(())
    }
    /// Verify that this summary is bound to the supplied configuration.
    pub fn validate_config(&self, config: &ScientificConfig) -> Result<(), IoError> {
        self.validate_semantics()?;
        if self.config_checksum != config.checksum()?
            || self.convention_schema != config.conventions.schema
        {
            return Err(IoError::InvalidData(
                "summary/config binding mismatch".into(),
            ));
        }
        Ok(())
    }
    /// Serialize as strict JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        strict_json(self)
    }
    /// Deserialize and validate strict JSON.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self = strict_from_json(value)?;
        if parsed.schema_version != SUMMARY_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: SUMMARY_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Serialize as strict YAML.
    pub fn to_yaml(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_yaml_ng::to_string(self).map_err(|error| IoError::Serialization(error.to_string()))
    }
    /// Deserialize and validate strict YAML.
    pub fn from_yaml(value: &str) -> Result<Self, IoError> {
        let parsed: Self = serde_yaml_ng::from_str(value)
            .map_err(|error| IoError::Serialization(error.to_string()))?;
        if parsed.schema_version != SUMMARY_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: SUMMARY_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
}
impl ScientificConfig {
    /// Construct a current versioned configuration.
    pub fn new(
        run_id: &str,
        conventions: ConventionMetadata,
        model: ModelMetadata,
        solver: SolverMetadata,
    ) -> Self {
        Self {
            schema_version: CONFIG_SCHEMA_VERSION.into(),
            run_id: run_id.into(),
            conventions,
            model,
            solver,
        }
    }
    /// Serialize as strict JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_json::to_string_pretty(self).map_err(|e| IoError::Serialization(e.to_string()))
    }
    /// Deserialize and validate strict JSON.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let raw: serde_json::Value =
            serde_json::from_str(value).map_err(|e| IoError::Serialization(e.to_string()))?;
        if raw.get("schema_version").is_none() {
            return Err(IoError::LegacyUnsupported {
                detected: "unversioned scientific configuration".into(),
                migration: "resolve the legacy input through a named compatibility adapter before loading qslib-config-v1",
            });
        }
        let parsed: Self =
            serde_json::from_value(raw).map_err(|e| IoError::Serialization(e.to_string()))?;
        if parsed.schema_version != CONFIG_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: CONFIG_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Serialize as strict YAML.
    pub fn to_yaml(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_yaml_ng::to_string(self).map_err(|e| IoError::Serialization(e.to_string()))
    }
    /// Deserialize and validate strict YAML.
    pub fn from_yaml(value: &str) -> Result<Self, IoError> {
        let parsed: Self =
            serde_yaml_ng::from_str(value).map_err(|e| IoError::Serialization(e.to_string()))?;
        if parsed.schema_version != CONFIG_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: CONFIG_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Return the BLAKE3 checksum of canonical JSON bytes.
    pub fn checksum(&self) -> Result<String, IoError> {
        Ok(checksum(self.to_json()?.as_bytes()))
    }
    /// Validate physical identifiers, dimensions, and finite numeric values.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.run_id.trim().is_empty()
            || self.conventions.schema != CONVENTION_SCHEMA
            || !matches!(self.conventions.site_order.as_str(), "row_major" | "custom")
            || self.conventions.byte_order != "little_endian"
            || !matches!(self.conventions.basis.as_str(), "x" | "y" | "z")
            || self.solver.backend.trim().is_empty()
            || self.solver.qslib_version.trim().is_empty()
            || self.solver.software_revision.trim().is_empty()
            || self.solver.dtype != "f64"
            || self.solver.rng_algorithm != RNG_ALGORITHM
            || self.solver.seed_scheme != SEED_SCHEME
            || self
                .solver
                .tolerances
                .iter()
                .any(|(name, value)| name.trim().is_empty() || !value.is_finite() || *value < 0.0)
        {
            return Err(IoError::InvalidData(
                "scientific configuration violates canonical semantics".into(),
            ));
        }
        self.model
            .validate_semantics(&self.conventions.site_order, &self.conventions.basis)
    }
}

/// One immutable artifact entry in a manifest.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArtifactEntry {
    /// Relative artifact path.
    pub path: String,
    /// BLAKE3 digest.
    pub checksum: String,
    /// Byte length.
    pub bytes: u64,
}
impl ArtifactEntry {
    /// Construct an artifact entry.
    pub fn new(path: &str, checksum: String, bytes: u64) -> Self {
        Self {
            path: path.into(),
            checksum,
            bytes,
        }
    }
}

/// Atomic manifest describing a reconstructible run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArtifactManifest {
    /// Schema version.
    pub schema_version: String,
    /// Scientific convention schema used by every listed artifact.
    pub convention_schema: String,
    /// Configuration checksum.
    pub config_checksum: String,
    /// Entries.
    pub entries: Vec<ArtifactEntry>,
    /// Completion marker.
    pub complete: bool,
}
impl ArtifactManifest {
    /// Construct a complete manifest candidate.
    pub fn new(config_checksum: String, entries: Vec<ArtifactEntry>) -> Result<Self, IoError> {
        let value = Self {
            schema_version: MANIFEST_SCHEMA_VERSION.into(),
            convention_schema: CONVENTION_SCHEMA.into(),
            config_checksum,
            entries,
            complete: true,
        };
        value.validate_semantics()?;
        Ok(value)
    }
    /// Validate one artifact's bytes against the manifest.
    pub fn validate_artifact(&self, path: &str, bytes: &[u8]) -> Result<(), IoError> {
        let entry = self
            .entries
            .iter()
            .find(|entry| entry.path == path)
            .ok_or_else(|| IoError::InvalidData("artifact is not listed".into()))?;
        if entry.bytes != bytes.len() as u64 || entry.checksum != checksum(bytes) {
            return Err(IoError::InvalidData(
                "artifact checksum or length mismatch".into(),
            ));
        }
        Ok(())
    }
    /// Serialize the manifest as JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_json::to_string_pretty(self).map_err(|e| IoError::Serialization(e.to_string()))
    }
    /// Deserialize and validate a manifest.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self =
            serde_json::from_str(value).map_err(|e| IoError::Serialization(e.to_string()))?;
        if parsed.schema_version != MANIFEST_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: MANIFEST_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Validate schema, checksums, paths, and duplicate identities.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        let mut paths = std::collections::BTreeSet::new();
        if !self.complete
            || self.convention_schema != CONVENTION_SCHEMA
            || self.config_checksum.len() != 64
            || !self
                .config_checksum
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
            || self.entries.iter().any(|entry| {
                !safe_relative_path(&entry.path)
                    || entry.checksum.len() != 64
                    || !entry.checksum.bytes().all(|byte| byte.is_ascii_hexdigit())
                    || !paths.insert(entry.path.clone())
            })
        {
            return Err(IoError::InvalidData("invalid manifest semantics".into()));
        }
        Ok(())
    }
    /// Validate that the manifest is bound to this resolved configuration.
    pub fn validate_config(&self, config: &ScientificConfig) -> Result<(), IoError> {
        if self.convention_schema == config.conventions.schema
            && self.config_checksum == config.checksum()?
        {
            Ok(())
        } else {
            Err(IoError::InvalidData(
                "manifest/config checksum mismatch".into(),
            ))
        }
    }
}

/// Integration method recorded at an accepted evolution boundary.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum EvolutionMethod {
    /// First-order forward Euler.
    Euler,
    /// Second-order explicit Heun method.
    Heun,
}

/// Adaptive error metric recorded at an accepted evolution boundary.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum EvolutionErrorMetric {
    /// Euclidean predictor-corrector error.
    Euclidean,
    /// Quantum-geometric-tensor error.
    Qgt,
}

/// Complete trajectory-changing controls needed to restore an evolution
/// driver without consulting the original configuration or regenerating a
/// seed.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EvolutionControls {
    /// Integration method.
    pub method: EvolutionMethod,
    /// Adaptive error metric.
    pub error_metric: EvolutionErrorMetric,
    /// Whether adaptive acceptance is enabled.
    pub adaptive: bool,
    /// Adaptive error tolerance.
    pub step_tolerance: f64,
    /// Lower step bound.
    pub dt_min: f64,
    /// Upper step bound.
    pub dt_max: f64,
    /// Adaptive safety factor.
    pub safety_factor: f64,
    /// Evolution master seed.
    pub seed: u64,
    /// Version of deterministic stage-seed derivation.
    pub seed_algorithm_version: u32,
}
impl Default for EvolutionControls {
    fn default() -> Self {
        Self {
            method: EvolutionMethod::Heun,
            error_metric: EvolutionErrorMetric::Euclidean,
            adaptive: false,
            step_tolerance: 1.0e-6,
            dt_min: 1.0e-8,
            dt_max: 1.0,
            safety_factor: 0.9,
            seed: 0,
            seed_algorithm_version: 1,
        }
    }
}
impl EvolutionControls {
    /// Validate trajectory-changing controls.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.adaptive && self.method == EvolutionMethod::Euler
            || !self.step_tolerance.is_finite()
            || self.step_tolerance <= 0.0
            || !self.dt_min.is_finite()
            || self.dt_min <= 0.0
            || !self.dt_max.is_finite()
            || self.dt_max < self.dt_min
            || !self.safety_factor.is_finite()
            || self.safety_factor <= 0.0
            || self.safety_factor > 1.0
            || self.seed_algorithm_version != 1
        {
            return Err(IoError::InvalidData("invalid evolution controls".into()));
        }
        Ok(())
    }
}

/// Typed accepted-boundary metadata for restartable evolution.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AcceptedStateMetadata {
    /// Versioned accepted-state schema.
    pub schema_version: String,
    /// Physical time at the accepted boundary.
    pub physical_time: f64,
    /// Proposed next integration step.
    pub next_step: f64,
    /// Number of accepted steps.
    pub accepted_steps: u64,
    /// Number of rejected proposals.
    pub rejected_steps: u64,
    /// Parameter-layout fingerprint shared with the envelope.
    pub parameter_layout_fingerprint: String,
    /// Shape of the named parameter array.
    pub parameter_shape: Vec<usize>,
    /// Parameter dtype.
    pub dtype: String,
    /// Parameter memory order.
    pub order: String,
    /// Complete trajectory-changing integration controls.
    pub evolution: EvolutionControls,
}
impl AcceptedStateMetadata {
    /// Construct typed accepted-boundary metadata.
    pub fn new(
        physical_time: f64,
        next_step: f64,
        accepted_steps: u64,
        rejected_steps: u64,
        parameter_layout_fingerprint: &str,
        parameter_shape: Vec<usize>,
        evolution: EvolutionControls,
    ) -> Result<Self, IoError> {
        let value = Self {
            schema_version: ACCEPTED_STATE_SCHEMA.into(),
            physical_time,
            next_step,
            accepted_steps,
            rejected_steps,
            parameter_layout_fingerprint: parameter_layout_fingerprint.into(),
            parameter_shape,
            dtype: "f64".into(),
            order: "C".into(),
            evolution,
        };
        value.validate_semantics()?;
        Ok(value)
    }
    /// Validate restart metadata.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        let product = self
            .parameter_shape
            .iter()
            .try_fold(1usize, |product, &size| product.checked_mul(size));
        if self.schema_version != ACCEPTED_STATE_SCHEMA
            || !self.physical_time.is_finite()
            || !self.next_step.is_finite()
            || self.next_step <= 0.0
            || self.parameter_layout_fingerprint.trim().is_empty()
            || self.parameter_shape.is_empty()
            || self.parameter_shape.contains(&0)
            || product.is_none()
            || self.dtype != "f64"
            || self.order != "C"
        {
            return Err(IoError::InvalidData(
                "invalid accepted-state metadata".into(),
            ));
        }
        self.evolution.validate_semantics()?;
        if self.next_step < self.evolution.dt_min || self.next_step > self.evolution.dt_max {
            return Err(IoError::InvalidData(
                "accepted next step is outside evolution bounds".into(),
            ));
        }
        Ok(())
    }
    /// Encode the metadata as a portable payload.
    pub fn to_bytes(&self) -> Result<Vec<u8>, IoError> {
        self.validate_semantics()?;
        serde_json::to_vec(self).map_err(|error| IoError::Serialization(error.to_string()))
    }
    /// Decode and validate a portable accepted-state payload.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, IoError> {
        let value: Self = serde_json::from_slice(bytes)
            .map_err(|error| IoError::Serialization(error.to_string()))?;
        value.validate_semantics()?;
        Ok(value)
    }
}

/// Versioned metadata for a serialized ChaCha20 stream state.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RngStateMetadata {
    /// Versioned RNG-state schema.
    pub schema_version: String,
    /// Algorithm identifier.
    pub algorithm: String,
    /// Master seed for the logical stream.
    pub seed: [u8; 32],
    /// Number of complete 64-byte ChaCha20 blocks consumed by the stream.
    pub position_blocks: u64,
}
impl RngStateMetadata {
    /// Construct typed RNG-state metadata from a ChaCha20 block position.
    pub fn new(seed: [u8; 32], position_blocks: u64) -> Self {
        Self {
            schema_version: RNG_STATE_SCHEMA.into(),
            algorithm: RNG_ALGORITHM.into(),
            seed,
            position_blocks,
        }
    }
    /// Return the equivalent little-endian 32-bit word position used by
    /// ChaCha20 stream implementations (16 words per 64-byte block).
    pub fn word_position(&self) -> u128 {
        self.position_blocks as u128 * 16
    }
    /// Validate the algorithm and schema identifier.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.schema_version != RNG_STATE_SCHEMA || self.algorithm != RNG_ALGORITHM {
            return Err(IoError::InvalidData("invalid RNG-state metadata".into()));
        }
        Ok(())
    }
    /// Encode the metadata as a portable payload.
    pub fn to_bytes(&self) -> Result<Vec<u8>, IoError> {
        self.validate_semantics()?;
        serde_json::to_vec(self).map_err(|error| IoError::Serialization(error.to_string()))
    }
    /// Decode and validate a portable RNG-state payload.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, IoError> {
        let value: Self = serde_json::from_slice(bytes)
            .map_err(|error| IoError::Serialization(error.to_string()))?;
        value.validate_semantics()?;
        Ok(value)
    }
}

/// Restartable accepted-boundary checkpoint metadata.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Checkpoint {
    /// Schema version.
    pub schema_version: String,
    /// Scientific convention schema used by the checkpoint.
    pub convention_schema: String,
    /// Checksum of the resolved configuration this boundary belongs to.
    pub config_checksum: String,
    /// Run identifier.
    pub run_id: String,
    /// Accepted step.
    pub accepted_step: u64,
    /// Parameter-layout fingerprint.
    pub parameter_layout_fingerprint: String,
    /// RNG state checksum.
    pub rng_state_checksum: String,
    /// Relative path of the serialized RNG state payload.
    pub rng_state_path: String,
    /// Payload checksum.
    pub payload_checksum: String,
    /// Relative path of the opaque checkpoint payload.
    pub payload_path: String,
    /// Typed accepted-boundary metadata.
    pub accepted_state: Option<AcceptedStateMetadata>,
    /// Typed RNG-state metadata.
    pub rng_state: Option<RngStateMetadata>,
    /// Named typed array payload metadata.
    pub arrays: Vec<CheckpointArray>,
}
/// Metadata for one checkpoint array payload.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CheckpointArray {
    /// Stable logical array name.
    pub name: String,
    /// Relative NPY payload path.
    pub path: String,
    /// Portable dtype identifier.
    pub dtype: String,
    /// Logical shape.
    pub shape: Vec<usize>,
    /// Memory order identifier.
    pub order: String,
    /// Checksum of the stored NPY payload.
    pub checksum: String,
}
impl CheckpointArray {
    /// Construct NPY array metadata.
    pub fn new(name: &str, dtype: &str, shape: Vec<usize>, order: &str, payload: &[u8]) -> Self {
        Self {
            name: name.into(),
            path: format!("{name}.npy"),
            dtype: dtype.into(),
            shape,
            order: order.into(),
            checksum: checksum(payload),
        }
    }
    /// Validate canonical shape and payload metadata.
    pub fn validate(&self, payload: &[u8]) -> Result<(), IoError> {
        let expected_values = self
            .shape
            .iter()
            .try_fold(1usize, |product, &size| product.checked_mul(size));
        if self.name.trim().is_empty()
            || self.path.is_empty()
            || self.path.starts_with('/')
            || self.path.split('/').any(|component| component == "..")
            || self.dtype != "f64"
            || self.order != "C"
            || self.shape.contains(&0)
            || expected_values.is_none()
            || payload.len() < 10
            || checksum(payload) != self.checksum
        {
            return Err(IoError::InvalidData(
                "invalid checkpoint array metadata or checksum".into(),
            ));
        }
        let (values, shape) = decode_npy(payload)?;
        if shape != self.shape || values.iter().any(|value| !value.is_finite()) {
            return Err(IoError::InvalidData(
                "checkpoint array shape or numeric values mismatch".into(),
            ));
        }
        Ok(())
    }
    /// Write a little-endian C-order NumPy `.npy` payload atomically.
    pub fn write_npy(&self, path: &Path, values: &[f64]) -> Result<(), IoError> {
        let payload = npy_bytes(&self.shape, values)?;
        atomic_write(path, &payload)
    }
    /// Construct metadata whose checksum matches encoded NPY values.
    pub fn from_values(
        name: &str,
        shape: Vec<usize>,
        values: &[f64],
    ) -> Result<(Self, Vec<u8>), IoError> {
        let payload = npy_bytes(&shape, values)?;
        Ok((
            Self {
                name: name.into(),
                path: format!("{name}.npy"),
                dtype: "f64".into(),
                shape,
                order: "C".into(),
                checksum: checksum(&payload),
            },
            payload,
        ))
    }
    /// Read and validate a little-endian C-order NumPy `.npy` payload.
    pub fn read_npy(path: &Path) -> Result<(Vec<f64>, Vec<usize>), IoError> {
        decode_npy(&fs::read(path)?)
    }
}

/// Materialized checkpoint data returned after validating every checksum and shape.
#[derive(Clone, Debug, PartialEq)]
pub struct CheckpointBundle {
    /// Checkpoint envelope.
    pub checkpoint: Checkpoint,
    /// Opaque accepted-state payload.
    pub payload: Vec<u8>,
    /// Serialized RNG state.
    pub rng_state: Vec<u8>,
    /// Named array metadata and decoded values.
    pub arrays: Vec<(CheckpointArray, Vec<f64>)>,
}

impl Checkpoint {
    /// Construct checkpoint metadata from a payload.
    pub fn new(
        run_id: &str,
        accepted_step: u64,
        layout: &str,
        rng_state_checksum: String,
        payload: &[u8],
    ) -> Self {
        Self {
            schema_version: CHECKPOINT_SCHEMA_VERSION.into(),
            convention_schema: String::new(),
            config_checksum: String::new(),
            run_id: run_id.into(),
            accepted_step,
            parameter_layout_fingerprint: layout.into(),
            rng_state_checksum,
            rng_state_path: "rng-state.bin".into(),
            payload_checksum: checksum(payload),
            payload_path: "payload.bin".into(),
            accepted_state: None,
            rng_state: None,
            arrays: Vec::new(),
        }
    }
    /// Attach the typed accepted-boundary and RNG state descriptions.
    pub fn with_state_metadata(
        mut self,
        accepted_state: AcceptedStateMetadata,
        rng_state: RngStateMetadata,
    ) -> Result<Self, IoError> {
        accepted_state.validate_semantics()?;
        rng_state.validate_semantics()?;
        if accepted_state.parameter_layout_fingerprint != self.parameter_layout_fingerprint
            || accepted_state.accepted_steps != self.accepted_step
        {
            return Err(IoError::InvalidData(
                "checkpoint state metadata disagrees with envelope".into(),
            ));
        }
        self.accepted_state = Some(accepted_state);
        self.rng_state = Some(rng_state);
        Ok(self)
    }
    /// Bind the checkpoint to a resolved configuration before publication.
    pub fn bind_config(mut self, config: &ScientificConfig) -> Result<Self, IoError> {
        self.convention_schema = config.conventions.schema.clone();
        self.config_checksum = config.checksum()?;
        Ok(self)
    }
    /// Validate the configuration and convention binding.
    pub fn validate_config(&self, config: &ScientificConfig) -> Result<(), IoError> {
        self.validate_semantics()?;
        if self.config_checksum != config.checksum()?
            || self.convention_schema != config.conventions.schema
        {
            return Err(IoError::InvalidData(
                "checkpoint/config binding mismatch".into(),
            ));
        }
        Ok(())
    }
    /// Attach named typed NPY array metadata.
    pub fn with_arrays(mut self, arrays: Vec<CheckpointArray>) -> Result<Self, IoError> {
        let mut names = BTreeSet::new();
        let mut paths = BTreeSet::new();
        if arrays.iter().any(|array| {
            let shape_product = array
                .shape
                .iter()
                .try_fold(1usize, |product, &size| product.checked_mul(size));
            array.name.trim().is_empty()
                || !names.insert(array.name.clone())
                || array.path.is_empty()
                || !paths.insert(array.path.clone())
                || array.path.starts_with('/')
                || array.path.split('/').any(|component| component == "..")
                || matches!(
                    array.path.as_str(),
                    "checkpoint.json" | "payload.bin" | "rng-state.bin"
                )
                || array.dtype != "f64"
                || array.order != "C"
                || array.shape.is_empty()
                || array.shape.contains(&0)
                || shape_product.is_none()
                || array.checksum.len() != 64
                || !is_hex(&array.checksum)
        }) {
            return Err(IoError::InvalidData(
                "invalid or duplicate checkpoint array metadata".into(),
            ));
        }
        self.arrays = arrays;
        Ok(self)
    }
    /// Validate checkpoint envelope fields without reading its payloads.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        let mut candidate = self.clone();
        let arrays = std::mem::take(&mut candidate.arrays);
        if self.schema_version != CHECKPOINT_SCHEMA_VERSION
            || self.convention_schema != CONVENTION_SCHEMA
            || self.config_checksum.len() != 64
            || !is_hex(&self.config_checksum)
            || self.run_id.trim().is_empty()
            || self.parameter_layout_fingerprint.trim().is_empty()
            || self.rng_state_path != "rng-state.bin"
            || self.payload_path != "payload.bin"
            || self.payload_checksum.len() != 64
            || !is_hex(&self.payload_checksum)
            || self.rng_state_checksum.len() != 64
            || !is_hex(&self.rng_state_checksum)
            || self.accepted_state.is_none()
            || self.rng_state.is_none()
        {
            return Err(IoError::InvalidData("invalid checkpoint metadata".into()));
        }
        candidate.arrays = arrays;
        let arrays = candidate.arrays.clone();
        candidate.with_arrays(arrays)?;
        let accepted_state = self.accepted_state.as_ref().ok_or_else(|| {
            IoError::InvalidData("checkpoint accepted-state metadata is missing".into())
        })?;
        accepted_state.validate_semantics()?;
        if accepted_state.parameter_layout_fingerprint != self.parameter_layout_fingerprint
            || accepted_state.accepted_steps != self.accepted_step
        {
            return Err(IoError::InvalidData(
                "accepted-state metadata disagrees with checkpoint envelope".into(),
            ));
        }
        let rng_state = self.rng_state.as_ref().ok_or_else(|| {
            IoError::InvalidData("checkpoint RNG-state metadata is missing".into())
        })?;
        rng_state.validate_semantics()?;
        if !self.arrays.iter().any(|array| {
            array.name == "parameters" && array.shape == accepted_state.parameter_shape
        }) {
            return Err(IoError::InvalidData(
                "checkpoint parameter array does not match accepted-state metadata".into(),
            ));
        }
        Ok(())
    }
    /// Validate an opaque checkpoint payload checksum.
    pub fn validate_payload(&self, payload: &[u8]) -> Result<(), IoError> {
        if checksum(payload) == self.payload_checksum {
            Ok(())
        } else {
            Err(IoError::InvalidData(
                "checkpoint payload checksum mismatch".into(),
            ))
        }
    }
    /// Validate a serialized RNG-state checksum.
    pub fn validate_rng_state(&self, rng_state: &[u8]) -> Result<(), IoError> {
        if checksum(rng_state) == self.rng_state_checksum {
            Ok(())
        } else {
            Err(IoError::InvalidData(
                "checkpoint RNG-state checksum mismatch".into(),
            ))
        }
    }
    /// Serialize checkpoint metadata as JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        strict_json(self)
    }
    /// Deserialize and validate checkpoint metadata.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self = strict_from_json(value)?;
        if parsed.schema_version != CHECKPOINT_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: CHECKPOINT_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
}

fn decode_npy(bytes: &[u8]) -> Result<(Vec<f64>, Vec<usize>), IoError> {
    if bytes.len() < 12 || &bytes[..6] != b"\x93NUMPY" || bytes[6..8] != [1, 0] {
        return Err(IoError::InvalidData("unsupported NPY header".into()));
    }
    let header_length = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
    let header_end = 10usize
        .checked_add(header_length)
        .ok_or_else(|| IoError::InvalidData("NPY header overflow".into()))?;
    if header_end > bytes.len() {
        return Err(IoError::InvalidData("truncated NPY header".into()));
    }
    let header = std::str::from_utf8(&bytes[10..header_end])
        .map_err(|_| IoError::InvalidData("NPY header is not UTF-8".into()))?;
    if !header.contains("'descr': '<f8'") || !header.contains("'fortran_order': False") {
        return Err(IoError::InvalidData(
            "unsupported NPY dtype or order".into(),
        ));
    }
    let shape_start = header
        .find("'shape': (")
        .ok_or_else(|| IoError::InvalidData("NPY shape is missing".into()))?
        + 10;
    let shape_end = header[shape_start..]
        .find(')')
        .ok_or_else(|| IoError::InvalidData("NPY shape is malformed".into()))?
        + shape_start;
    let shape = header[shape_start..shape_end]
        .split(',')
        .filter(|value| !value.trim().is_empty())
        .map(|value| {
            value
                .trim()
                .parse::<usize>()
                .map_err(|_| IoError::InvalidData("NPY shape contains a non-integer".into()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if shape.is_empty() || shape.contains(&0) {
        return Err(IoError::InvalidData(
            "NPY shape is empty or zero-sized".into(),
        ));
    }
    let expected = shape
        .iter()
        .try_fold(1usize, |product, &size| product.checked_mul(size))
        .ok_or_else(|| IoError::InvalidData("NPY shape overflow".into()))?;
    let data = &bytes[header_end..];
    if data.len()
        != expected
            .checked_mul(8)
            .ok_or_else(|| IoError::InvalidData("NPY payload overflow".into()))?
    {
        return Err(IoError::InvalidData("NPY payload length mismatch".into()));
    }
    let mut values = Vec::with_capacity(expected);
    for chunk in data.chunks_exact(8) {
        let bytes: [u8; 8] = chunk
            .try_into()
            .map_err(|_| IoError::InvalidData("NPY scalar width mismatch".into()))?;
        values.push(f64::from_le_bytes(bytes));
    }
    Ok((values, shape))
}

/// One accepted trajectory row.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrajectoryRow {
    /// Accepted step.
    pub step: u64,
    /// Physical time.
    pub time: f64,
    /// Observable energy.
    pub energy: f64,
}
impl TrajectoryRow {
    /// Construct a trajectory row.
    pub fn new(step: u64, time: f64, energy: f64) -> Self {
        Self { step, time, energy }
    }
}

/// JSON columnar trajectory representation selected for the initial IO API.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ColumnarTrajectory {
    /// Schema version.
    pub schema_version: String,
    /// Scientific convention schema stored in Parquet metadata.
    pub convention_schema: String,
    /// Checksum of the resolved configuration stored in Parquet metadata.
    pub config_checksum: String,
    /// Accepted-step column.
    pub steps: Vec<u64>,
    /// Time column.
    pub times: Vec<f64>,
    /// Energy column.
    pub energies: Vec<f64>,
}
impl ColumnarTrajectory {
    /// Build columns from rows after checking finite values.
    pub fn new(rows: Vec<TrajectoryRow>) -> Result<Self, IoError> {
        Self::from_columns(
            rows.iter().map(|r| r.step).collect(),
            rows.iter().map(|r| r.time).collect(),
            rows.iter().map(|r| r.energy).collect(),
        )
    }
    /// Build columns and validate equal lengths and finite values.
    pub fn from_columns(
        steps: Vec<u64>,
        times: Vec<f64>,
        energies: Vec<f64>,
    ) -> Result<Self, IoError> {
        if steps.len() != times.len()
            || steps.len() != energies.len()
            || times.iter().chain(energies.iter()).any(|v| !v.is_finite())
        {
            return Err(IoError::InvalidData(
                "trajectory columns have unequal length or non-finite value".into(),
            ));
        }
        Ok(Self {
            schema_version: TRAJECTORY_SCHEMA_VERSION.into(),
            convention_schema: String::new(),
            config_checksum: String::new(),
            steps,
            times,
            energies,
        })
    }
    /// Number of rows.
    pub fn len(&self) -> usize {
        self.steps.len()
    }
    /// Whether the table contains no rows.
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }
    /// Bind the trajectory to a resolved configuration before publication.
    pub fn bind_config(mut self, config: &ScientificConfig) -> Result<Self, IoError> {
        self.convention_schema = config.conventions.schema.clone();
        self.config_checksum = config.checksum()?;
        self.validate_semantics()?;
        Ok(self)
    }
    /// Validate columns and their scientific configuration binding.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        if self.schema_version != TRAJECTORY_SCHEMA_VERSION
            || self.convention_schema != CONVENTION_SCHEMA
            || self.config_checksum.len() != 64
            || !is_hex(&self.config_checksum)
            || self.steps.len() != self.times.len()
            || self.steps.len() != self.energies.len()
            || self
                .times
                .iter()
                .chain(self.energies.iter())
                .any(|value| !value.is_finite())
        {
            return Err(IoError::InvalidData("invalid trajectory semantics".into()));
        }
        Ok(())
    }
    /// Serialize the columnar table as JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        serde_json::to_string_pretty(self).map_err(|e| IoError::Serialization(e.to_string()))
    }
    /// Deserialize and validate the columnar table.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self =
            serde_json::from_str(value).map_err(|e| IoError::Serialization(e.to_string()))?;
        if parsed.schema_version != TRAJECTORY_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: TRAJECTORY_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        let mut trajectory = Self::from_columns(parsed.steps, parsed.times, parsed.energies)?;
        trajectory.convention_schema = parsed.convention_schema;
        trajectory.config_checksum = parsed.config_checksum;
        trajectory.validate_semantics()?;
        Ok(trajectory)
    }
    /// Write the trajectory as an immutable Apache Parquet part file.
    pub fn write_parquet(&self, path: &Path) -> Result<(), IoError> {
        use arrow_array::{ArrayRef, Float64Array, RecordBatch, UInt64Array};
        use arrow_schema::{DataType, Field, Schema};
        use parquet::arrow::ArrowWriter;
        self.validate_semantics()?;
        let mut metadata = HashMap::new();
        metadata.insert(
            "qslib_convention_schema".into(),
            self.convention_schema.clone(),
        );
        metadata.insert("qslib_config_checksum".into(), self.config_checksum.clone());
        let schema = Arc::new(
            Schema::new(vec![
                Field::new("step", DataType::UInt64, false),
                Field::new("time", DataType::Float64, false),
                Field::new("energy", DataType::Float64, false),
            ])
            .with_metadata(metadata),
        );
        let columns: Vec<ArrayRef> = vec![
            Arc::new(UInt64Array::from(self.steps.clone())),
            Arc::new(Float64Array::from(self.times.clone())),
            Arc::new(Float64Array::from(self.energies.clone())),
        ];
        let batch = RecordBatch::try_new(schema.clone(), columns)
            .map_err(|error| IoError::InvalidData(error.to_string()))?;
        let parent = path
            .parent()
            .ok_or_else(|| IoError::InvalidData("parquet target has no parent".into()))?;
        if path.exists() {
            return Err(IoError::InvalidData(
                "Parquet part already exists and is immutable".into(),
            ));
        }
        fs::create_dir_all(parent)?;
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|_| IoError::InvalidData("system clock before epoch".into()))?
            .as_nanos();
        let temporary = parent.join(format!(
            ".{}.part-{nonce}",
            path.file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("trajectory")
        ));
        let result = (|| {
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)?;
            let mut writer = ArrowWriter::try_new(file, schema, None)
                .map_err(|error| IoError::InvalidData(error.to_string()))?;
            writer
                .write(&batch)
                .map_err(|error| IoError::InvalidData(error.to_string()))?;
            writer
                .close()
                .map_err(|error| IoError::InvalidData(error.to_string()))?;
            fs::File::open(&temporary)?.sync_all()?;
            fs::rename(&temporary, path)?;
            sync_directory(parent)?;
            Ok::<(), IoError>(())
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temporary);
        }
        result
    }
    /// Read a Parquet part through Arrow and reconstruct validated columns.
    pub fn read_parquet(path: &Path) -> Result<Self, IoError> {
        use arrow_array::{Array, Float64Array, UInt64Array};
        use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
        let file = fs::File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|error| IoError::InvalidData(error.to_string()))?;
        let metadata = builder.schema().metadata().clone();
        let reader = builder
            .build()
            .map_err(|error| IoError::InvalidData(error.to_string()))?;
        let mut steps = Vec::new();
        let mut times = Vec::new();
        let mut energies = Vec::new();
        for batch in reader {
            let batch = batch.map_err(|error| IoError::InvalidData(error.to_string()))?;
            if batch.num_columns() != 3
                || batch.schema().field(0).name() != "step"
                || batch.schema().field(1).name() != "time"
                || batch.schema().field(2).name() != "energy"
            {
                return Err(IoError::InvalidData(
                    "unexpected Parquet trajectory schema".into(),
                ));
            }
            let step = batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| IoError::InvalidData("step column dtype mismatch".into()))?;
            let time = batch
                .column(1)
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| IoError::InvalidData("time column dtype mismatch".into()))?;
            let energy = batch
                .column(2)
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| IoError::InvalidData("energy column dtype mismatch".into()))?;
            steps.extend(step.values());
            times.extend(time.values());
            energies.extend(energy.values());
        }
        let mut trajectory = Self::from_columns(steps, times, energies)?;
        trajectory.convention_schema = metadata
            .get("qslib_convention_schema")
            .cloned()
            .unwrap_or_default();
        trajectory.config_checksum = metadata
            .get("qslib_config_checksum")
            .cloned()
            .unwrap_or_default();
        trajectory.validate_semantics()?;
        Ok(trajectory)
    }
}

/// Atomic manifest for an append-only set of immutable Parquet trajectory parts.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParquetDatasetManifest {
    /// Dataset manifest schema version.
    pub schema_version: String,
    /// Scientific convention schema used by every part.
    pub convention_schema: String,
    /// Checksum of the resolved scientific configuration.
    pub config_checksum: String,
    /// Immutable Parquet parts in append order.
    pub parts: Vec<ArtifactEntry>,
    /// Whether the completion marker has been published.
    pub complete: bool,
}

impl ParquetDatasetManifest {
    /// Create an incomplete dataset manifest bound to a configuration checksum.
    pub fn new(config_checksum: String) -> Result<Self, IoError> {
        let value = Self {
            schema_version: PARQUET_DATASET_SCHEMA_VERSION.into(),
            convention_schema: CONVENTION_SCHEMA.into(),
            config_checksum,
            parts: Vec::new(),
            complete: false,
        };
        value.validate_semantics()?;
        Ok(value)
    }
    /// Append one immutable Parquet part and record its exact bytes.
    pub fn append_part(
        &mut self,
        directory: &Path,
        filename: &str,
        trajectory: &ColumnarTrajectory,
    ) -> Result<ArtifactEntry, IoError> {
        if self.complete || !safe_relative_path(filename) || !filename.ends_with(".parquet") {
            return Err(IoError::InvalidData(
                "invalid or completed Parquet dataset".into(),
            ));
        }
        if self.parts.iter().any(|part| part.path == filename) {
            return Err(IoError::InvalidData("duplicate Parquet part".into()));
        }
        let path = directory.join(filename);
        trajectory.write_parquet(&path)?;
        let bytes = fs::read(&path)?;
        let entry = ArtifactEntry::new(filename, checksum(&bytes), bytes.len() as u64);
        self.parts.push(entry.clone());
        self.validate_semantics()?;
        Ok(entry)
    }
    /// Write the current incomplete or complete manifest atomically.
    pub fn write_manifest(&self, directory: &Path) -> Result<(), IoError> {
        self.validate_semantics()?;
        atomic_write(&directory.join("manifest.json"), self.to_json()?.as_bytes())
    }
    /// Publish the complete manifest and then its completion marker.
    pub fn finish(&mut self, directory: &Path) -> Result<(), IoError> {
        if self.complete {
            return Err(IoError::InvalidData(
                "Parquet dataset is already complete".into(),
            ));
        }
        // Validate all part files before changing the in-memory state. This
        // permits a retry if marker publication fails.
        self.validate_parts(directory)?;
        let mut candidate = self.clone();
        candidate.complete = true;
        candidate.validate_semantics()?;
        candidate.write_manifest(directory)?;
        let marker = atomic_write(
            &directory.join(DATASET_COMPLETE_MARKER),
            b"qslib-dataset-complete-v1\n",
        );
        if marker.is_ok() {
            self.complete = true;
        }
        marker
    }
    /// Serialize the manifest as strict JSON.
    pub fn to_json(&self) -> Result<String, IoError> {
        self.validate_semantics()?;
        strict_json(self)
    }
    /// Deserialize a manifest without assuming the dataset is present.
    pub fn from_json(value: &str) -> Result<Self, IoError> {
        let parsed: Self = strict_from_json(value)?;
        if parsed.schema_version != PARQUET_DATASET_SCHEMA_VERSION {
            return Err(IoError::UnsupportedVersion {
                expected: PARQUET_DATASET_SCHEMA_VERSION,
                actual: parsed.schema_version,
            });
        }
        parsed.validate_semantics()?;
        Ok(parsed)
    }
    /// Load a complete dataset and validate marker, paths, sizes, and checksums.
    pub fn load(directory: &Path) -> Result<Self, IoError> {
        let value = Self::from_json(&fs::read_to_string(directory.join("manifest.json"))?)?;
        if !value.complete {
            return Err(IoError::InvalidData(
                "Parquet dataset has no completed marker".into(),
            ));
        }
        if fs::read(directory.join(DATASET_COMPLETE_MARKER))
            .ok()
            .as_deref()
            == Some(b"qslib-dataset-complete-v1\n")
        {
            value.validate_parts(directory)?;
            Ok(value)
        } else {
            Self::recover(directory)
        }
    }
    /// Inspect a complete dataset without repairing or otherwise mutating it.
    ///
    /// Unlike [`Self::load`], this method reports a missing or incorrect
    /// completion marker as an error and never republishes one. It is intended
    /// for read-only CLI and audit operations.
    pub fn inspect(directory: &Path) -> Result<Self, IoError> {
        let value = Self::from_json(&fs::read_to_string(directory.join("manifest.json"))?)?;
        if !value.complete {
            return Err(IoError::InvalidData(
                "Parquet dataset is not marked complete".into(),
            ));
        }
        if fs::read(directory.join(DATASET_COMPLETE_MARKER))?.as_slice()
            != b"qslib-dataset-complete-v1\n"
        {
            return Err(IoError::InvalidData(
                "Parquet dataset completion marker is missing or invalid".into(),
            ));
        }
        value.validate_parts(directory)?;
        Ok(value)
    }
    /// Recover a dataset whose complete manifest was durable but whose marker
    /// was interrupted or lost. Recovery revalidates every immutable part,
    /// then atomically republishes the exact marker before returning.
    pub fn recover(directory: &Path) -> Result<Self, IoError> {
        let value = Self::from_json(&fs::read_to_string(directory.join("manifest.json"))?)?;
        if !value.complete {
            return Err(IoError::InvalidData(
                "cannot recover an incomplete Parquet dataset".into(),
            ));
        }
        value.validate_parts(directory)?;
        atomic_write(
            &directory.join(DATASET_COMPLETE_MARKER),
            b"qslib-dataset-complete-v1\n",
        )?;
        Ok(value)
    }
    /// Validate that this manifest is bound to a resolved configuration.
    pub fn validate_config(&self, config: &ScientificConfig) -> Result<(), IoError> {
        if self.convention_schema == config.conventions.schema
            && self.config_checksum == config.checksum()?
        {
            Ok(())
        } else {
            Err(IoError::InvalidData(
                "Parquet/config binding mismatch".into(),
            ))
        }
    }
    fn validate_parts(&self, directory: &Path) -> Result<(), IoError> {
        for part in &self.parts {
            let bytes = fs::read(directory.join(&part.path))?;
            if part.bytes != bytes.len() as u64 || part.checksum != checksum(&bytes) {
                return Err(IoError::InvalidData(
                    "Parquet part checksum or length mismatch".into(),
                ));
            }
            ColumnarTrajectory::read_parquet(&directory.join(&part.path))?;
        }
        Ok(())
    }
    /// Validate schema, path safety, checksums, and duplicate parts.
    pub fn validate_semantics(&self) -> Result<(), IoError> {
        let mut paths = BTreeSet::new();
        if self.schema_version != PARQUET_DATASET_SCHEMA_VERSION
            || self.convention_schema != CONVENTION_SCHEMA
            || self.config_checksum.len() != 64
            || !is_hex(&self.config_checksum)
            || self.parts.iter().any(|part| {
                !safe_relative_path(&part.path)
                    || !part.path.ends_with(".parquet")
                    || part.bytes == 0
                    || part.checksum.len() != 64
                    || !is_hex(&part.checksum)
                    || !paths.insert(part.path.clone())
            })
        {
            return Err(IoError::InvalidData(
                "invalid Parquet dataset manifest".into(),
            ));
        }
        Ok(())
    }
}

fn npy_bytes(shape: &[usize], values: &[f64]) -> Result<Vec<u8>, IoError> {
    let expected = shape
        .iter()
        .try_fold(1usize, |product, &size| product.checked_mul(size))
        .ok_or_else(|| IoError::InvalidData("checkpoint array shape overflow".into()))?;
    if expected != values.len()
        || shape.is_empty()
        || shape.contains(&0)
        || values.iter().any(|value| !value.is_finite())
    {
        return Err(IoError::InvalidData(
            "NPY payload does not match checkpoint array metadata".into(),
        ));
    }
    let shape_text = if shape.len() == 1 {
        format!("({},)", shape[0])
    } else {
        format!(
            "({})",
            shape
                .iter()
                .map(|size| size.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    let mut header =
        format!("{{'descr': '<f8', 'fortran_order': False, 'shape': {shape_text}, }}").into_bytes();
    let padding = (16 - ((10 + header.len() + 1) % 16)) % 16;
    header.extend(std::iter::repeat_n(b' ', padding));
    header.push(b'\n');
    if header.len() > u16::MAX as usize {
        return Err(IoError::InvalidData(
            "NPY header exceeds version 1 limit".into(),
        ));
    }
    let mut payload = b"\x93NUMPY\x01\x00".to_vec();
    payload.extend((header.len() as u16).to_le_bytes());
    payload.extend(header);
    for value in values {
        payload.extend(value.to_le_bytes());
    }
    Ok(payload)
}

fn strict_json<T: Serialize>(value: &T) -> Result<String, IoError> {
    serde_json::to_string_pretty(value).map_err(|error| IoError::Serialization(error.to_string()))
}

fn strict_from_json<T: DeserializeOwned>(value: &str) -> Result<T, IoError> {
    serde_json::from_str(value).map_err(|error| IoError::Serialization(error.to_string()))
}

fn is_hex(value: &str) -> bool {
    value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn is_canonical_channel(value: &str) -> bool {
    matches!(
        value,
        "ising_zz" | "heisenberg_exchange" | "rydberg_density_density"
    ) || value
        .strip_prefix("generic:")
        .is_some_and(|generic| !generic.is_empty() && !generic.chars().any(char::is_whitespace))
}

fn canonical_interaction_identity(interaction: &InteractionMetadata) -> String {
    let name = interaction.name.as_deref().unwrap_or("");
    format!(
        "bond:{},{};image:{},{};source:{};direction:{},{};channel:{}:{};name:{}:{}",
        interaction.first,
        interaction.second,
        interaction.image_x,
        interaction.image_y,
        interaction.source,
        interaction.direction_x,
        interaction.direction_y,
        interaction.channel.len(),
        interaction.channel,
        name.len(),
        name
    )
}

fn safe_relative_path(value: &str) -> bool {
    !value.is_empty()
        && !value.starts_with('/')
        && !value.contains('\\')
        && value
            .split('/')
            .all(|component| !component.is_empty() && component != "." && component != "..")
}

/// Compute a lowercase BLAKE3 checksum.
pub fn checksum(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

/// Atomically write bytes by completing a sibling temporary file before rename.
pub fn atomic_write(path: &Path, bytes: &[u8]) -> Result<(), IoError> {
    let parent = path
        .parent()
        .ok_or_else(|| IoError::InvalidData("target has no parent".into()))?;
    fs::create_dir_all(parent)?;
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|_| IoError::InvalidData("system clock before epoch".into()))?
        .as_nanos();
    let temporary = parent.join(format!(
        ".{}.tmp-{}-{nonce}",
        path.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("artifact"),
        std::process::id()
    ));
    let result = (|| {
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        file.write_all(bytes)?;
        file.sync_all()?;
        fs::rename(&temporary, path)?;
        sync_directory(parent)?;
        Ok::<(), IoError>(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

/// Flush a directory entry after an atomic rename where the platform exposes
/// a portable directory handle. Windows does not permit opening a directory
/// as a synchronizable `File`, while the rename itself still provides the
/// supported replacement boundary there.
fn sync_directory(path: &Path) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        fs::File::open(path)?.sync_all()
    }
    #[cfg(not(unix))]
    {
        let _ = path;
        Ok(())
    }
}

/// Atomically publish a checkpoint envelope and its named NPY arrays.
pub fn write_checkpoint_bundle(
    directory: &Path,
    checkpoint: &Checkpoint,
    payload: &[u8],
    rng_state: &[u8],
    arrays: &[(&CheckpointArray, &[f64])],
) -> Result<(), IoError> {
    checkpoint.validate_semantics()?;
    checkpoint.validate_payload(payload)?;
    checkpoint.validate_rng_state(rng_state)?;
    let accepted_state = AcceptedStateMetadata::from_bytes(payload)?;
    let declared_state = checkpoint.accepted_state.as_ref().ok_or_else(|| {
        IoError::InvalidData("checkpoint accepted-state metadata is missing".into())
    })?;
    if &accepted_state != declared_state {
        return Err(IoError::InvalidData(
            "accepted-state payload disagrees with checkpoint metadata".into(),
        ));
    }
    let decoded_rng = RngStateMetadata::from_bytes(rng_state)?;
    let declared_rng = checkpoint
        .rng_state
        .as_ref()
        .ok_or_else(|| IoError::InvalidData("checkpoint RNG-state metadata is missing".into()))?;
    if &decoded_rng != declared_rng {
        return Err(IoError::InvalidData(
            "RNG-state payload disagrees with checkpoint metadata".into(),
        ));
    }
    if arrays.len() != checkpoint.arrays.len() {
        return Err(IoError::InvalidData(
            "checkpoint array count mismatch".into(),
        ));
    }
    let mut supplied_names = BTreeSet::new();
    let mut supplied_paths = BTreeSet::new();
    if arrays.iter().any(|(metadata, _)| {
        !supplied_names.insert(metadata.name.clone())
            || !supplied_paths.insert(metadata.path.clone())
            || !checkpoint
                .arrays
                .iter()
                .any(|declared| declared == *metadata)
    }) {
        return Err(IoError::InvalidData(
            "checkpoint supplied arrays are duplicated or undeclared".into(),
        ));
    }
    if checkpoint
        .arrays
        .iter()
        .any(|declared| !supplied_names.contains(&declared.name))
    {
        return Err(IoError::InvalidData(
            "checkpoint supplied array set is incomplete".into(),
        ));
    }
    let parent = directory
        .parent()
        .ok_or_else(|| IoError::InvalidData("checkpoint directory has no parent".into()))?;
    fs::create_dir_all(parent)?;
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|_| IoError::InvalidData("system clock before epoch".into()))?
        .as_nanos();
    let temporary = parent.join(format!(".checkpoint-{nonce}.tmp-{}", std::process::id()));
    if directory.exists() {
        return Err(IoError::InvalidData(
            "checkpoint target already exists; choose a new accepted boundary".into(),
        ));
    }
    let result = (|| {
        fs::create_dir(&temporary)?;
        atomic_write(&temporary.join(&checkpoint.payload_path), payload)?;
        atomic_write(&temporary.join(&checkpoint.rng_state_path), rng_state)?;
        for (metadata, values) in arrays {
            let expected = checkpoint
                .arrays
                .iter()
                .find(|candidate| candidate.name == metadata.name)
                .ok_or_else(|| IoError::InvalidData("checkpoint array is not declared".into()))?;
            if expected != *metadata {
                return Err(IoError::InvalidData(
                    "checkpoint array metadata does not match envelope".into(),
                ));
            }
            let encoded = npy_bytes(&metadata.shape, values)?;
            metadata.validate(&encoded)?;
            atomic_write(&temporary.join(&metadata.path), &encoded)?;
        }
        atomic_write(
            &temporary.join("checkpoint.json"),
            checkpoint.to_json()?.as_bytes(),
        )?;
        sync_directory(&temporary)?;
        fs::rename(&temporary, directory)?;
        sync_directory(parent)?;
        Ok::<(), IoError>(())
    })();
    if result.is_err() {
        let _ = fs::remove_dir_all(&temporary);
    }
    result
}

/// Read and validate a complete checkpoint bundle.
pub fn read_checkpoint_bundle(directory: &Path) -> Result<CheckpointBundle, IoError> {
    let checkpoint =
        Checkpoint::from_json(&fs::read_to_string(directory.join("checkpoint.json"))?)?;
    let payload = fs::read(directory.join(&checkpoint.payload_path))?;
    checkpoint.validate_payload(&payload)?;
    if AcceptedStateMetadata::from_bytes(&payload)?
        != checkpoint.accepted_state.clone().ok_or_else(|| {
            IoError::InvalidData("checkpoint accepted-state metadata is missing".into())
        })?
    {
        return Err(IoError::InvalidData(
            "accepted-state payload disagrees with checkpoint metadata".into(),
        ));
    }
    let rng_state = fs::read(directory.join(&checkpoint.rng_state_path))?;
    checkpoint.validate_rng_state(&rng_state)?;
    if RngStateMetadata::from_bytes(&rng_state)?
        != checkpoint.rng_state.clone().ok_or_else(|| {
            IoError::InvalidData("checkpoint RNG-state metadata is missing".into())
        })?
    {
        return Err(IoError::InvalidData(
            "RNG-state payload disagrees with checkpoint metadata".into(),
        ));
    }
    let mut arrays = Vec::with_capacity(checkpoint.arrays.len());
    for metadata in &checkpoint.arrays {
        let encoded = fs::read(directory.join(&metadata.path))?;
        metadata.validate(&encoded)?;
        let (values, shape) = decode_npy(&encoded)?;
        if shape != metadata.shape {
            return Err(IoError::InvalidData(
                "checkpoint array shape differs from envelope".into(),
            ));
        }
        arrays.push((metadata.clone(), values));
    }
    Ok(CheckpointBundle {
        checkpoint,
        payload,
        rng_state,
        arrays,
    })
}
