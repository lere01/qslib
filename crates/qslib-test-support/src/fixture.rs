use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt::{self, Display, Formatter};

use serde::Deserialize;
use serde::de::DeserializeOwned;
use serde_json::Value;

/// Schema identifier for every neutral qslib conformance fixture.
pub const FIXTURE_SCHEMA: &str = "qslib-conformance-fixture-v1";

/// Scientific convention identifier required by the first fixture schema.
pub const CONVENTION_SCHEMA: &str = "qslib-conventions-v1";

const MANIFEST_SCHEMA: &str = "qslib-conformance-manifest-v1";

/// A category of independently derived conformance evidence.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd)]
#[serde(rename_all = "snake_case")]
pub enum FixtureKind {
    /// Rectangular site indexing and nearest-neighbour bonds.
    RectangularIndexing,
    /// Little-endian bit packing and local physical values.
    BitPacking,
    /// A two-site transverse-field Ising matrix.
    TfimOneBond,
    /// A two-site isotropic Heisenberg matrix.
    HeisenbergOneBond,
    /// State action with pair-dependent Heisenberg couplings.
    HeisenbergHeterogeneous,
    /// A two-site driven Rydberg matrix.
    RydbergTwoSite,
    /// Product-state observable normalizations.
    ObservableNormalization,
    /// Spectrum parity under a Hadamard basis rotation.
    BasisRotationSpectrum,
}

/// Provenance for an expected result that is independent of qslib production code.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Oracle {
    /// Independent calculation method, such as explicit Pauli matrices.
    pub method: String,
    /// Stable locator for the governing scientific statement.
    pub source: String,
    /// Revision of the source used for this derivation.
    pub source_revision: String,
    /// Short reproducible derivation of the expected result.
    pub derivation: String,
    /// Whether the expected value was produced without the implementation under test.
    pub independent_of_production: bool,
}

/// Authorship and review provenance for a neutral fixture.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Authorship {
    /// Role or group that prepared the fixture.
    pub prepared_by: String,
    /// Independent reviewer role or identity.
    pub reviewed_by: String,
    /// Review date in ISO `YYYY-MM-DD` form.
    pub reviewed_on: String,
}

/// One validated, language-neutral scientific conformance fixture.
#[derive(Clone, Debug)]
pub struct Fixture {
    /// Fixture envelope schema identifier.
    pub fixture_schema: String,
    /// Normative scientific convention schema identifier.
    pub convention_schema: String,
    /// Stable identifier unique within the fixture set.
    pub id: String,
    /// Case-specific fixture category.
    pub kind: FixtureKind,
    /// Physical or representational claim made by the fixture.
    pub claim: String,
    /// Independent source and derivation of expected values.
    pub oracle: Oracle,
    /// Preparation and independent-review provenance.
    pub authorship: Authorship,
    /// Explicit convention fields needed to interpret the payload.
    pub conventions: Value,
    /// Case-specific expected values and physical inputs.
    pub data: Value,
    /// Quantity-specific comparison policy for expected values.
    pub comparison: ComparisonPolicy,
}

/// A structured failure while parsing or validating neutral fixture evidence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FixtureError {
    message: String,
}

impl FixtureError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl Display for FixtureError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl Error for FixtureError {}

impl Fixture {
    /// Parse and validate one fixture from UTF-8 JSON.
    pub fn from_json(json: &str) -> Result<Self, FixtureError> {
        let value: Value = serde_json::from_str(json)
            .map_err(|error| FixtureError::new(format!("invalid fixture JSON: {error}")))?;
        validate_top_level_keys(&value)?;
        let raw: RawFixture = serde_json::from_value(value)
            .map_err(|error| FixtureError::new(format!("invalid fixture envelope: {error}")))?;
        let fixture = Self {
            fixture_schema: raw.fixture_schema,
            convention_schema: raw.convention_schema,
            id: raw.id,
            kind: raw.kind,
            claim: raw.claim,
            oracle: raw.oracle,
            authorship: raw.authorship,
            conventions: raw.conventions,
            data: raw.data,
            comparison: raw.comparison,
        };
        fixture.validate()?;
        Ok(fixture)
    }

    /// Validate schema, provenance, comparison policy, and case-specific structure.
    pub fn validate(&self) -> Result<(), FixtureError> {
        ensure(
            self.fixture_schema == FIXTURE_SCHEMA,
            "unsupported fixture schema",
        )?;
        ensure(
            self.convention_schema == CONVENTION_SCHEMA,
            "unsupported convention schema",
        )?;
        ensure_nonempty(&self.id, "fixture id")?;
        ensure_nonempty(&self.claim, "fixture claim")?;
        ensure_nonempty(&self.oracle.method, "oracle method")?;
        ensure_nonempty(&self.oracle.source, "oracle source")?;
        ensure_nonempty(&self.oracle.source_revision, "oracle source revision")?;
        ensure_nonempty(&self.oracle.derivation, "oracle derivation")?;
        ensure(
            self.oracle.independent_of_production,
            "oracle must be independent of production code",
        )?;
        ensure_nonempty(&self.authorship.prepared_by, "fixture preparer")?;
        ensure_nonempty(&self.authorship.reviewed_by, "fixture reviewer")?;
        validate_date(&self.authorship.reviewed_on)?;
        ensure_nonempty_object(&self.conventions, "resolved conventions")?;
        ensure_nonempty_object(&self.data, "fixture data")?;
        self.comparison.validate()?;
        validate_case(self.kind, &self.data, self.comparison.tolerance())
            .map_err(|error| FixtureError::new(format!("fixture {}: {error}", self.id)))
    }

    /// Return the quantity-specific comparison policy.
    pub const fn comparison(&self) -> &ComparisonPolicy {
        &self.comparison
    }
}

/// Return the complete set of fixture kinds required at the Milestone 1 gate.
pub fn required_fixture_kinds() -> Vec<FixtureKind> {
    vec![
        FixtureKind::RectangularIndexing,
        FixtureKind::BitPacking,
        FixtureKind::TfimOneBond,
        FixtureKind::HeisenbergOneBond,
        FixtureKind::HeisenbergHeterogeneous,
        FixtureKind::RydbergTwoSite,
        FixtureKind::ObservableNormalization,
        FixtureKind::BasisRotationSpectrum,
    ]
}

/// Validate an explicit manifest and its named fixture byte sequences.
///
/// This lower-level entry point supports negative harness tests and other
/// language adapters. Paths are logical manifest paths, not filesystem paths.
pub fn validate_fixture_set(
    manifest_json: &str,
    files: &[(&str, &[u8])],
) -> Result<Vec<Fixture>, FixtureError> {
    let manifest: FixtureManifest = serde_json::from_str(manifest_json)
        .map_err(|error| FixtureError::new(format!("invalid fixture manifest: {error}")))?;
    manifest.validate()?;

    let mut available = BTreeMap::new();
    for (path, bytes) in files {
        ensure(
            available.insert((*path).to_owned(), *bytes).is_none(),
            format!("duplicate fixture path: {path}"),
        )?;
    }
    ensure(
        available.len() == manifest.entries.len(),
        "manifest and supplied fixture counts differ",
    )?;

    let mut fixtures = Vec::with_capacity(manifest.entries.len());
    let mut ids = BTreeSet::new();
    let mut kinds = BTreeSet::new();
    for entry in &manifest.entries {
        let bytes = available
            .remove(&entry.path)
            .ok_or_else(|| FixtureError::new(format!("missing fixture file: {}", entry.path)))?;
        let actual_digest = blake3::hash(bytes).to_hex().to_string();
        ensure(
            actual_digest == entry.blake3,
            format!("BLAKE3 mismatch for {}", entry.path),
        )?;
        let json = std::str::from_utf8(bytes)
            .map_err(|error| FixtureError::new(format!("{} is not UTF-8: {error}", entry.path)))?;
        let fixture = Fixture::from_json(json)?;
        ensure(fixture.id == entry.id, "manifest fixture id mismatch")?;
        ensure(fixture.kind == entry.kind, "manifest fixture kind mismatch")?;
        ensure(ids.insert(fixture.id.clone()), "duplicate fixture id")?;
        ensure(kinds.insert(fixture.kind), "duplicate fixture kind")?;
        fixtures.push(fixture);
    }
    ensure(available.is_empty(), "unlisted fixture file supplied")?;
    let required: BTreeSet<_> = required_fixture_kinds().into_iter().collect();
    ensure(
        kinds == required,
        "fixture set is incomplete or contains an extra kind",
    )?;
    Ok(fixtures)
}

/// Load and validate the complete fixture set embedded in this crate.
pub fn load_conformance_fixtures() -> Result<Vec<Fixture>, FixtureError> {
    const MANIFEST: &str = include_str!("../../../fixtures/conformance/v1/manifest.json");
    const FILES: &[(&str, &[u8])] = &[
        (
            "basis_rotation_spectrum.json",
            include_bytes!("../../../fixtures/conformance/v1/basis_rotation_spectrum.json"),
        ),
        (
            "bit_packing.json",
            include_bytes!("../../../fixtures/conformance/v1/bit_packing.json"),
        ),
        (
            "heisenberg_heterogeneous.json",
            include_bytes!("../../../fixtures/conformance/v1/heisenberg_heterogeneous.json"),
        ),
        (
            "heisenberg_one_bond.json",
            include_bytes!("../../../fixtures/conformance/v1/heisenberg_one_bond.json"),
        ),
        (
            "observable_normalization.json",
            include_bytes!("../../../fixtures/conformance/v1/observable_normalization.json"),
        ),
        (
            "rectangular_indexing.json",
            include_bytes!("../../../fixtures/conformance/v1/rectangular_indexing.json"),
        ),
        (
            "rydberg_two_site.json",
            include_bytes!("../../../fixtures/conformance/v1/rydberg_two_site.json"),
        ),
        (
            "tfim_one_bond.json",
            include_bytes!("../../../fixtures/conformance/v1/tfim_one_bond.json"),
        ),
    ];
    validate_fixture_set(MANIFEST, FILES)
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawFixture {
    fixture_schema: String,
    convention_schema: String,
    id: String,
    kind: FixtureKind,
    claim: String,
    oracle: Oracle,
    authorship: Authorship,
    conventions: Value,
    comparison: ComparisonPolicy,
    data: Value,
}

/// Comparison policy recorded by a conformance fixture.
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ComparisonPolicy {
    /// Require exact equality for the recorded quantity.
    Exact {
        /// Quantity being compared.
        quantity: String,
        /// Reason exact comparison is valid.
        rationale: String,
    },
    /// Compare using absolute plus relative tolerances.
    FloatingPoint {
        /// Numeric representation used by the fixture.
        dtype: String,
        /// Quantity being compared.
        quantity: String,
        /// Absolute tolerance.
        absolute_tolerance: f64,
        /// Relative tolerance.
        relative_tolerance: f64,
        /// Reason the tolerances are appropriate.
        rationale: String,
    },
}

impl ComparisonPolicy {
    fn validate(&self) -> Result<(), FixtureError> {
        match self {
            Self::Exact {
                quantity,
                rationale,
            } => {
                ensure_nonempty(quantity, "comparison quantity")?;
                ensure_nonempty(rationale, "comparison rationale")
            }
            Self::FloatingPoint {
                dtype,
                quantity,
                absolute_tolerance,
                relative_tolerance,
                rationale,
            } => {
                ensure(dtype == "f64", "only f64 fixture comparisons are supported")?;
                ensure_nonempty(quantity, "comparison quantity")?;
                ensure_nonempty(rationale, "comparison rationale")?;
                ensure(
                    absolute_tolerance.is_finite() && *absolute_tolerance > 0.0,
                    "absolute tolerance must be finite and positive",
                )?;
                ensure(
                    relative_tolerance.is_finite() && *relative_tolerance >= 0.0,
                    "relative tolerance must be finite and non-negative",
                )
            }
        }
    }

    fn tolerance(&self) -> f64 {
        match self {
            Self::Exact { .. } => 0.0,
            Self::FloatingPoint {
                absolute_tolerance, ..
            } => *absolute_tolerance,
        }
    }

    /// Return whether the comparison is exact.
    pub const fn is_exact(&self) -> bool {
        matches!(self, Self::Exact { .. })
    }

    /// Return the absolute tolerance, or zero for exact comparisons.
    pub const fn absolute_tolerance(&self) -> f64 {
        match self {
            Self::Exact { .. } => 0.0,
            Self::FloatingPoint {
                absolute_tolerance, ..
            } => *absolute_tolerance,
        }
    }

    /// Return the relative tolerance, or zero for exact comparisons.
    pub const fn relative_tolerance(&self) -> f64 {
        match self {
            Self::Exact { .. } => 0.0,
            Self::FloatingPoint {
                relative_tolerance, ..
            } => *relative_tolerance,
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FixtureManifest {
    manifest_schema: String,
    fixture_schema: String,
    convention_schema: String,
    digest_algorithm: String,
    entries: Vec<ManifestEntry>,
}

impl FixtureManifest {
    fn validate(&self) -> Result<(), FixtureError> {
        ensure(
            self.manifest_schema == MANIFEST_SCHEMA,
            "unsupported manifest schema",
        )?;
        ensure(
            self.fixture_schema == FIXTURE_SCHEMA,
            "manifest fixture schema mismatch",
        )?;
        ensure(
            self.convention_schema == CONVENTION_SCHEMA,
            "manifest convention schema mismatch",
        )?;
        ensure(
            self.digest_algorithm == "blake3",
            "unsupported digest algorithm",
        )?;
        ensure(!self.entries.is_empty(), "fixture manifest is empty")?;
        let mut previous = None;
        for entry in &self.entries {
            ensure(
                entry.path.ends_with(".json")
                    && !entry.path.contains('/')
                    && !entry.path.contains('\\'),
                "fixture manifest paths must be plain JSON file names",
            )?;
            if let Some(previous) = previous {
                ensure(
                    previous < entry.path.as_str(),
                    "manifest paths must be sorted",
                )?;
            }
            previous = Some(entry.path.as_str());
            ensure_nonempty(&entry.id, "manifest fixture id")?;
            ensure(
                entry.blake3.len() == 64
                    && entry.blake3.bytes().all(|byte| byte.is_ascii_hexdigit()),
                "manifest BLAKE3 digest must contain 64 hexadecimal characters",
            )?;
        }
        Ok(())
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ManifestEntry {
    path: String,
    id: String,
    kind: FixtureKind,
    blake3: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CoordinateSite {
    coordinate: [u32; 2],
    site: u32,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RectangularData {
    lx: u32,
    ly: u32,
    coordinate_to_site: Vec<CoordinateSite>,
    open_bonds: Vec<[u32; 2]>,
    periodic_bonds: Vec<[u32; 2]>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BitPackingData {
    bits: Vec<u8>,
    packed_mask: u64,
    pauli_z_values: Vec<i8>,
    rydberg_occupancies: Vec<u8>,
    hamming_weight: u32,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ComplexEntry {
    re: f64,
    im: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Matrix {
    shape: [usize; 2],
    dtype: String,
    layout: String,
    entry_definition: String,
    entries: Vec<ComplexEntry>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WeightedInteraction {
    interaction_id: String,
    sites: Vec<u32>,
    channel: String,
    coefficient: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TfimData {
    site_count: u32,
    bond: [u32; 2],
    j: f64,
    fields: Vec<f64>,
    weighted_interactions: Vec<WeightedInteraction>,
    matrix: Matrix,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct HeisenbergData {
    site_count: u32,
    bond: [u32; 2],
    j: f64,
    weighted_interactions: Vec<WeightedInteraction>,
    matrix: Matrix,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WeightedBond {
    sites: [u32; 2],
    j: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Connection {
    output_bits: Vec<u8>,
    output_mask: u64,
    matrix_element: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct HeterogeneousData {
    site_count: u32,
    input_bits: Vec<u8>,
    input_mask: u64,
    weighted_bonds: Vec<WeightedBond>,
    diagonal_energy: f64,
    connections: Vec<Connection>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RydbergData {
    site_count: u32,
    omega: f64,
    detunings: Vec<f64>,
    interaction_v: f64,
    weighted_interactions: Vec<WeightedInteraction>,
    matrix: Matrix,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ObservableData {
    site_count: u32,
    bits: Vec<u8>,
    sigma_z_per_site: f64,
    raw_zz_for_every_pair: f64,
    connected_zz_for_every_pair: f64,
    uniform_pauli_generator_fisher_density: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BasisRotationData {
    site_count: u32,
    j: f64,
    h: f64,
    weighted_interactions: Vec<WeightedInteraction>,
    z_basis_matrix: Matrix,
    x_basis_matrix: Matrix,
    sorted_spectrum: Vec<f64>,
}

fn validate_case(kind: FixtureKind, data: &Value, tolerance: f64) -> Result<(), FixtureError> {
    match kind {
        FixtureKind::RectangularIndexing => validate_rectangle(decode(data)?),
        FixtureKind::BitPacking => validate_bit_packing(decode(data)?),
        FixtureKind::TfimOneBond => validate_tfim(decode(data)?, tolerance),
        FixtureKind::HeisenbergOneBond => validate_heisenberg(decode(data)?, tolerance),
        FixtureKind::HeisenbergHeterogeneous => validate_heterogeneous(decode(data)?),
        FixtureKind::RydbergTwoSite => validate_rydberg(decode(data)?, tolerance),
        FixtureKind::ObservableNormalization => validate_observable(decode(data)?),
        FixtureKind::BasisRotationSpectrum => validate_basis_rotation(decode(data)?, tolerance),
    }
}

fn decode<T: DeserializeOwned>(data: &Value) -> Result<T, FixtureError> {
    serde_json::from_value(data.clone())
        .map_err(|error| FixtureError::new(format!("invalid case payload: {error}")))
}

fn validate_rectangle(data: RectangularData) -> Result<(), FixtureError> {
    ensure(
        data.lx > 0 && data.ly > 0,
        "rectangle dimensions must be positive",
    )?;
    let count = data
        .lx
        .checked_mul(data.ly)
        .ok_or_else(|| FixtureError::new("rectangle site count overflow"))?;
    ensure(
        data.coordinate_to_site.len() == count as usize,
        "coordinate map length differs from Lx*Ly",
    )?;
    for entry in &data.coordinate_to_site {
        ensure(
            entry.coordinate[0] < data.lx && entry.coordinate[1] < data.ly,
            "coordinate lies outside the rectangle",
        )?;
        ensure(entry.site < count, "site id lies outside the rectangle")?;
    }
    validate_bonds(&data.open_bonds, count)?;
    validate_bonds(&data.periodic_bonds, count)
}

fn validate_bonds(bonds: &[[u32; 2]], site_count: u32) -> Result<(), FixtureError> {
    let mut previous = None;
    for bond in bonds {
        ensure(
            bond[0] < bond[1] && bond[1] < site_count,
            "bond endpoints must be canonical distinct valid sites",
        )?;
        if let Some(previous) = previous {
            ensure(
                previous < *bond,
                "bonds must be unique and lexicographically sorted",
            )?;
        }
        previous = Some(*bond);
    }
    Ok(())
}

fn validate_bit_packing(data: BitPackingData) -> Result<(), FixtureError> {
    validate_bits(&data.bits)?;
    ensure(
        data.bits.len() == data.pauli_z_values.len()
            && data.bits.len() == data.rydberg_occupancies.len(),
        "local-value vectors must match the bit count",
    )?;
    let mut mask = 0_u64;
    for (site, bit) in data.bits.iter().copied().enumerate() {
        ensure(site < u64::BITS as usize, "test mask exceeds u64 width")?;
        mask |= u64::from(bit) << site;
        ensure(
            data.pauli_z_values[site] == 1 - 2 * bit as i8,
            "Pauli-Z value does not match the bit",
        )?;
        ensure(
            data.rydberg_occupancies[site] == bit,
            "Rydberg occupation does not match the bit",
        )?;
    }
    ensure(mask == data.packed_mask, "packed mask does not match bits")?;
    let weight = data.bits.iter().map(|bit| u32::from(*bit)).sum::<u32>();
    ensure(
        weight == data.hamming_weight,
        "Hamming weight does not match bits",
    )
}

fn validate_tfim(data: TfimData, tolerance: f64) -> Result<(), FixtureError> {
    ensure(
        data.site_count == 2,
        "one-bond TFIM fixture must have two sites",
    )?;
    ensure(
        data.bond == [0, 1],
        "one-bond TFIM fixture must use bond 0-1",
    )?;
    ensure(data.j.is_finite(), "TFIM coupling must be finite")?;
    validate_finite(&data.fields, "TFIM fields")?;
    ensure(
        data.fields.len() == 2,
        "TFIM fixture requires one field per site",
    )?;
    validate_interactions(&data.weighted_interactions, 2)?;
    validate_matrix(&data.matrix, 4, tolerance)
}

fn validate_heisenberg(data: HeisenbergData, tolerance: f64) -> Result<(), FixtureError> {
    ensure(
        data.site_count == 2,
        "one-bond Heisenberg fixture must have two sites",
    )?;
    ensure(data.bond == [0, 1], "Heisenberg fixture must use bond 0-1")?;
    ensure(data.j.is_finite(), "Heisenberg coupling must be finite")?;
    validate_interactions(&data.weighted_interactions, 2)?;
    validate_matrix(&data.matrix, 4, tolerance)
}

fn validate_heterogeneous(data: HeterogeneousData) -> Result<(), FixtureError> {
    ensure(
        data.site_count == 3,
        "heterogeneous fixture must have three sites",
    )?;
    validate_bits(&data.input_bits)?;
    ensure(
        data.input_bits.len() == 3,
        "heterogeneous input must have three bits",
    )?;
    ensure(
        data.input_mask < 8,
        "heterogeneous input mask is out of range",
    )?;
    ensure(
        data.diagonal_energy.is_finite(),
        "diagonal energy must be finite",
    )?;
    ensure(
        data.weighted_bonds.len() >= 2,
        "fixture must contain heterogeneous bonds",
    )?;
    let mut coefficients = BTreeSet::new();
    for bond in &data.weighted_bonds {
        ensure(
            bond.sites[0] < bond.sites[1] && bond.sites[1] < data.site_count,
            "weighted bond sites are invalid",
        )?;
        ensure(
            bond.j.is_finite(),
            "weighted bond coefficient must be finite",
        )?;
        coefficients.insert(bond.j.to_bits());
    }
    ensure(
        coefficients.len() > 1,
        "fixture couplings must be heterogeneous",
    )?;
    for connection in &data.connections {
        validate_bits(&connection.output_bits)?;
        ensure(
            connection.output_bits.len() == 3,
            "connection must have three bits",
        )?;
        ensure(
            connection.output_mask < 8,
            "connection mask is out of range",
        )?;
        ensure(
            connection.matrix_element.is_finite(),
            "matrix element must be finite",
        )?;
    }
    Ok(())
}

fn validate_rydberg(data: RydbergData, tolerance: f64) -> Result<(), FixtureError> {
    ensure(data.site_count == 2, "Rydberg fixture must have two sites")?;
    ensure(data.omega.is_finite(), "Rydberg drive must be finite")?;
    ensure(
        data.interaction_v.is_finite(),
        "Rydberg interaction must be finite",
    )?;
    validate_finite(&data.detunings, "Rydberg detunings")?;
    ensure(
        data.detunings.len() == 2,
        "Rydberg fixture requires one detuning per site",
    )?;
    validate_interactions(&data.weighted_interactions, 2)?;
    validate_matrix(&data.matrix, 4, tolerance)
}

fn validate_observable(data: ObservableData) -> Result<(), FixtureError> {
    ensure(data.site_count > 0, "observable fixture must contain sites")?;
    validate_bits(&data.bits)?;
    ensure(
        data.bits.len() == data.site_count as usize,
        "observable bit count differs from site count",
    )?;
    validate_finite(
        &[
            data.sigma_z_per_site,
            data.raw_zz_for_every_pair,
            data.connected_zz_for_every_pair,
            data.uniform_pauli_generator_fisher_density,
        ],
        "observable values",
    )
}

fn validate_basis_rotation(data: BasisRotationData, tolerance: f64) -> Result<(), FixtureError> {
    ensure(
        data.site_count == 2,
        "basis-rotation fixture must have two sites",
    )?;
    validate_finite(&[data.j, data.h], "basis-rotation parameters")?;
    validate_interactions(&data.weighted_interactions, 2)?;
    validate_matrix(&data.z_basis_matrix, 4, tolerance)?;
    validate_matrix(&data.x_basis_matrix, 4, tolerance)?;
    validate_finite(&data.sorted_spectrum, "reference spectrum")?;
    ensure(
        data.sorted_spectrum.len() == 4,
        "reference spectrum must have four values",
    )?;
    ensure(
        data.sorted_spectrum
            .windows(2)
            .all(|pair| pair[0] <= pair[1]),
        "reference spectrum must be sorted",
    )
}

fn validate_interactions(
    interactions: &[WeightedInteraction],
    site_count: u32,
) -> Result<(), FixtureError> {
    ensure(
        !interactions.is_empty(),
        "resolved interaction table is empty",
    )?;
    let mut identities = BTreeSet::new();
    for interaction in interactions {
        ensure_nonempty(&interaction.interaction_id, "interaction id")?;
        ensure_nonempty(&interaction.channel, "interaction channel")?;
        ensure(
            interaction.coefficient.is_finite(),
            "interaction coefficient must be finite",
        )?;
        ensure(
            !interaction.sites.is_empty(),
            "interaction support is empty",
        )?;
        ensure(
            interaction.sites.iter().all(|site| *site < site_count),
            "interaction site is out of range",
        )?;
        ensure(
            identities.insert(interaction.interaction_id.as_str()),
            "duplicate interaction identity",
        )?;
    }
    Ok(())
}

fn validate_matrix(
    matrix: &Matrix,
    expected_dimension: usize,
    tolerance: f64,
) -> Result<(), FixtureError> {
    ensure(
        matrix.shape == [expected_dimension, expected_dimension],
        "matrix shape does not match the basis dimension",
    )?;
    ensure(
        matrix.dtype == "complex_f64",
        "matrix dtype must be complex_f64",
    )?;
    ensure(
        matrix.layout == "row_major",
        "matrix layout must be row_major",
    )?;
    ensure_nonempty(&matrix.entry_definition, "matrix entry definition")?;
    ensure(
        matrix.entries.len() == expected_dimension * expected_dimension,
        "matrix entry count does not match its shape",
    )?;
    for entry in &matrix.entries {
        ensure(
            entry.re.is_finite() && entry.im.is_finite(),
            "matrix entry must be finite",
        )?;
    }
    for row in 0..expected_dimension {
        for column in 0..expected_dimension {
            let value = &matrix.entries[row * expected_dimension + column];
            let adjoint = &matrix.entries[column * expected_dimension + row];
            ensure(
                (value.re - adjoint.re).abs() <= tolerance
                    && (value.im + adjoint.im).abs() <= tolerance,
                "matrix is not Hermitian within its declared tolerance",
            )?;
        }
    }
    Ok(())
}

fn validate_bits(bits: &[u8]) -> Result<(), FixtureError> {
    ensure(!bits.is_empty(), "bit vector must not be empty")?;
    ensure(
        bits.iter().all(|bit| *bit <= 1),
        "bit vector contains a non-binary value",
    )
}

fn validate_finite(values: &[f64], name: &str) -> Result<(), FixtureError> {
    ensure(
        values.iter().all(|value| value.is_finite()),
        format!("{name} must be finite"),
    )
}

fn validate_top_level_keys(value: &Value) -> Result<(), FixtureError> {
    let object = value
        .as_object()
        .ok_or_else(|| FixtureError::new("fixture must be a JSON object"))?;
    const EXPECTED: [&str; 10] = [
        "authorship",
        "claim",
        "comparison",
        "convention_schema",
        "conventions",
        "data",
        "fixture_schema",
        "id",
        "kind",
        "oracle",
    ];
    let actual: BTreeSet<_> = object.keys().map(String::as_str).collect();
    let expected: BTreeSet<_> = EXPECTED.into_iter().collect();
    ensure(
        actual == expected,
        "fixture has missing or unknown top-level fields",
    )
}

fn validate_date(date: &str) -> Result<(), FixtureError> {
    let bytes = date.as_bytes();
    ensure(
        bytes.len() == 10
            && bytes[4] == b'-'
            && bytes[7] == b'-'
            && bytes
                .iter()
                .enumerate()
                .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit()),
        "review date must use YYYY-MM-DD",
    )
}

fn ensure_nonempty(value: &str, name: &str) -> Result<(), FixtureError> {
    ensure(
        !value.trim().is_empty(),
        format!("{name} must not be empty"),
    )
}

fn ensure_nonempty_object(value: &Value, name: &str) -> Result<(), FixtureError> {
    ensure(
        value.as_object().is_some_and(|object| !object.is_empty()),
        format!("{name} must be a non-empty object"),
    )
}

fn ensure(condition: bool, message: impl Into<String>) -> Result<(), FixtureError> {
    if condition {
        Ok(())
    } else {
        Err(FixtureError::new(message))
    }
}
