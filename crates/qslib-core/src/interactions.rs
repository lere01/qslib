use crate::{Bond, Real, SiteCount, SiteId, derive_seed};
use blake3::Hasher;
use rand_chacha::ChaCha20Rng;
use rand_core::{Rng, SeedableRng};
use std::fmt::{self, Display, Formatter};

/// Errors raised while resolving weighted interactions and coupling inputs.
#[derive(Clone, Debug, PartialEq)]
pub enum InteractionError {
    /// A coefficient or disorder bound was non-finite.
    NonFiniteCoefficient {
        /// Supplied non-finite value.
        value: Real,
    },
    /// Disorder bounds were reversed.
    InvalidDisorderBounds {
        /// Lower bound.
        lower: Real,
        /// Upper bound.
        upper: Real,
    },
    /// An interaction endpoint is outside the declared system.
    SiteOutOfRange {
        /// Invalid endpoint.
        site: SiteId,
        /// Declared number of sites.
        site_count: usize,
    },
    /// Two terms have the same endpoint, channel, and image identity.
    DuplicateIdentity {
        /// Repeated bond identity.
        bond: Bond,
        /// Repeated channel.
        channel: InteractionChannel,
    },
    /// A dense coupling array is not square with the declared site count.
    InvalidDenseShape {
        /// Required number of entries.
        expected: usize,
        /// Supplied number of entries.
        actual: usize,
    },
    /// A dense undirected coupling matrix is not symmetric.
    NonSymmetric {
        /// Row index.
        row: usize,
        /// Column index.
        column: usize,
        /// Upper-triangle value.
        upper: Real,
        /// Lower-triangle value.
        lower: Real,
    },
    /// A dense matrix has a non-zero diagonal.
    NonZeroDiagonal {
        /// Diagonal index.
        index: usize,
        /// Supplied diagonal value.
        value: Real,
    },
    /// A sparse table repeats an unordered pair.
    DuplicatePair {
        /// Canonical first endpoint.
        first: SiteId,
        /// Canonical second endpoint.
        second: SiteId,
    },
    /// A generic channel name is empty or not a stable identifier.
    InvalidChannel {
        /// Invalid channel value.
        value: String,
    },
    /// A symmetry tolerance is negative or non-finite.
    InvalidTolerance {
        /// Invalid tolerance value.
        value: Real,
    },
}

impl Display for InteractionError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteCoefficient { value } => {
                write!(formatter, "interaction coefficient {value:?} is non-finite")
            }
            Self::InvalidDisorderBounds { lower, upper } => {
                write!(formatter, "disorder bounds [{lower}, {upper}] are invalid")
            }
            Self::SiteOutOfRange { site, site_count } => write!(
                formatter,
                "site {} is outside a {site_count}-site interaction table",
                site.get()
            ),
            Self::DuplicateIdentity { bond, channel } => write!(
                formatter,
                "duplicate interaction identity ({}, {}) on {channel:?}",
                bond.first().get(),
                bond.second().get()
            ),
            Self::InvalidDenseShape { expected, actual } => write!(
                formatter,
                "dense coupling matrix needs {expected} entries, received {actual}"
            ),
            Self::NonSymmetric {
                row,
                column,
                upper,
                lower,
            } => write!(
                formatter,
                "dense coupling matrix differs at ({row}, {column}): {upper} versus {lower}"
            ),
            Self::NonZeroDiagonal { index, value } => write!(
                formatter,
                "dense coupling diagonal {index} is {value}, not zero"
            ),
            Self::DuplicatePair { first, second } => write!(
                formatter,
                "sparse coupling pair ({}, {}) is duplicated",
                first.get(),
                second.get()
            ),
            Self::InvalidChannel { value } => write!(
                formatter,
                "interaction channel {value:?} is not a stable non-empty identifier"
            ),
            Self::InvalidTolerance { value } => {
                write!(formatter, "symmetry tolerance {value:?} is invalid")
            }
        }
    }
}

impl std::error::Error for InteractionError {}

/// Operator channel attached to one weighted interaction.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum InteractionChannel {
    /// Ising `Z_i Z_j` bond channel.
    IsingZZ,
    /// Isotropic Heisenberg exchange channel.
    HeisenbergExchange,
    /// Rydberg occupation-density pair channel.
    RydbergDensityDensity,
    /// Generic named channel for model-specific operators.
    Generic(String),
}

impl InteractionChannel {
    /// Construct a generic channel from a stable non-empty identifier.
    pub fn generic(value: impl Into<String>) -> Result<Self, InteractionError> {
        let value = value.into();
        if value.trim().is_empty() || value.chars().any(char::is_whitespace) {
            return Err(InteractionError::InvalidChannel { value });
        }
        Ok(Self::Generic(value))
    }
}

/// Full identity of a physical interaction term.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct InteractionIdentity {
    bond: Bond,
    channel: InteractionChannel,
    name: Option<String>,
}

impl InteractionIdentity {
    /// Construct an unnamed interaction identity.
    pub fn new(bond: Bond, channel: InteractionChannel) -> Self {
        Self {
            bond,
            channel,
            name: None,
        }
    }

    /// Construct a named identity for overlapping shells or physical terms.
    pub fn named(
        bond: Bond,
        channel: InteractionChannel,
        name: impl Into<String>,
    ) -> Result<Self, InteractionError> {
        let name = name.into();
        if name.trim().is_empty() || name.chars().any(char::is_whitespace) {
            return Err(InteractionError::InvalidChannel { value: name });
        }
        Ok(Self {
            bond,
            channel,
            name: Some(name),
        })
    }

    /// Return the canonical bond.
    pub const fn bond(&self) -> Bond {
        self.bond
    }

    /// Return the operator channel.
    pub fn channel(&self) -> &InteractionChannel {
        &self.channel
    }

    /// Return the optional stable physical-term name.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
}

/// A resolved coefficient attached to one canonical interaction identity.
#[derive(Clone, Debug, PartialEq)]
pub struct WeightedInteraction {
    identity: InteractionIdentity,
    coefficient: Real,
}

impl WeightedInteraction {
    /// Construct a finite weighted interaction.
    pub fn new(
        bond: Bond,
        channel: InteractionChannel,
        coefficient: Real,
    ) -> Result<Self, InteractionError> {
        Self::from_identity(InteractionIdentity::new(bond, channel), coefficient)
    }

    /// Construct a finite interaction from its complete identity.
    pub fn from_identity(
        identity: InteractionIdentity,
        coefficient: Real,
    ) -> Result<Self, InteractionError> {
        if let InteractionChannel::Generic(value) = identity.channel() {
            if value.trim().is_empty() || value.chars().any(char::is_whitespace) {
                return Err(InteractionError::InvalidChannel {
                    value: value.clone(),
                });
            }
        }
        if !coefficient.is_finite() {
            return Err(InteractionError::NonFiniteCoefficient { value: coefficient });
        }
        Ok(Self {
            identity,
            coefficient,
        })
    }

    /// Return the canonical bond identity.
    pub const fn bond(&self) -> Bond {
        self.identity.bond()
    }

    /// Return the operator channel.
    pub fn channel(&self) -> &InteractionChannel {
        self.identity.channel()
    }

    /// Return the complete interaction identity.
    pub fn identity(&self) -> &InteractionIdentity {
        &self.identity
    }

    /// Return the resolved coefficient.
    pub const fn coefficient(&self) -> Real {
        self.coefficient
    }
}

/// Deterministic collection of resolved weighted interactions.
#[derive(Clone, Debug, PartialEq)]
pub struct InteractionTable {
    site_count: SiteCount,
    interactions: Vec<WeightedInteraction>,
}

impl InteractionTable {
    /// Validate, canonicalize, and sort weighted interactions.
    pub fn new(
        site_count: SiteCount,
        terms: Vec<(Bond, InteractionChannel, Real)>,
    ) -> Result<Self, InteractionError> {
        let identities = terms
            .into_iter()
            .map(|(bond, channel, coefficient)| {
                (InteractionIdentity::new(bond, channel), coefficient)
            })
            .collect();
        Self::new_with_identities(site_count, identities)
    }

    /// Validate, canonicalize, and sort terms with named identities.
    pub fn new_with_identities(
        site_count: SiteCount,
        terms: Vec<(InteractionIdentity, Real)>,
    ) -> Result<Self, InteractionError> {
        let mut interactions = Vec::with_capacity(terms.len());
        for (identity, coefficient) in terms {
            let bond = identity.bond();
            if let InteractionChannel::Generic(value) = identity.channel() {
                if value.trim().is_empty() || value.chars().any(char::is_whitespace) {
                    return Err(InteractionError::InvalidChannel {
                        value: value.clone(),
                    });
                }
            }
            site_count
                .validate(bond.first())
                .map_err(|_| InteractionError::SiteOutOfRange {
                    site: bond.first(),
                    site_count: site_count.get(),
                })?;
            site_count
                .validate(bond.second())
                .map_err(|_| InteractionError::SiteOutOfRange {
                    site: bond.second(),
                    site_count: site_count.get(),
                })?;
            interactions.push(WeightedInteraction::from_identity(identity, coefficient)?);
        }
        interactions.sort_by(|a, b| a.identity().cmp(b.identity()));
        for pair in interactions.windows(2) {
            if pair[0].identity() == pair[1].identity() {
                return Err(InteractionError::DuplicateIdentity {
                    bond: pair[0].bond(),
                    channel: pair[0].channel().clone(),
                });
            }
        }
        Ok(Self {
            site_count,
            interactions,
        })
    }

    /// Return the declared site count.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Return canonical resolved interactions.
    pub fn interactions(&self) -> &[WeightedInteraction] {
        &self.interactions
    }

    /// Return all declared terms, including exact-zero provenance entries.
    pub fn declared_interactions(&self) -> &[WeightedInteraction] {
        &self.interactions
    }

    /// Return terms that contribute numerically to an evaluated Hamiltonian.
    pub fn active_interactions(&self) -> Vec<&WeightedInteraction> {
        self.interactions
            .iter()
            .filter(|term| term.coefficient() != 0.0)
            .collect()
    }

    /// Generate a deterministic realized uniform disorder table from identities.
    pub fn realize_uniform_disorder(
        &self,
        seed: [u8; 32],
        lower: Real,
        upper: Real,
    ) -> Result<DisorderRealization, InteractionError> {
        self.realize_uniform_disorder_at(seed, 0, lower, upper)
    }

    /// Generate a realization with an explicit logical realization index.
    pub fn realize_uniform_disorder_at(
        &self,
        seed: [u8; 32],
        realization_index: u64,
        lower: Real,
        upper: Real,
    ) -> Result<DisorderRealization, InteractionError> {
        if !lower.is_finite() || !upper.is_finite() || lower > upper {
            return Err(InteractionError::InvalidDisorderBounds { lower, upper });
        }
        let interactions = self
            .interactions
            .iter()
            .map(|term| {
                let identity_digest = canonical_identity_digest(term);
                let indices = [
                    realization_index,
                    term.bond().first().get() as u64,
                    term.bond().second().get() as u64,
                    term.bond().image_translation().0 as i64 as u64,
                    term.bond().image_translation().1 as i64 as u64,
                    term.bond().source().get() as u64,
                    term.bond().direction().0 as i64 as u64,
                    term.bond().direction().1 as i64 as u64,
                    digest_word(&identity_digest[0..8]),
                    digest_word(&identity_digest[8..16]),
                    digest_word(&identity_digest[16..24]),
                    digest_word(&identity_digest[24..32]),
                ];
                let child_seed = derive_seed(&seed, "disorder", &indices);
                let mut rng = ChaCha20Rng::from_seed(child_seed);
                let unit = rng.next_u64() as Real / u64::MAX as Real;
                WeightedInteraction {
                    identity: term.identity().clone(),
                    coefficient: lower + (upper - lower) * unit,
                }
            })
            .collect();
        Ok(DisorderRealization {
            site_count: self.site_count,
            interactions,
            provenance: DisorderProvenance {
                seed,
                distribution: format!("uniform[{lower},{upper}]"),
                lower,
                upper,
                semantics: "replacement".to_owned(),
                rng_algorithm: "chacha20".to_owned(),
                seed_scheme: "qslib-seed-v1".to_owned(),
                domain: "disorder".to_owned(),
                realization_index,
            },
        })
    }
}

fn canonical_identity_digest(term: &WeightedInteraction) -> [u8; 32] {
    let mut hasher = Hasher::new();
    hasher.update(b"qslib-interaction-identity-v1\0");
    hasher.update(&term.bond().first().get().to_le_bytes());
    hasher.update(&term.bond().second().get().to_le_bytes());
    hasher.update(&term.bond().image_translation().0.to_le_bytes());
    hasher.update(&term.bond().image_translation().1.to_le_bytes());
    hasher.update(&term.bond().source().get().to_le_bytes());
    hasher.update(&term.bond().direction().0.to_le_bytes());
    hasher.update(&term.bond().direction().1.to_le_bytes());
    match term.channel() {
        InteractionChannel::IsingZZ => {
            hasher.update(b"ising_zz");
        }
        InteractionChannel::HeisenbergExchange => {
            hasher.update(b"heisenberg_exchange");
        }
        InteractionChannel::RydbergDensityDensity => {
            hasher.update(b"rydberg_density_density");
        }
        InteractionChannel::Generic(value) => {
            hasher.update(b"generic");
            hasher.update(&(value.len() as u32).to_le_bytes());
            hasher.update(value.as_bytes());
        }
    }
    let name = term.identity().name().unwrap_or("");
    hasher.update(&(name.len() as u32).to_le_bytes());
    hasher.update(name.as_bytes());
    *hasher.finalize().as_bytes()
}

fn digest_word(bytes: &[u8]) -> u64 {
    let mut raw = [0_u8; 8];
    raw.copy_from_slice(bytes);
    u64::from_le_bytes(raw)
}

/// A validated dense symmetric coupling matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct DenseCouplings {
    site_count: SiteCount,
    values: Vec<Real>,
}

impl DenseCouplings {
    /// Validate a row-major `[N, N]` matrix with zero diagonal.
    pub fn new(site_count: SiteCount, values: Vec<Real>) -> Result<Self, InteractionError> {
        Self::new_with_tolerance(site_count, values, 1.0e-12)
    }

    /// Validate with an explicit absolute symmetry tolerance.
    pub fn new_with_tolerance(
        site_count: SiteCount,
        values: Vec<Real>,
        tolerance: Real,
    ) -> Result<Self, InteractionError> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(InteractionError::InvalidTolerance { value: tolerance });
        }
        let n = site_count.get();
        let expected = n
            .checked_mul(n)
            .ok_or(InteractionError::InvalidDenseShape {
                expected: usize::MAX,
                actual: values.len(),
            })?;
        if values.len() != expected {
            return Err(InteractionError::InvalidDenseShape {
                expected,
                actual: values.len(),
            });
        }
        for i in 0..n {
            let diagonal = values[i * n + i];
            if !diagonal.is_finite() {
                return Err(InteractionError::NonFiniteCoefficient { value: diagonal });
            }
            if diagonal != 0.0 {
                return Err(InteractionError::NonZeroDiagonal {
                    index: i,
                    value: diagonal,
                });
            }
            for j in (i + 1)..n {
                let upper = values[i * n + j];
                let lower = values[j * n + i];
                if !upper.is_finite() {
                    return Err(InteractionError::NonFiniteCoefficient { value: upper });
                }
                if !lower.is_finite() {
                    return Err(InteractionError::NonFiniteCoefficient { value: lower });
                }
                if (upper - lower).abs() > tolerance {
                    return Err(InteractionError::NonSymmetric {
                        row: i,
                        column: j,
                        upper,
                        lower,
                    });
                }
            }
        }
        Ok(Self { site_count, values })
    }

    /// Return the number of sites represented by the matrix.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Resolve upper-triangle entries into canonical interactions.
    pub fn to_interactions(
        &self,
        channel: InteractionChannel,
    ) -> Result<Vec<WeightedInteraction>, InteractionError> {
        let n = self.site_count.get();
        let mut result = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let bond = Bond::new(
                    SiteId::try_from_usize(i).map_err(|_| InteractionError::SiteOutOfRange {
                        site: SiteId::new(u32::MAX),
                        site_count: n,
                    })?,
                    SiteId::try_from_usize(j).map_err(|_| InteractionError::SiteOutOfRange {
                        site: SiteId::new(u32::MAX),
                        site_count: n,
                    })?,
                )
                .map_err(|_| InteractionError::SiteOutOfRange {
                    site: SiteId::new(u32::MAX),
                    site_count: n,
                })?;
                result.push(WeightedInteraction::new(
                    bond,
                    channel.clone(),
                    self.values[i * n + j],
                )?);
            }
        }
        Ok(result)
    }
}

/// A validated sparse pair coupling table.
#[derive(Clone, Debug, PartialEq)]
pub struct SparseCouplings {
    site_count: SiteCount,
    pairs: Vec<(Bond, Real)>,
}

impl SparseCouplings {
    /// Validate and canonicalize sparse unordered pairs.
    pub fn new(
        site_count: SiteCount,
        entries: Vec<(SiteId, SiteId, Real)>,
    ) -> Result<Self, InteractionError> {
        let mut pairs = Vec::with_capacity(entries.len());
        for (a, b, coefficient) in entries {
            site_count
                .validate(a)
                .map_err(|_| InteractionError::SiteOutOfRange {
                    site: a,
                    site_count: site_count.get(),
                })?;
            site_count
                .validate(b)
                .map_err(|_| InteractionError::SiteOutOfRange {
                    site: b,
                    site_count: site_count.get(),
                })?;
            let bond = Bond::new(a, b).map_err(|_| InteractionError::DuplicatePair {
                first: a,
                second: b,
            })?;
            if !coefficient.is_finite() {
                return Err(InteractionError::NonFiniteCoefficient { value: coefficient });
            }
            pairs.push((bond, coefficient));
        }
        pairs.sort_by_key(|entry| entry.0);
        if let Some(window) = pairs.windows(2).find(|pair| pair[0].0 == pair[1].0) {
            return Err(InteractionError::DuplicatePair {
                first: window[0].0.first(),
                second: window[0].0.second(),
            });
        }
        Ok(Self { site_count, pairs })
    }

    /// Resolve sparse entries into weighted interactions.
    pub fn to_interactions(
        &self,
        channel: InteractionChannel,
    ) -> Result<Vec<WeightedInteraction>, InteractionError> {
        self.pairs
            .iter()
            .map(|(bond, coefficient)| {
                WeightedInteraction::new(*bond, channel.clone(), *coefficient)
            })
            .collect()
    }
}

/// Durable provenance for a realized disorder table.
#[derive(Clone, Debug, PartialEq)]
pub struct DisorderProvenance {
    seed: [u8; 32],
    distribution: String,
    lower: Real,
    upper: Real,
    semantics: String,
    rng_algorithm: String,
    seed_scheme: String,
    domain: String,
    realization_index: u64,
}

impl DisorderProvenance {
    /// Return the master seed used to derive the realization.
    pub const fn seed(&self) -> [u8; 32] {
        self.seed
    }

    /// Return the named distribution and bounds.
    pub fn distribution(&self) -> &str {
        &self.distribution
    }

    /// Return the lower bound used by the realization.
    pub const fn lower(&self) -> Real {
        self.lower
    }

    /// Return the upper bound used by the realization.
    pub const fn upper(&self) -> Real {
        self.upper
    }

    /// Return the coefficient semantics, currently `replacement`.
    pub fn semantics(&self) -> &str {
        &self.semantics
    }

    /// Return the deterministic RNG algorithm identifier.
    pub fn rng_algorithm(&self) -> &str {
        &self.rng_algorithm
    }

    /// Return the child-seed derivation schema identifier.
    pub fn seed_scheme(&self) -> &str {
        &self.seed_scheme
    }

    /// Return the logical derivation domain.
    pub fn domain(&self) -> &str {
        &self.domain
    }

    /// Return the logical realization index.
    pub const fn realization_index(&self) -> u64 {
        self.realization_index
    }
}

/// A realized, persisted per-interaction disorder table.
#[derive(Clone, Debug, PartialEq)]
pub struct DisorderRealization {
    site_count: SiteCount,
    interactions: Vec<WeightedInteraction>,
    provenance: DisorderProvenance,
}

impl DisorderRealization {
    /// Return the validated source topology's site count.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Return resolved coefficients; this table is authoritative for reruns.
    pub fn interactions(&self) -> &[WeightedInteraction] {
        &self.interactions
    }

    /// Return generation provenance retained alongside the realized values.
    pub const fn provenance(&self) -> &DisorderProvenance {
        &self.provenance
    }
}
