//! Canonical qslib SSE terms and sign-safe TFIM/Rydberg decompositions.

use qslib_core::{BasisBit, InteractionChannel, WeightedInteraction};
use std::fmt::{self, Display, Formatter};

/// Identity, diagonal, or off-diagonal classification of an SSE vertex.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OperatorKind {
    /// Unused padded position.
    Identity,
    /// Vertex diagonal in the canonical bit basis.
    Diagonal,
    /// Vertex that changes one or more basis bits.
    OffDiagonal,
}

/// A non-negative local operator in the expansion `H=E_shift-sum(B_a)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SseTerm {
    /// `B = J (1 + z_i z_j)` for TFIM.
    TfimBond {
        /// First endpoint.
        site_i: u32,
        /// Second endpoint.
        site_j: u32,
        /// Coupling.
        coupling: f64,
        /// Non-negative shift.
        shift: f64,
    },
    /// Constant partner of a transverse-field vertex.
    SiteConstant {
        /// Site.
        site: u32,
        /// Matrix element.
        amplitude: f64,
    },
    /// `B = h sigma_x`.
    SpinFlip {
        /// Site.
        site: u32,
        /// Matrix element.
        amplitude: f64,
    },
    /// `B = shift + detuning*n_i`.
    RydbergDetuning {
        /// Site.
        site: u32,
        /// Signed detuning.
        detuning: f64,
        /// Shift.
        shift: f64,
    },
    /// `B = shift - interaction*n_i*n_j`.
    RydbergInteraction {
        /// First endpoint.
        site_i: u32,
        /// Second endpoint.
        site_j: u32,
        /// Interaction.
        interaction: f64,
        /// Shift.
        shift: f64,
    },
}
impl SseTerm {
    /// Return the vertex classification.
    pub fn operator_kind(self) -> OperatorKind {
        match self {
            Self::SpinFlip { .. } => OperatorKind::OffDiagonal,
            _ => OperatorKind::Diagonal,
        }
    }
    /// Return the referenced sites.
    pub fn sites(self) -> (u32, Option<u32>) {
        match self {
            Self::TfimBond { site_i, site_j, .. }
            | Self::RydbergInteraction { site_i, site_j, .. } => (site_i, Some(site_j)),
            Self::SiteConstant { site, .. }
            | Self::SpinFlip { site, .. }
            | Self::RydbergDetuning { site, .. } => (site, None),
        }
    }
}

/// Errors from SSE decomposition and propagation contracts.
#[derive(Clone, Debug, PartialEq)]
pub enum SseModelError {
    /// A sign-safe constructor received a negative or non-finite coupling.
    UnsupportedSign {
        /// Coefficient name.
        name: &'static str,
        /// Received value.
        value: f64,
    },
    /// A site is outside the model.
    InvalidSite {
        /// Invalid site identifier.
        site: u32,
        /// Model size.
        num_sites: usize,
    },
    /// A pair has identical endpoints.
    SelfPair,
    /// A vector has an invalid length.
    InvalidLength {
        /// Required length.
        expected: usize,
        /// Received length.
        actual: usize,
    },
    /// A term index is out of range.
    InvalidTermIndex {
        /// Invalid term index.
        term_index: usize,
        /// Number of terms.
        num_terms: usize,
    },
    /// A vertex classification disagrees with its term.
    InvalidOperatorKind,
    /// A matrix element is materially negative or zero.
    NonPositiveMatrixElement {
        /// Term that produced the value.
        term_index: usize,
        /// Matrix element.
        value: f64,
    },
    /// Propagation did not close the trace.
    TraceNotClosed,
    /// A non-finite numeric input or result occurred.
    NonFinite(&'static str),
    /// A canonical interaction channel is not supported by this decomposition.
    UnsupportedInteractionChannel,
}
impl Display for SseModelError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSign { name, value } => {
                write!(f, "{name} has unsupported sign/value {value}")
            }
            Self::InvalidSite { site, num_sites } => {
                write!(f, "site {site} is invalid for {num_sites} sites")
            }
            Self::SelfPair => f.write_str("an SSE pair cannot use the same endpoint twice"),
            Self::InvalidLength { expected, actual } => {
                write!(f, "expected length {expected}, got {actual}")
            }
            Self::InvalidTermIndex {
                term_index,
                num_terms,
            } => write!(f, "term {term_index} is invalid for {num_terms} terms"),
            Self::InvalidOperatorKind => {
                f.write_str("operator kind does not match the referenced term")
            }
            Self::NonPositiveMatrixElement { term_index, value } => write!(
                f,
                "term {term_index} has non-positive matrix element {value}"
            ),
            Self::TraceNotClosed => f.write_str("SSE operator string does not close the trace"),
            Self::NonFinite(name) => write!(f, "non-finite SSE {name}"),
            Self::UnsupportedInteractionChannel => {
                f.write_str("canonical interaction channel is not supported by TFIM SSE")
            }
        }
    }
}
impl std::error::Error for SseModelError {}

/// Interface consumed by states and samplers.
pub trait SseModel {
    /// Number of canonical sites.
    fn num_sites(&self) -> usize;
    /// Number of local terms.
    fn num_terms(&self) -> usize;
    /// Decomposition shift.
    fn energy_shift(&self) -> f64;
    /// Borrow a term.
    fn term(&self, index: usize) -> Option<&SseTerm>;
    /// Diagonal term indices.
    fn diagonal_term_indices(&self) -> &[u32];
    /// Whether every vertex supports the deterministic transverse-field Ising
    /// linked-cluster breakup implemented by the sampler.
    ///
    /// The default is conservative: models that do not opt in are rejected by
    /// the cluster update rather than sampled incorrectly.
    fn supports_tfim_cluster_update(&self) -> bool {
        false
    }
    /// Matching diagonal/off-diagonal single-site partner vertex, when one
    /// exists with an identical matrix element on the same site.
    fn transverse_partner(&self, index: usize) -> Option<u32> {
        let _ = index;
        None
    }
    /// Required kind for a term.
    fn operator_kind(&self, index: usize) -> Result<OperatorKind, SseModelError>;
    /// Evaluate a term on a canonical bit state.
    fn matrix_element(&self, index: usize, bits: &[BasisBit]) -> Result<f64, SseModelError>;
    /// Apply an off-diagonal term.
    fn apply_off_diagonal(&self, index: usize, bits: &mut [BasisBit]) -> Result<(), SseModelError>;
}

/// Sign-safe local TFIM or Rydberg decomposition.
#[derive(Clone, Debug)]
pub struct LocalSseModel {
    num_sites: usize,
    terms: Vec<SseTerm>,
    diagonal: Vec<u32>,
    transverse_partners: Vec<Option<u32>>,
    cluster_capable: bool,
    shift: f64,
}
impl LocalSseModel {
    /// Build `-J sum z_i z_j - h sum sigma_x` with explicit bonds.
    pub fn tfim(
        num_sites: usize,
        bonds: &[(u32, u32)],
        coupling: f64,
        field: f64,
    ) -> Result<Self, SseModelError> {
        nonnegative("TFIM coupling", coupling)?;
        nonnegative("TFIM field", field)?;
        let weighted_bonds = bonds
            .iter()
            .map(|&(site_i, site_j)| (site_i, site_j, coupling))
            .collect::<Vec<_>>();
        Self::tfim_weighted(num_sites, &weighted_bonds, &vec![field; num_sites])
    }
    /// Build TFIM with independently resolved bond couplings and site fields.
    pub fn tfim_weighted(
        num_sites: usize,
        bonds: &[(u32, u32, f64)],
        fields: &[f64],
    ) -> Result<Self, SseModelError> {
        if fields.len() != num_sites {
            return Err(SseModelError::InvalidLength {
                expected: num_sites,
                actual: fields.len(),
            });
        }
        let mut terms = Vec::with_capacity(bonds.len() + 2 * num_sites);
        let mut shift = 0.0;
        for &(site_i, site_j, coupling) in bonds {
            validate_pair(site_i, site_j, num_sites)?;
            nonnegative("TFIM bond coupling", coupling)?;
            terms.push(SseTerm::TfimBond {
                site_i,
                site_j,
                coupling,
                shift: 1.0,
            });
            shift += coupling;
        }
        for (site, field) in fields.iter().copied().enumerate() {
            nonnegative("TFIM site field", field)?;
            let site = site as u32;
            terms.push(SseTerm::SiteConstant {
                site,
                amplitude: field,
            });
            terms.push(SseTerm::SpinFlip {
                site,
                amplitude: field,
            });
            shift += field;
        }
        Self::new(num_sites, terms, shift)
    }
    /// Build TFIM from canonical resolved Ising interactions.
    ///
    /// This is the preferred adapter for callers that already use qslib-core's
    /// checked interaction vocabulary. The tuple constructor remains as a
    /// convenience for small scripts and legacy integrations.
    pub fn tfim_resolved(
        num_sites: usize,
        interactions: &[WeightedInteraction],
        fields: &[f64],
    ) -> Result<Self, SseModelError> {
        let bonds = interactions
            .iter()
            .map(|interaction| {
                if interaction.channel() != &InteractionChannel::IsingZZ {
                    return Err(SseModelError::UnsupportedInteractionChannel);
                }
                Ok((
                    interaction.bond().first().get(),
                    interaction.bond().second().get(),
                    interaction.coefficient(),
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::tfim_weighted(num_sites, &bonds, fields)
    }
    /// Build the canonical Rydberg decomposition with bit-one occupation.
    pub fn rydberg(
        num_sites: usize,
        detunings: &[f64],
        interactions: &[(u32, u32, f64)],
        omega: f64,
    ) -> Result<Self, SseModelError> {
        if detunings.len() != num_sites {
            return Err(SseModelError::InvalidLength {
                expected: num_sites,
                actual: detunings.len(),
            });
        }
        nonnegative("Rydberg omega", omega)?;
        if detunings.iter().any(|value| !value.is_finite()) {
            return Err(SseModelError::NonFinite("detuning"));
        }
        let mut terms = Vec::with_capacity(3 * num_sites + interactions.len());
        let mut shift = 0.0;
        for (index, detuning) in detunings.iter().copied().enumerate() {
            let onsite_shift = (-detuning).max(0.0);
            terms.push(SseTerm::SiteConstant {
                site: index as u32,
                amplitude: 0.5 * omega,
            });
            terms.push(SseTerm::SpinFlip {
                site: index as u32,
                amplitude: 0.5 * omega,
            });
            terms.push(SseTerm::RydbergDetuning {
                site: index as u32,
                detuning,
                shift: onsite_shift,
            });
            shift += 0.5 * omega + onsite_shift;
        }
        for &(site_i, site_j, interaction) in interactions {
            validate_pair(site_i, site_j, num_sites)?;
            if !interaction.is_finite() {
                return Err(SseModelError::NonFinite("interaction"));
            }
            let pair_shift = interaction.max(0.0);
            terms.push(SseTerm::RydbergInteraction {
                site_i,
                site_j,
                interaction,
                shift: pair_shift,
            });
            shift += pair_shift;
        }
        Self::new(num_sites, terms, shift)
    }
    /// Build Rydberg SSE from canonical resolved density-density interactions.
    ///
    /// Detunings remain explicit per-site coefficients while each canonical
    /// `RydbergDensityDensity` interaction supplies its own pair strength.
    pub fn rydberg_resolved(
        num_sites: usize,
        detunings: &[f64],
        interactions: &[WeightedInteraction],
        omega: f64,
    ) -> Result<Self, SseModelError> {
        let pairs = interactions
            .iter()
            .map(|interaction| {
                if interaction.channel() != &InteractionChannel::RydbergDensityDensity {
                    return Err(SseModelError::UnsupportedInteractionChannel);
                }
                Ok((
                    interaction.bond().first().get(),
                    interaction.bond().second().get(),
                    interaction.coefficient(),
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::rydberg(num_sites, detunings, &pairs, omega)
    }
    fn new(num_sites: usize, terms: Vec<SseTerm>, shift: f64) -> Result<Self, SseModelError> {
        if num_sites == 0 || terms.is_empty() {
            return Err(SseModelError::InvalidLength {
                expected: 1,
                actual: 0,
            });
        }
        let diagonal = terms
            .iter()
            .enumerate()
            .filter_map(|(index, term)| {
                (term.operator_kind() == OperatorKind::Diagonal).then_some(index as u32)
            })
            .collect();
        let transverse_partners = compute_transverse_partners(&terms);
        let cluster_capable = terms
            .iter()
            .zip(&transverse_partners)
            .all(|(term, partner)| match term {
                SseTerm::TfimBond { .. } => true,
                SseTerm::SiteConstant { .. } | SseTerm::SpinFlip { .. } => partner.is_some(),
                SseTerm::RydbergDetuning { .. } | SseTerm::RydbergInteraction { .. } => false,
            });
        Ok(Self {
            num_sites,
            terms,
            diagonal,
            transverse_partners,
            cluster_capable,
            shift,
        })
    }
    /// Borrow terms.
    pub fn terms(&self) -> &[SseTerm] {
        &self.terms
    }
    /// Return all off-diagonal term indices in stable order.
    pub fn off_diagonal_term_indices(&self) -> Vec<u32> {
        self.terms
            .iter()
            .enumerate()
            .filter_map(|(index, term)| {
                (term.operator_kind() == OperatorKind::OffDiagonal).then_some(index as u32)
            })
            .collect()
    }
    /// Return explicit bit-one occupation.
    pub fn occupation(&self, bits: &[BasisBit], site: usize) -> Result<f64, SseModelError> {
        if bits.len() != self.num_sites {
            return Err(SseModelError::InvalidLength {
                expected: self.num_sites,
                actual: bits.len(),
            });
        }
        bits.get(site)
            .map(|bit| (*bit == BasisBit::One) as u8 as f64)
            .ok_or(SseModelError::InvalidSite {
                site: site as u32,
                num_sites: self.num_sites,
            })
    }
}
fn compute_transverse_partners(terms: &[SseTerm]) -> Vec<Option<u32>> {
    let single_site = |term: &SseTerm| match *term {
        SseTerm::SiteConstant { site, amplitude } => {
            Some((site, amplitude, OperatorKind::Diagonal))
        }
        SseTerm::SpinFlip { site, amplitude } => Some((site, amplitude, OperatorKind::OffDiagonal)),
        _ => None,
    };
    let mut partners = vec![None; terms.len()];
    for (left_index, left) in terms.iter().enumerate() {
        let Some((left_site, left_amplitude, left_kind)) = single_site(left) else {
            continue;
        };
        partners[left_index] = terms.iter().enumerate().find_map(|(right_index, right)| {
            let (right_site, right_amplitude, right_kind) = single_site(right)?;
            (left_site == right_site
                && left_amplitude == right_amplitude
                && left_kind != right_kind)
                .then_some(right_index as u32)
        });
    }
    partners
}
fn nonnegative(name: &'static str, value: f64) -> Result<(), SseModelError> {
    if !value.is_finite() || value < 0.0 {
        Err(SseModelError::UnsupportedSign { name, value })
    } else {
        Ok(())
    }
}
fn validate_site(site: u32, num_sites: usize) -> Result<(), SseModelError> {
    if site as usize >= num_sites {
        Err(SseModelError::InvalidSite { site, num_sites })
    } else {
        Ok(())
    }
}
fn validate_pair(site_i: u32, site_j: u32, num_sites: usize) -> Result<(), SseModelError> {
    validate_site(site_i, num_sites)?;
    validate_site(site_j, num_sites)?;
    if site_i == site_j {
        Err(SseModelError::SelfPair)
    } else {
        Ok(())
    }
}
impl SseModel for LocalSseModel {
    fn num_sites(&self) -> usize {
        self.num_sites
    }
    fn num_terms(&self) -> usize {
        self.terms.len()
    }
    fn energy_shift(&self) -> f64 {
        self.shift
    }
    fn term(&self, index: usize) -> Option<&SseTerm> {
        self.terms.get(index)
    }
    fn diagonal_term_indices(&self) -> &[u32] {
        &self.diagonal
    }
    fn supports_tfim_cluster_update(&self) -> bool {
        self.cluster_capable
    }
    fn transverse_partner(&self, index: usize) -> Option<u32> {
        self.transverse_partners.get(index).copied().flatten()
    }
    fn operator_kind(&self, index: usize) -> Result<OperatorKind, SseModelError> {
        self.terms
            .get(index)
            .map(|term| term.operator_kind())
            .ok_or(SseModelError::InvalidTermIndex {
                term_index: index,
                num_terms: self.terms.len(),
            })
    }
    fn matrix_element(&self, index: usize, bits: &[BasisBit]) -> Result<f64, SseModelError> {
        if bits.len() != self.num_sites {
            return Err(SseModelError::InvalidLength {
                expected: self.num_sites,
                actual: bits.len(),
            });
        }
        let term = *self
            .terms
            .get(index)
            .ok_or(SseModelError::InvalidTermIndex {
                term_index: index,
                num_terms: self.terms.len(),
            })?;
        let z = |site: u32| bits[site as usize].pauli_eigenvalue() as f64;
        let n = |site: u32| (bits[site as usize] == BasisBit::One) as u8 as f64;
        let value = match term {
            SseTerm::TfimBond {
                site_i,
                site_j,
                coupling,
                shift,
            } => coupling * (shift + z(site_i) * z(site_j)),
            SseTerm::SiteConstant { amplitude, .. } | SseTerm::SpinFlip { amplitude, .. } => {
                amplitude
            }
            SseTerm::RydbergDetuning {
                site,
                detuning,
                shift,
            } => shift + detuning * n(site),
            SseTerm::RydbergInteraction {
                site_i,
                site_j,
                interaction,
                shift,
            } => shift - interaction * n(site_i) * n(site_j),
        };
        if !value.is_finite() || value < -1.0e-12 {
            return Err(SseModelError::NonPositiveMatrixElement {
                term_index: index,
                value,
            });
        }
        Ok(value.max(0.0))
    }
    fn apply_off_diagonal(&self, index: usize, bits: &mut [BasisBit]) -> Result<(), SseModelError> {
        if bits.len() != self.num_sites {
            return Err(SseModelError::InvalidLength {
                expected: self.num_sites,
                actual: bits.len(),
            });
        }
        let term = *self
            .terms
            .get(index)
            .ok_or(SseModelError::InvalidTermIndex {
                term_index: index,
                num_terms: self.terms.len(),
            })?;
        match term {
            SseTerm::SpinFlip { site, .. } => {
                bits[site as usize] = match bits[site as usize] {
                    BasisBit::Zero => BasisBit::One,
                    BasisBit::One => BasisBit::Zero,
                };
                Ok(())
            }
            _ => Err(SseModelError::InvalidOperatorKind),
        }
    }
}
