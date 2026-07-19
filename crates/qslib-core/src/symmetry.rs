use crate::{
    BasisState, Bond, Complex64, GeometryError, InteractionTable, LatticeKind, Pauli, PauliString,
    RectangularGeometry, ResolvedModel, SiteCount, SiteId, XMajorAdapter,
};
use std::collections::BTreeSet;
use std::fmt::{self, Display, Formatter};

/// Errors raised while constructing or applying a symmetry action.
#[derive(Clone, Debug, PartialEq)]
pub enum SymmetryError {
    /// A permutation does not contain exactly one source for every destination.
    InvalidPermutation {
        /// Number of declared sites.
        site_count: usize,
    },
    /// Two permutations act on different site counts.
    SiteCountMismatch {
        /// Left site count.
        left: usize,
        /// Right site count.
        right: usize,
    },
    /// A state or amplitude vector has an incompatible size.
    StateLength {
        /// Expected number of sites or amplitudes.
        expected: usize,
        /// Supplied size.
        actual: usize,
    },
    /// A generated geometric transformation is undefined for the boundaries.
    InvalidGeometryAction(String),
    /// A finite group is missing closure or identity.
    InvalidGroup(&'static str),
    /// A character has the wrong number of values or is not unitary.
    InvalidCharacter(&'static str),
    /// A mapped operator site is outside its declared system.
    OperatorSite(SiteId),
    /// A diagonal gauge is invalid.
    InvalidGauge(&'static str),
    /// A resolved Hamiltonian term is not preserved by an action.
    BrokenHamiltonian {
        /// Term index in the resolved Hamiltonian.
        term_index: usize,
        /// Mapped-term explanation.
        reason: String,
    },
}

impl Display for SymmetryError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPermutation { site_count } => {
                write!(formatter, "invalid permutation for {site_count} sites")
            }
            Self::SiteCountMismatch { left, right } => {
                write!(
                    formatter,
                    "symmetry site counts differ: {left} versus {right}"
                )
            }
            Self::StateLength { expected, actual } => {
                write!(formatter, "symmetry expects {expected}, received {actual}")
            }
            Self::InvalidGeometryAction(message) => formatter.write_str(message),
            Self::InvalidGroup(message) => formatter.write_str(message),
            Self::InvalidCharacter(message) => formatter.write_str(message),
            Self::OperatorSite(site) => {
                write!(formatter, "operator site {} is invalid", site.get())
            }
            Self::InvalidGauge(message) => formatter.write_str(message),
            Self::BrokenHamiltonian { term_index, reason } => {
                write!(
                    formatter,
                    "Hamiltonian term {term_index} breaks symmetry: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for SymmetryError {}

impl From<GeometryError> for SymmetryError {
    fn from(error: GeometryError) -> Self {
        Self::InvalidGeometryAction(error.to_string())
    }
}

/// A gather-direction site permutation storing `source_for_destination`.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Permutation {
    site_count: SiteCount,
    source_for_destination: Vec<SiteId>,
}

impl Permutation {
    /// Construct and validate a gather-direction permutation.
    pub fn new(
        site_count: SiteCount,
        source_for_destination: Vec<SiteId>,
    ) -> Result<Self, SymmetryError> {
        if source_for_destination.len() != site_count.get() {
            return Err(SymmetryError::InvalidPermutation {
                site_count: site_count.get(),
            });
        }
        let mut seen = BTreeSet::new();
        for &source in &source_for_destination {
            if site_count.validate(source).is_err() || !seen.insert(source) {
                return Err(SymmetryError::InvalidPermutation {
                    site_count: site_count.get(),
                });
            }
        }
        Ok(Self {
            site_count,
            source_for_destination,
        })
    }

    /// Construct the identity permutation.
    pub fn identity(site_count: SiteCount) -> Result<Self, SymmetryError> {
        let source_for_destination = (0..site_count.get())
            .map(|index| {
                SiteId::try_from_usize(index).map_err(|_| SymmetryError::InvalidPermutation {
                    site_count: site_count.get(),
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            site_count,
            source_for_destination,
        })
    }

    /// Return the number of sites acted on.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Return the gather source for each destination site.
    pub fn source_for_destination(&self) -> &[SiteId] {
        &self.source_for_destination
    }

    /// Return the inverse gather permutation.
    pub fn inverse(&self) -> Result<Self, SymmetryError> {
        let mut inverse = vec![SiteId::new(0); self.site_count.get()];
        for (destination, &source) in self.source_for_destination.iter().enumerate() {
            inverse[source.get() as usize] = SiteId::try_from_usize(destination).map_err(|_| {
                SymmetryError::InvalidPermutation {
                    site_count: self.site_count.get(),
                }
            })?;
        }
        Ok(Self {
            site_count: self.site_count,
            source_for_destination: inverse,
        })
    }

    /// Compose this action after `other`, preserving gather semantics.
    pub fn compose(&self, other: &Self) -> Result<Self, SymmetryError> {
        ensure_same_sites(self, other)?;
        let sources = self
            .source_for_destination
            .iter()
            .map(|&source| other.source_for_destination[source.get() as usize])
            .collect();
        Self::new(self.site_count, sources)
    }

    /// Apply the gather action `(g b)_j = b_{p_g(j)}` to a dense state.
    pub fn apply_state(&self, state: &BasisState) -> Result<BasisState, SymmetryError> {
        if state.len() != self.site_count.get() {
            return Err(SymmetryError::StateLength {
                expected: self.site_count.get(),
                actual: state.len(),
            });
        }
        let bits = self
            .source_for_destination
            .iter()
            .map(|&source| state.bits()[source.get() as usize])
            .collect::<Vec<_>>();
        BasisState::from_bits(&bits).map_err(|_| SymmetryError::StateLength {
            expected: self.site_count.get(),
            actual: bits.len(),
        })
    }

    /// Map a Pauli support under the corresponding physical site action.
    pub fn map_pauli_string(&self, operator: &PauliString) -> Result<PauliString, SymmetryError> {
        let inverse = self.inverse()?;
        let factors = operator
            .factors()
            .iter()
            .map(|&(site, pauli)| {
                self.site_count
                    .validate(site)
                    .map_err(|_| SymmetryError::OperatorSite(site))?;
                Ok((inverse.source_for_destination[site.get() as usize], pauli))
            })
            .collect::<Result<Vec<_>, SymmetryError>>()?;
        PauliString::new(factors).map_err(|_| SymmetryError::OperatorSite(SiteId::new(u32::MAX)))
    }

    /// Map a simple bond under this action.
    pub fn map_bond(&self, bond: Bond) -> Result<Bond, SymmetryError> {
        if bond.image_translation() != (0, 0) {
            return Err(SymmetryError::InvalidGeometryAction(
                "symmetry mapping of periodic-image bonds is not implicit".into(),
            ));
        }
        let inverse = self.inverse()?;
        Bond::new(
            inverse.source_for_destination[bond.first().get() as usize],
            inverse.source_for_destination[bond.second().get() as usize],
        )
        .map_err(|_| SymmetryError::InvalidGeometryAction("mapped bond is invalid".into()))
    }
}

fn ensure_same_sites(left: &Permutation, right: &Permutation) -> Result<(), SymmetryError> {
    if left.site_count != right.site_count {
        Err(SymmetryError::SiteCountMismatch {
            left: left.site_count.get(),
            right: right.site_count.get(),
        })
    } else {
        Ok(())
    }
}

/// A validated finite group of site permutations.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGroup {
    elements: Vec<Permutation>,
}

impl FiniteGroup {
    /// Construct a group, validating identity, inverse, and closure.
    pub fn new(mut elements: Vec<Permutation>) -> Result<Self, SymmetryError> {
        if elements.is_empty() {
            return Err(SymmetryError::InvalidGroup(
                "a finite group cannot be empty",
            ));
        }
        elements.sort();
        elements.dedup();
        let identity = Permutation::identity(elements[0].site_count())?;
        if !elements.contains(&identity) {
            return Err(SymmetryError::InvalidGroup("finite group lacks identity"));
        }
        for element in &elements {
            if !elements.contains(&element.inverse()?) {
                return Err(SymmetryError::InvalidGroup("finite group lacks an inverse"));
            }
            for other in &elements {
                if !elements.contains(&element.compose(other)?) {
                    return Err(SymmetryError::InvalidGroup("finite group is not closed"));
                }
            }
        }
        Ok(Self { elements })
    }

    /// Return the group order.
    pub fn order(&self) -> usize {
        self.elements.len()
    }

    /// Return the canonical sorted group elements.
    pub fn elements(&self) -> &[Permutation] {
        &self.elements
    }

    /// Return the orbit of one dense basis state without duplicates.
    pub fn orbit(&self, state: &BasisState) -> Result<Vec<BasisState>, SymmetryError> {
        let mut orbit = Vec::new();
        for element in &self.elements {
            let transformed = element.apply_state(state)?;
            if !orbit.contains(&transformed) {
                orbit.push(transformed);
            }
        }
        orbit.sort_by_key(state_key);
        Ok(orbit)
    }

    /// Return the lexicographically least packed-state orbit representative.
    pub fn canonical_representative(
        &self,
        state: &BasisState,
    ) -> Result<BasisState, SymmetryError> {
        self.orbit(state)?
            .into_iter()
            .min_by_key(state_key)
            .ok_or(SymmetryError::InvalidGroup("empty orbit"))
    }

    /// Project amplitudes using the convention `chi(g)* psi(g^-1 b)`.
    pub fn project_amplitudes(
        &self,
        amplitudes: &[Complex64],
        character: &SymmetryCharacter,
    ) -> Result<Vec<Complex64>, SymmetryError> {
        if character.values.len() != self.order() {
            return Err(SymmetryError::InvalidCharacter(
                "character length differs from group order",
            ));
        }
        if !character.validated {
            return Err(SymmetryError::InvalidCharacter(
                "character must be validated against this group",
            ));
        }
        if let Some(validated_elements) = &character.group_elements {
            if validated_elements != &self.elements {
                return Err(SymmetryError::InvalidCharacter(
                    "character was validated against a different group",
                ));
            }
        }
        let site_count = self.elements[0].site_count().get();
        let dimension =
            1usize
                .checked_shl(site_count as u32)
                .ok_or(SymmetryError::StateLength {
                    expected: usize::MAX,
                    actual: amplitudes.len(),
                })?;
        if amplitudes.len() != dimension {
            return Err(SymmetryError::StateLength {
                expected: dimension,
                actual: amplitudes.len(),
            });
        }
        let mut projected = vec![Complex64::new(0.0, 0.0); dimension];
        for (mask, projected_value) in projected.iter_mut().enumerate() {
            let state = state_from_mask(mask, site_count)?;
            for (element, character_value) in self.elements.iter().zip(&character.values) {
                let transformed = element.inverse()?.apply_state(&state)?;
                *projected_value += character_value.conj() * amplitudes[state_mask(&transformed)];
            }
            *projected_value /= self.order() as f64;
        }
        Ok(projected)
    }
}

/// A validated one-dimensional unitary character aligned with group elements.
#[derive(Clone, Debug, PartialEq)]
pub struct SymmetryCharacter {
    values: Vec<Complex64>,
    validated: bool,
    group_elements: Option<Vec<Permutation>>,
}

impl SymmetryCharacter {
    /// Construct the trivial character for a group of the given order.
    pub fn trivial(order: usize) -> Self {
        Self {
            values: vec![Complex64::new(1.0, 0.0); order],
            validated: true,
            group_elements: None,
        }
    }

    /// Construct a character after validating unit modulus values.
    pub fn new(values: Vec<Complex64>) -> Result<Self, SymmetryError> {
        if values.iter().any(|value| {
            !value.re.is_finite() || !value.im.is_finite() || (value.norm() - 1.0).abs() > 1.0e-12
        }) {
            return Err(SymmetryError::InvalidCharacter(
                "character values must have unit modulus",
            ));
        }
        Ok(Self {
            values,
            validated: false,
            group_elements: None,
        })
    }

    /// Construct a group-aligned character and validate its homomorphism law.
    pub fn new_for_group(
        group: &FiniteGroup,
        values: Vec<Complex64>,
    ) -> Result<Self, SymmetryError> {
        if values.len() != group.order() {
            return Err(SymmetryError::InvalidCharacter(
                "character length differs from group order",
            ));
        }
        let mut character = Self::new(values)?;
        character.validated = true;
        character.group_elements = Some(group.elements.clone());
        for (left_index, left) in group.elements.iter().enumerate() {
            for (right_index, right) in group.elements.iter().enumerate() {
                let product = left.compose(right)?;
                let product_index = group
                    .elements
                    .iter()
                    .position(|candidate| *candidate == product)
                    .ok_or(SymmetryError::InvalidGroup("group product is absent"))?;
                if (character.values[product_index]
                    - character.values[left_index] * character.values[right_index])
                    .norm()
                    > 1.0e-12
                {
                    return Err(SymmetryError::InvalidCharacter(
                        "character does not satisfy the group homomorphism law",
                    ));
                }
            }
        }
        Ok(character)
    }

    /// Return values in canonical group-element order.
    pub fn values(&self) -> &[Complex64] {
        &self.values
    }
}

/// Global bitwise spin inversion.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SpinInversion {
    site_count: SiteCount,
}

/// A diagonal gauge with `phi(b) = pi * sum_i a_i b_i`.
#[derive(Clone, Debug, PartialEq)]
pub struct DiagonalGauge {
    coefficients: Vec<f64>,
}

impl DiagonalGauge {
    /// Construct a finite diagonal gauge from site coefficients.
    pub fn new(coefficients: Vec<f64>) -> Result<Self, SymmetryError> {
        if coefficients.is_empty() {
            return Err(SymmetryError::InvalidGauge(
                "a diagonal gauge needs at least one site",
            ));
        }
        if coefficients.iter().any(|value| !value.is_finite()) {
            return Err(SymmetryError::InvalidGauge(
                "gauge coefficients must be finite",
            ));
        }
        Ok(Self { coefficients })
    }

    /// Return the site coefficients.
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Return the phase `exp(i pi sum_i a_i b_i)` for a state.
    pub fn phase(&self, state: &BasisState) -> Result<Complex64, SymmetryError> {
        if state.len() != self.coefficients.len() {
            return Err(SymmetryError::StateLength {
                expected: self.coefficients.len(),
                actual: state.len(),
            });
        }
        let exponent = state
            .bits()
            .iter()
            .zip(&self.coefficients)
            .map(|(bit, coefficient)| {
                if bit.as_u8() == 0 {
                    0.0
                } else {
                    coefficient.rem_euclid(2.0)
                }
            })
            .sum::<f64>()
            .rem_euclid(2.0);
        if self.is_sign_only() {
            return Ok(if exponent.round() as i64 % 2 == 0 {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(-1.0, 0.0)
            });
        }
        let angle = std::f64::consts::PI * exponent;
        Ok(Complex64::new(angle.cos(), angle.sin()))
    }

    /// Apply the gauge phase to one wavefunction amplitude.
    pub fn apply_amplitude(
        &self,
        state: &BasisState,
        amplitude: Complex64,
    ) -> Result<Complex64, SymmetryError> {
        Ok(self.phase(state)? * amplitude)
    }

    /// Return whether all coefficients are integers and therefore sign-only.
    pub fn is_sign_only(&self) -> bool {
        self.coefficients
            .iter()
            .all(|coefficient| coefficient.fract() == 0.0)
    }
}

/// Construct the checkerboard sublattice sign gauge on a square geometry.
pub fn sublattice_gauge(geometry: &RectangularGeometry) -> Result<DiagonalGauge, SymmetryError> {
    if geometry.kind() != LatticeKind::Square {
        return Err(SymmetryError::InvalidGauge(
            "checkerboard sublattice gauge requires square geometry",
        ));
    }
    let coefficients = (0..geometry.site_count().get())
        .map(|site| {
            let x = site % geometry.lx();
            let y = site / geometry.lx();
            if (x + y) % 2 == 1 { 1.0 } else { 0.0 }
        })
        .collect();
    DiagonalGauge::new(coefficients)
}

impl SpinInversion {
    /// Construct global inversion on a declared site count.
    pub const fn new(site_count: SiteCount) -> Self {
        Self { site_count }
    }

    /// Apply `b_i -> 1-b_i`.
    pub fn apply(&self, state: &BasisState) -> Result<BasisState, SymmetryError> {
        if state.len() != self.site_count.get() {
            return Err(SymmetryError::StateLength {
                expected: self.site_count.get(),
                actual: state.len(),
            });
        }
        let bits = state
            .bits()
            .iter()
            .map(|bit| match bit {
                crate::BasisBit::Zero => crate::BasisBit::One,
                crate::BasisBit::One => crate::BasisBit::Zero,
            })
            .collect::<Vec<_>>();
        BasisState::from_bits(&bits).map_err(|_| SymmetryError::StateLength {
            expected: self.site_count.get(),
            actual: bits.len(),
        })
    }
}

/// Return whether a weighted interaction table is invariant under a site action.
pub fn is_interaction_symmetry(table: &InteractionTable, permutation: &Permutation) -> bool {
    if table.site_count() != permutation.site_count() {
        return false;
    }
    table.interactions().iter().all(|term| {
        let Ok(mapped_bond) = permutation.map_bond(term.bond()) else {
            return false;
        };
        table.interactions().iter().any(|candidate| {
            candidate.bond() == mapped_bond
                && candidate.channel() == term.channel()
                && candidate.identity().name() == term.identity().name()
                && candidate.coefficient() == term.coefficient()
        })
    })
}

/// Validate a full resolved model, including onsite and mapped Pauli terms.
pub fn validate_model_symmetry(
    model: &ResolvedModel,
    permutation: &Permutation,
) -> Result<(), SymmetryError> {
    if model.site_count() != permutation.site_count() {
        return Err(SymmetryError::SiteCountMismatch {
            left: model.site_count().get(),
            right: permutation.site_count().get(),
        });
    }
    for (term_index, (coefficient, operator)) in model.terms().iter().enumerate() {
        let mapped = permutation.map_pauli_string(operator)?;
        if !model
            .terms()
            .iter()
            .any(|(candidate_coefficient, candidate)| {
                *candidate_coefficient == *coefficient && *candidate == mapped
            })
        {
            return Err(SymmetryError::BrokenHamiltonian {
                term_index,
                reason: format!(
                    "mapped operator {mapped:?} with coefficient {coefficient:?} is absent"
                ),
            });
        }
    }
    Ok(())
}

/// Validate global spin inversion against a resolved model in its selected basis.
pub fn validate_spin_inversion(model: &ResolvedModel) -> Result<(), SymmetryError> {
    for (term_index, (coefficient, operator)) in model.terms().iter().enumerate() {
        let sign = operator
            .factors()
            .iter()
            .map(|(_, pauli)| match pauli {
                Pauli::X | Pauli::I => 1.0,
                Pauli::Y | Pauli::Z => -1.0,
            })
            .product::<f64>();
        let expected = *coefficient * Complex64::new(sign, 0.0);
        if !model
            .terms()
            .iter()
            .any(|(candidate_coefficient, candidate)| {
                *candidate_coefficient == expected && *candidate == *operator
            })
        {
            return Err(SymmetryError::BrokenHamiltonian {
                term_index,
                reason: format!("spin inversion changes coefficient to {expected:?}"),
            });
        }
    }
    Ok(())
}

/// Build a rectangular periodic or open translation.
pub fn translation(
    geometry: &RectangularGeometry,
    dx: isize,
    dy: isize,
) -> Result<Permutation, SymmetryError> {
    let source_dx = dx.checked_neg().ok_or_else(|| {
        SymmetryError::InvalidGeometryAction("translation shift cannot be negated".into())
    })?;
    let source_dy = dy.checked_neg().ok_or_else(|| {
        SymmetryError::InvalidGeometryAction("translation shift cannot be negated".into())
    })?;
    let mut sources = Vec::with_capacity(geometry.site_count().get());
    for destination in 0..geometry.site_count().get() {
        let destination = SiteId::try_from_usize(destination)
            .map_err(|_| SymmetryError::InvalidGeometryAction("site identifier overflow".into()))?;
        let coordinate = geometry.coordinate(destination)?;
        let x = shifted_coordinate(
            coordinate.x() as isize,
            source_dx,
            geometry.lx(),
            geometry.boundary_x(),
        )?;
        let y = shifted_coordinate(
            coordinate.y() as isize,
            source_dy,
            geometry.ly(),
            geometry.boundary_y(),
        )?;
        sources.push(geometry.site_id(x, y)?);
    }
    Permutation::new(geometry.site_count(), sources)
}

/// Build the translation group of a periodic rectangular geometry.
pub fn translation_group(geometry: &RectangularGeometry) -> Result<FiniteGroup, SymmetryError> {
    let mut elements = Vec::new();
    for dx in 0..geometry.lx() {
        for dy in 0..geometry.ly() {
            elements.push(translation(geometry, dx as isize, dy as isize)?);
        }
    }
    FiniteGroup::new(elements)
}

/// Build the explicit legacy x-major to canonical row-major site permutation.
pub fn legacy_x_major_permutation(lx: usize, ly: usize) -> Result<Permutation, SymmetryError> {
    let adapter = XMajorAdapter::new(lx, ly)?;
    let site_count = SiteCount::new(lx.checked_mul(ly).ok_or_else(|| {
        SymmetryError::InvalidGeometryAction("legacy geometry site count overflow".into())
    })?)
    .map_err(|_| SymmetryError::InvalidGeometryAction("empty legacy geometry".into()))?;
    let sources = (0..site_count.get())
        .map(|canonical| {
            let canonical = SiteId::try_from_usize(canonical).map_err(|_| {
                SymmetryError::InvalidPermutation {
                    site_count: site_count.get(),
                }
            })?;
            adapter
                .from_canonical(canonical)
                .map_err(SymmetryError::from)
        })
        .collect::<Result<Vec<_>, SymmetryError>>()?;
    Permutation::new(site_count, sources)
}

/// Build the four-element rectangle point group.
pub fn rectangle_point_group(geometry: &RectangularGeometry) -> Result<FiniteGroup, SymmetryError> {
    if geometry.kind() != LatticeKind::Square {
        return Err(SymmetryError::InvalidGeometryAction(
            "rectangle point group requires square rectangular embedding".into(),
        ));
    }
    (0..4)
        .map(|kind| coordinate_permutation(geometry, kind, false))
        .collect::<Result<Vec<_>, _>>()
        .and_then(FiniteGroup::new)
}

/// Build the eight-element square dihedral point group.
pub fn square_point_group(geometry: &RectangularGeometry) -> Result<FiniteGroup, SymmetryError> {
    if geometry.lx() != geometry.ly()
        || geometry.kind() != LatticeKind::Square
        || geometry.boundary_x() != geometry.boundary_y()
    {
        return Err(SymmetryError::InvalidGeometryAction(
            "square point group requires equal extents".into(),
        ));
    }
    (0..8)
        .map(|kind| coordinate_permutation(geometry, kind, true))
        .collect::<Result<Vec<_>, _>>()
        .and_then(FiniteGroup::new)
}

fn shifted_coordinate(
    value: isize,
    shift: isize,
    extent: usize,
    boundary: crate::Boundary,
) -> Result<usize, SymmetryError> {
    if boundary == crate::Boundary::Periodic {
        let extent = isize::try_from(extent)
            .map_err(|_| SymmetryError::InvalidGeometryAction("extent exceeds isize".into()))?;
        let reduced = shift.checked_rem_euclid(extent).ok_or_else(|| {
            SymmetryError::InvalidGeometryAction("translation shift cannot be reduced".into())
        })?;
        let candidate = value.checked_add(reduced).ok_or_else(|| {
            SymmetryError::InvalidGeometryAction("translation overflows coordinate".into())
        })?;
        Ok(candidate.rem_euclid(extent) as usize)
    } else {
        let candidate = value.checked_add(shift).ok_or_else(|| {
            SymmetryError::InvalidGeometryAction(
                "open-boundary translation overflows coordinate".into(),
            )
        })?;
        if (0..extent as isize).contains(&candidate) {
            Ok(candidate as usize)
        } else {
            Err(SymmetryError::InvalidGeometryAction(
                "open-boundary translation leaves geometry".into(),
            ))
        }
    }
}

fn coordinate_permutation(
    geometry: &RectangularGeometry,
    kind: usize,
    square: bool,
) -> Result<Permutation, SymmetryError> {
    let mut sources = Vec::with_capacity(geometry.site_count().get());
    for destination in 0..geometry.site_count().get() {
        let destination = SiteId::try_from_usize(destination)
            .map_err(|_| SymmetryError::InvalidGeometryAction("site identifier overflow".into()))?;
        let index = destination.get() as usize;
        let x = index % geometry.lx();
        let y = index / geometry.lx();
        let (x, y) = match (square, kind) {
            (false, 0) => (x, y),
            (false, 1) => (geometry.lx() - 1 - x, y),
            (false, 2) => (x, geometry.ly() - 1 - y),
            (false, 3) => (geometry.lx() - 1 - x, geometry.ly() - 1 - y),
            (true, 0) => (x, y),
            (true, 1) => (geometry.lx() - 1 - y, x),
            (true, 2) => (geometry.lx() - 1 - x, geometry.ly() - 1 - y),
            (true, 3) => (y, geometry.ly() - 1 - x),
            (true, 4) => (geometry.lx() - 1 - x, y),
            (true, 5) => (x, geometry.ly() - 1 - y),
            (true, 6) => (y, x),
            (true, 7) => (geometry.lx() - 1 - y, geometry.ly() - 1 - x),
            _ => {
                return Err(SymmetryError::InvalidGeometryAction(
                    "unknown point-group operation".into(),
                ));
            }
        };
        sources.push(geometry.site_id(x, y)?);
    }
    Permutation::new(geometry.site_count(), sources)
}

fn state_key(state: &BasisState) -> Vec<u8> {
    state.bits().iter().rev().map(|bit| bit.as_u8()).collect()
}

fn state_mask(state: &BasisState) -> usize {
    state
        .bits()
        .iter()
        .enumerate()
        .map(|(site, bit)| (bit.as_u8() as usize) << site)
        .sum()
}

fn state_from_mask(mask: usize, site_count: usize) -> Result<BasisState, SymmetryError> {
    let bits = (0..site_count)
        .map(|site| ((mask >> site) & 1) as u8)
        .collect::<Vec<_>>();
    BasisState::from_raw_bits(&bits).map_err(|_| SymmetryError::StateLength {
        expected: site_count,
        actual: bits.len(),
    })
}
