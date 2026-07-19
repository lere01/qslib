use crate::{
    BondMultiplicity, Complex64, DenseCouplings, Hamiltonian, InteractionChannel,
    InteractionIdentity, InteractionTable, LatticeKind, OperatorError, Pauli, PauliString, Real,
    RectangularGeometry, ShellTolerance, SimulationBasis, SiteCount, SiteId, WeightedInteraction,
};
use std::fmt::{self, Display, Formatter};
use std::ops::Deref;

/// Errors raised while constructing a physical model.
#[derive(Clone, Debug, PartialEq)]
pub enum ModelError {
    /// An operator constructor rejected a term.
    Operator(OperatorError),
    /// A site-dependent parameter has the wrong length.
    FieldLength {
        /// Required number of site parameters.
        expected: usize,
        /// Supplied number of site parameters.
        actual: usize,
    },
    /// A J1 or J2 shell received the wrong number of resolved coefficients.
    ShellLength {
        /// Shell name.
        shell: &'static str,
        /// Required number of shell bonds.
        expected: usize,
        /// Supplied number of coefficients.
        actual: usize,
    },
    /// A model was requested in a basis not supported by its conversion.
    UnsupportedBasis {
        /// Model name.
        model: &'static str,
        /// Requested simulation basis.
        basis: SimulationBasis,
    },
    /// An interaction table contains a channel for another model family.
    UnexpectedChannel {
        /// Model name.
        model: &'static str,
        /// Supplied channel.
        channel: InteractionChannel,
    },
    /// A physical scalar parameter is not finite.
    NonFiniteParameter {
        /// Parameter name.
        name: &'static str,
        /// Invalid value.
        value: Real,
    },
    /// Geometry or interaction validation failed while building a model.
    InvalidInput(String),
}

/// A validated Hamiltonian together with the physical inputs that resolved it.
///
/// The operator representation is useful for matrix and local-energy work,
/// while this wrapper keeps the model family, simulation basis, weighted pair
/// identities, and site parameters available for provenance and inspection.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedModel {
    hamiltonian: Hamiltonian,
    family: &'static str,
    basis: SimulationBasis,
    interactions: Vec<WeightedInteraction>,
    specification: ModelSpecification,
}

/// The original typed physical specification retained by a resolved model.
#[derive(Clone, Debug, PartialEq)]
pub enum ModelSpecification {
    /// Inhomogeneous transverse fields for TFIM, in site order.
    Tfim {
        /// Per-site transverse fields.
        fields: Vec<Real>,
    },
    /// Isotropic pair exchange with resolved weighted interactions.
    Heisenberg,
    /// Driven Rydberg fields, in site order.
    Rydberg {
        /// Per-site drive amplitudes.
        omega: Vec<Real>,
        /// Per-site detunings.
        detuning: Vec<Real>,
    },
    /// J1 and J2 shell coefficients, in canonical shell order.
    J1J2 {
        /// Rectangular geometry including extents, boundaries, and kind.
        geometry: RectangularGeometry,
        /// Nearest-neighbour coefficients.
        j1: Vec<Real>,
        /// Diagonal next-nearest-neighbour coefficients.
        j2: Vec<Real>,
    },
}

impl ResolvedModel {
    fn new(
        hamiltonian: Hamiltonian,
        family: &'static str,
        basis: SimulationBasis,
        interactions: Vec<WeightedInteraction>,
        specification: ModelSpecification,
    ) -> Self {
        Self {
            hamiltonian,
            family,
            basis,
            interactions,
            specification,
        }
    }

    /// Return the resolved operator Hamiltonian.
    pub const fn hamiltonian(&self) -> &Hamiltonian {
        &self.hamiltonian
    }

    /// Return the stable model-family identifier.
    pub const fn family(&self) -> &'static str {
        self.family
    }

    /// Return the simulation basis used by the builder.
    pub const fn basis(&self) -> SimulationBasis {
        self.basis
    }

    /// Return all resolved weighted interactions, including zero coefficients.
    pub fn interactions(&self) -> &[WeightedInteraction] {
        &self.interactions
    }

    /// Return the original typed physical specification.
    pub const fn specification(&self) -> &ModelSpecification {
        &self.specification
    }
}

impl Deref for ResolvedModel {
    type Target = Hamiltonian;

    fn deref(&self) -> &Self::Target {
        &self.hamiltonian
    }
}

impl Display for ModelError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Operator(error) => error.fmt(formatter),
            Self::FieldLength { expected, actual } => write!(
                formatter,
                "model needs {expected} site parameters, received {actual}"
            ),
            Self::ShellLength {
                shell,
                expected,
                actual,
            } => write!(
                formatter,
                "model shell {shell} needs {expected} coefficients, received {actual}"
            ),
            Self::UnsupportedBasis { model, basis } => write!(
                formatter,
                "model {model} does not support simulation basis {basis}"
            ),
            Self::UnexpectedChannel { model, channel } => write!(
                formatter,
                "model {model} cannot consume interaction channel {channel:?}"
            ),
            Self::NonFiniteParameter { name, value } => {
                write!(formatter, "model parameter {name}={value:?} is non-finite")
            }
            Self::InvalidInput(message) => formatter.write_str(message),
        }
    }
}
impl std::error::Error for ModelError {}
impl From<OperatorError> for ModelError {
    fn from(error: OperatorError) -> Self {
        Self::Operator(error)
    }
}

/// Build the physical transverse-field Ising Hamiltonian.
pub fn tfim(
    table: &InteractionTable,
    fields: &[Real],
    basis: SimulationBasis,
) -> Result<ResolvedModel, ModelError> {
    let n = table.site_count().get();
    if fields.len() != n {
        return Err(ModelError::FieldLength {
            expected: n,
            actual: fields.len(),
        });
    }
    let (bond_pauli, field_pauli) = match basis {
        SimulationBasis::Z => (Pauli::Z, Pauli::X),
        SimulationBasis::X => (Pauli::X, Pauli::Z),
        SimulationBasis::Y => {
            return Err(ModelError::UnsupportedBasis {
                model: "tfim",
                basis,
            });
        }
    };
    let mut terms = Vec::new();
    for term in table.interactions() {
        if *term.channel() != InteractionChannel::IsingZZ {
            return Err(ModelError::UnexpectedChannel {
                model: "tfim",
                channel: term.channel().clone(),
            });
        }
        terms.push((
            Complex64::new(-term.coefficient(), 0.0),
            PauliString::new(vec![
                (term.bond().first(), bond_pauli),
                (term.bond().second(), bond_pauli),
            ])?,
        ));
    }
    for (index, &field) in fields.iter().enumerate() {
        finite_parameter("field", field)?;
        terms.push((
            Complex64::new(-field, 0.0),
            PauliString::new(vec![(
                SiteId::try_from_usize(index).map_err(|_| ModelError::FieldLength {
                    expected: n,
                    actual: fields.len(),
                })?,
                field_pauli,
            )])?,
        ));
    }
    let hamiltonian = Hamiltonian::new_hermitian(
        SiteCount::new(n).map_err(|_| ModelError::InvalidInput("empty TFIM system".into()))?,
        Complex64::new(0.0, 0.0),
        terms,
    )
    .map_err(ModelError::from)?;
    let interactions = table.interactions().to_vec();
    Ok(ResolvedModel::new(
        hamiltonian,
        "tfim",
        basis,
        interactions,
        ModelSpecification::Tfim {
            fields: fields.to_vec(),
        },
    ))
}

/// Build isotropic pair-dependent Heisenberg exchange in any Pauli-axis basis.
pub fn heisenberg(
    table: &InteractionTable,
    basis: SimulationBasis,
) -> Result<ResolvedModel, ModelError> {
    let n = table.site_count().get();
    let mut terms = Vec::new();
    for term in table.interactions() {
        if *term.channel() != InteractionChannel::HeisenbergExchange {
            return Err(ModelError::UnexpectedChannel {
                model: "heisenberg",
                channel: term.channel().clone(),
            });
        }
        let coefficient = Complex64::new(term.coefficient() / 4.0, 0.0);
        for pauli in [Pauli::X, Pauli::Y, Pauli::Z] {
            terms.push((
                coefficient,
                PauliString::new(vec![
                    (term.bond().first(), pauli),
                    (term.bond().second(), pauli),
                ])?,
            ));
        }
    }
    let _ = basis;
    let hamiltonian = Hamiltonian::new_hermitian(
        SiteCount::new(n)
            .map_err(|_| ModelError::InvalidInput("empty Heisenberg system".into()))?,
        Complex64::new(0.0, 0.0),
        terms,
    )
    .map_err(ModelError::from)?;
    Ok(ResolvedModel::new(
        hamiltonian,
        "heisenberg",
        basis,
        table.interactions().to_vec(),
        ModelSpecification::Heisenberg,
    ))
}

/// Build the canonical driven Rydberg Hamiltonian in the z simulation basis.
pub fn rydberg(
    couplings: &DenseCouplings,
    omega: &[Real],
    detuning: &[Real],
    basis: SimulationBasis,
) -> Result<ResolvedModel, ModelError> {
    if basis != SimulationBasis::Z {
        return Err(ModelError::UnsupportedBasis {
            model: "rydberg",
            basis,
        });
    }
    let n = couplings.site_count().get();
    if omega.len() != n {
        return Err(ModelError::FieldLength {
            expected: n,
            actual: omega.len(),
        });
    }
    if detuning.len() != n {
        return Err(ModelError::FieldLength {
            expected: n,
            actual: detuning.len(),
        });
    }
    let mut constant = 0.0;
    let mut terms = Vec::new();
    for index in 0..n {
        let drive = omega[index];
        let delta = detuning[index];
        finite_parameter("omega", drive)?;
        finite_parameter("detuning", delta)?;
        constant -= delta / 2.0;
        let site = SiteId::try_from_usize(index)
            .map_err(|_| ModelError::InvalidInput("site identifier overflow".into()))?;
        terms.push((
            Complex64::new(-drive / 2.0, 0.0),
            PauliString::new(vec![(site, Pauli::X)])?,
        ));
        terms.push((
            Complex64::new(delta / 2.0, 0.0),
            PauliString::new(vec![(site, Pauli::Z)])?,
        ));
    }
    let interactions = couplings
        .to_interactions(InteractionChannel::RydbergDensityDensity)
        .map_err(|error| ModelError::InvalidInput(error.to_string()))?;
    for term in &interactions {
        let v = term.coefficient();
        constant += v / 4.0;
        let i = term.bond().first();
        let j = term.bond().second();
        terms.push((
            Complex64::new(-v / 4.0, 0.0),
            PauliString::new(vec![(i, Pauli::Z)])?,
        ));
        terms.push((
            Complex64::new(-v / 4.0, 0.0),
            PauliString::new(vec![(j, Pauli::Z)])?,
        ));
        terms.push((
            Complex64::new(v / 4.0, 0.0),
            PauliString::new(vec![(i, Pauli::Z), (j, Pauli::Z)])?,
        ));
    }
    let hamiltonian = Hamiltonian::new_hermitian(
        SiteCount::new(n).map_err(|_| ModelError::InvalidInput("empty Rydberg system".into()))?,
        Complex64::new(constant, 0.0),
        terms,
    )
    .map_err(ModelError::from)?;
    Ok(ResolvedModel::new(
        hamiltonian,
        "rydberg",
        basis,
        interactions,
        ModelSpecification::Rydberg {
            omega: omega.to_vec(),
            detuning: detuning.to_vec(),
        },
    ))
}

/// Build the square-lattice homogeneous J1-J2 Heisenberg shorthand.
pub fn j1j2(
    geometry: &RectangularGeometry,
    j1: Real,
    j2: Real,
    basis: SimulationBasis,
) -> Result<ResolvedModel, ModelError> {
    finite_parameter("j1", j1)?;
    finite_parameter("j2", j2)?;
    let j1_count = geometry
        .bonds(BondMultiplicity::Simple)
        .map_err(|error| ModelError::InvalidInput(error.to_string()))?
        .len();
    let j2_count = geometry
        .pairs_at_squared_distance(2.0, ShellTolerance::Absolute(0.0))
        .map_err(|error| ModelError::InvalidInput(error.to_string()))?
        .len();
    j1j2_disordered(geometry, &vec![j1; j1_count], &vec![j2; j2_count], basis)
}

/// Build a J1-J2 Heisenberg model with one coefficient per resolved shell bond.
///
/// Coefficients are consumed in the deterministic order returned by the
/// geometry. If a periodic geometry gives the same endpoint pair in both
/// shells, the distinct `j1` and `j2` identities are retained and their
/// operator contributions add by construction.
pub fn j1j2_disordered(
    geometry: &RectangularGeometry,
    j1: &[Real],
    j2: &[Real],
    basis: SimulationBasis,
) -> Result<ResolvedModel, ModelError> {
    finite_slice("j1", j1)?;
    finite_slice("j2", j2)?;
    if geometry.kind() != LatticeKind::Square {
        return Err(ModelError::InvalidInput(
            "j1-j2 shorthand requires square rectangular geometry".into(),
        ));
    }
    let j1_bonds = geometry
        .bonds(BondMultiplicity::Simple)
        .map_err(|error| ModelError::InvalidInput(error.to_string()))?;
    let j2_bonds = geometry
        .pairs_at_squared_distance(2.0, ShellTolerance::Absolute(0.0))
        .map_err(|error| ModelError::InvalidInput(error.to_string()))?;
    if j1.len() != j1_bonds.len() {
        return Err(ModelError::ShellLength {
            shell: "j1",
            expected: j1_bonds.len(),
            actual: j1.len(),
        });
    }
    if j2.len() != j2_bonds.len() {
        return Err(ModelError::ShellLength {
            shell: "j2",
            expected: j2_bonds.len(),
            actual: j2.len(),
        });
    }
    let mut terms = Vec::with_capacity(j1.len() + j2.len());
    for (bond, coefficient) in j1_bonds.into_iter().zip(j1.iter().copied()) {
        terms.push((
            InteractionIdentity::named(bond, InteractionChannel::HeisenbergExchange, "j1")
                .map_err(|error| ModelError::InvalidInput(error.to_string()))?,
            coefficient,
        ));
    }
    for (bond, coefficient) in j2_bonds.into_iter().zip(j2.iter().copied()) {
        terms.push((
            InteractionIdentity::named(bond, InteractionChannel::HeisenbergExchange, "j2")
                .map_err(|error| ModelError::InvalidInput(error.to_string()))?,
            coefficient,
        ));
    }
    let table = InteractionTable::new_with_identities(geometry.site_count(), terms)
        .map_err(|error| ModelError::InvalidInput(error.to_string()))?;
    let hamiltonian = heisenberg(&table, basis)?.hamiltonian.clone();
    Ok(ResolvedModel::new(
        hamiltonian,
        "j1j2",
        basis,
        table.interactions().to_vec(),
        ModelSpecification::J1J2 {
            geometry: geometry.clone(),
            j1: j1.to_vec(),
            j2: j2.to_vec(),
        },
    ))
}

fn finite_parameter(name: &'static str, value: Real) -> Result<(), ModelError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ModelError::NonFiniteParameter { name, value })
    }
}

fn finite_slice(name: &'static str, values: &[Real]) -> Result<(), ModelError> {
    values
        .iter()
        .copied()
        .try_for_each(|value| finite_parameter(name, value))
}
