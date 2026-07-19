use crate::{BasisBit, BasisState, Complex64, SiteCount, SiteId};
use std::fmt::{self, Display, Formatter};

/// Errors raised while constructing or applying local operators.
#[derive(Clone, Debug, PartialEq)]
pub enum OperatorError {
    /// A Pauli support contains the same site more than once.
    DuplicateSupport {
        /// Repeated site.
        site: SiteId,
    },
    /// An operator support is outside the Hamiltonian's site count.
    SiteOutOfRange {
        /// Invalid site.
        site: SiteId,
        /// Declared site count.
        site_count: usize,
    },
    /// A coefficient contains NaN or infinity.
    NonFiniteCoefficient {
        /// Invalid coefficient.
        value: Complex64,
    },
    /// The state length does not match the Hamiltonian.
    StateLength {
        /// Required state length.
        expected: usize,
        /// Supplied state length.
        actual: usize,
    },
    /// A Hermitian Hamiltonian received a complex coefficient.
    NonHermitianCoefficient {
        /// Non-real coefficient.
        value: Complex64,
    },
    /// The reference wavefunction amplitude is absent or zero.
    ZeroReferenceAmplitude,
    /// A connected state has no supplied wavefunction amplitude.
    MissingAmplitude {
        /// Connected state with no supplied amplitude.
        state: BasisState,
    },
}

impl Display for OperatorError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateSupport { site } => {
                write!(formatter, "operator support repeats site {}", site.get())
            }
            Self::SiteOutOfRange { site, site_count } => write!(
                formatter,
                "operator site {} is outside {site_count} sites",
                site.get()
            ),
            Self::NonFiniteCoefficient { value } => {
                write!(formatter, "operator coefficient {value:?} is non-finite")
            }
            Self::StateLength { expected, actual } => write!(
                formatter,
                "operator expects {expected} sites, received {actual}"
            ),
            Self::NonHermitianCoefficient { value } => write!(
                formatter,
                "Hermitian operator coefficient {value:?} is non-real"
            ),
            Self::ZeroReferenceAmplitude => {
                formatter.write_str("local-energy reference amplitude is zero or missing")
            }
            Self::MissingAmplitude { state } => write!(
                formatter,
                "local-energy amplitude table omits a connected {}-site state",
                state.len()
            ),
        }
    }
}

impl std::error::Error for OperatorError {}

/// One local Pauli factor in a simulation basis.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum Pauli {
    /// Identity factor.
    I,
    /// Pauli X factor.
    X,
    /// Pauli Y factor.
    Y,
    /// Pauli Z factor.
    Z,
}

/// A Pauli product with distinct, canonical site support.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct PauliString {
    factors: Vec<(SiteId, Pauli)>,
}

impl PauliString {
    /// Construct a product and reject repeated support sites.
    pub fn new(mut factors: Vec<(SiteId, Pauli)>) -> Result<Self, OperatorError> {
        factors.retain(|(_, pauli)| *pauli != Pauli::I);
        factors.sort_unstable_by_key(|factor| factor.0);
        for pair in factors.windows(2) {
            if pair[0].0 == pair[1].0 {
                return Err(OperatorError::DuplicateSupport { site: pair[0].0 });
            }
        }
        Ok(Self { factors })
    }

    /// Reduce an ordered product, returning its canonical Pauli string and phase.
    pub fn product(factors: Vec<(SiteId, Pauli)>) -> Result<(Self, Complex64), OperatorError> {
        let mut reduced: Vec<(SiteId, Pauli)> = Vec::new();
        let mut phase = Complex64::new(1.0, 0.0);
        for (site, pauli) in factors {
            if pauli == Pauli::I {
                continue;
            }
            if let Some((_, existing)) =
                reduced.iter_mut().find(|(candidate, _)| *candidate == site)
            {
                let (next, factor) = multiply_pauli(*existing, pauli);
                *existing = next;
                phase *= factor;
            } else {
                reduced.push((site, pauli));
            }
        }
        Ok((Self::new(reduced)?, phase))
    }

    /// Return canonical support factors.
    pub fn factors(&self) -> &[(SiteId, Pauli)] {
        &self.factors
    }

    /// Apply the Pauli product to one dense binary state.
    pub fn apply(&self, state: &BasisState) -> Result<(BasisState, Complex64), OperatorError> {
        let mut bits = state.bits().to_vec();
        let mut coefficient = Complex64::new(1.0, 0.0);
        for &(site, pauli) in &self.factors {
            let index = site.get() as usize;
            let bit = bits.get_mut(index).ok_or(OperatorError::SiteOutOfRange {
                site,
                site_count: state.len(),
            })?;
            match (pauli, *bit) {
                (Pauli::I, _) => {}
                (Pauli::X, value) => *bit = flip(value),
                (Pauli::Y, value) => {
                    coefficient *= if value == BasisBit::Zero {
                        Complex64::new(0.0, 1.0)
                    } else {
                        Complex64::new(0.0, -1.0)
                    };
                    *bit = flip(value);
                }
                (Pauli::Z, value) => {
                    coefficient *= Complex64::new(value.pauli_eigenvalue() as f64, 0.0)
                }
            }
        }
        BasisState::from_bits(&bits)
            .map(|result| (result, coefficient))
            .map_err(|_| OperatorError::StateLength {
                expected: state.len(),
                actual: 0,
            })
    }
}

fn flip(bit: BasisBit) -> BasisBit {
    match bit {
        BasisBit::Zero => BasisBit::One,
        BasisBit::One => BasisBit::Zero,
    }
}

/// A deterministic Hamiltonian term list with a physical scalar constant.
#[derive(Clone, Debug, PartialEq)]
pub struct Hamiltonian {
    site_count: SiteCount,
    constant: Complex64,
    terms: Vec<(Complex64, PauliString)>,
}

impl Hamiltonian {
    /// Validate and sort a Hamiltonian term list.
    pub fn new(
        site_count: SiteCount,
        mut constant: Complex64,
        terms: Vec<(Complex64, PauliString)>,
    ) -> Result<Self, OperatorError> {
        if !constant.re.is_finite() || !constant.im.is_finite() {
            return Err(OperatorError::NonFiniteCoefficient { value: constant });
        }
        for (coefficient, operator) in &terms {
            if !coefficient.re.is_finite() || !coefficient.im.is_finite() {
                return Err(OperatorError::NonFiniteCoefficient {
                    value: *coefficient,
                });
            }
            for &(site, _) in operator.factors() {
                site_count
                    .validate(site)
                    .map_err(|_| OperatorError::SiteOutOfRange {
                        site,
                        site_count: site_count.get(),
                    })?;
            }
        }
        let mut constants = vec![constant];
        let mut non_identity = Vec::new();
        for (coefficient, operator) in terms {
            if operator.factors().is_empty() {
                constants.push(coefficient);
                continue;
            }
            non_identity.push((coefficient, operator));
        }
        constants.sort_by(complex_total_cmp);
        constant = deterministic_complex_sum(&constants);
        non_identity.sort_by(
            |(left_coefficient, left_operator), (right_coefficient, right_operator)| {
                left_operator
                    .cmp(right_operator)
                    .then_with(|| complex_total_cmp(left_coefficient, right_coefficient))
            },
        );
        let mut combined: Vec<(Complex64, PauliString)> = Vec::new();
        let mut pending_coefficients: Vec<Complex64> = Vec::new();
        let mut pending_operator: Option<PauliString> = None;
        for (coefficient, operator) in non_identity {
            if pending_operator.as_ref() != Some(&operator) {
                if let Some(previous_operator) = pending_operator.take() {
                    let resolved = deterministic_complex_sum(&pending_coefficients);
                    if !resolved.re.is_finite() || !resolved.im.is_finite() {
                        return Err(OperatorError::NonFiniteCoefficient { value: resolved });
                    }
                    combined.push((resolved, previous_operator));
                    pending_coefficients.clear();
                }
                pending_operator = Some(operator);
            }
            pending_coefficients.push(coefficient);
        }
        if let Some(previous_operator) = pending_operator {
            let resolved = deterministic_complex_sum(&pending_coefficients);
            if !resolved.re.is_finite() || !resolved.im.is_finite() {
                return Err(OperatorError::NonFiniteCoefficient { value: resolved });
            }
            combined.push((resolved, previous_operator));
        }
        combined.retain(|(coefficient, _)| !is_exact_zero(*coefficient));
        if !constant.re.is_finite() || !constant.im.is_finite() {
            return Err(OperatorError::NonFiniteCoefficient { value: constant });
        }
        combined.sort_by(|a, b| a.1.cmp(&b.1));
        Ok(Self {
            site_count,
            constant,
            terms: combined,
        })
    }

    /// Construct a Hamiltonian and require every term coefficient to be real.
    pub fn new_hermitian(
        site_count: SiteCount,
        constant: Complex64,
        terms: Vec<(Complex64, PauliString)>,
    ) -> Result<Self, OperatorError> {
        let result = Self::new(site_count, constant, terms)?;
        if result.constant.im != 0.0 {
            return Err(OperatorError::NonHermitianCoefficient {
                value: result.constant,
            });
        }
        if let Some((coefficient, _)) = result
            .terms
            .iter()
            .find(|(coefficient, _)| coefficient.im != 0.0)
        {
            return Err(OperatorError::NonHermitianCoefficient {
                value: *coefficient,
            });
        }
        Ok(result)
    }

    /// Number of sites represented by the Hamiltonian.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Return the physical scalar constant.
    pub const fn constant(&self) -> Complex64 {
        self.constant
    }

    /// Return deterministic coefficient/operator terms.
    pub fn terms(&self) -> &[(Complex64, PauliString)] {
        &self.terms
    }

    /// Apply the Hamiltonian and combine duplicate connected states.
    pub fn apply(&self, state: &BasisState) -> Result<Vec<(BasisState, Complex64)>, OperatorError> {
        if state.len() != self.site_count.get() {
            return Err(OperatorError::StateLength {
                expected: self.site_count.get(),
                actual: state.len(),
            });
        }
        let mut result = vec![(state.clone(), self.constant)];
        for &(coefficient, ref operator) in &self.terms {
            let (connected, matrix_element) = operator.apply(state)?;
            add_transition(&mut result, connected, coefficient * matrix_element);
        }
        result.retain(|(_, coefficient)| !is_exact_zero(*coefficient));
        result.sort_by_key(|(basis, _)| basis_key(basis));
        Ok(result)
    }

    /// Evaluate `sum_b H_ab psi_b / psi_a` from an explicit amplitude table.
    pub fn local_energy(
        &self,
        state: &BasisState,
        amplitudes: &[(BasisState, Complex64)],
    ) -> Result<Complex64, OperatorError> {
        let reference = amplitudes
            .iter()
            .find(|(candidate, _)| candidate == state)
            .map(|(_, amplitude)| *amplitude)
            .ok_or(OperatorError::ZeroReferenceAmplitude)?;
        if reference == Complex64::new(0.0, 0.0) {
            return Err(OperatorError::ZeroReferenceAmplitude);
        }
        let mut row = vec![(state.clone(), self.constant)];
        for &(coefficient, ref operator) in &self.terms {
            let (connected, matrix_element) = operator.apply(state)?;
            // Pauli strings are Hermitian, so P_ab = conjugate(P_ba),
            // while a general Hamiltonian coefficient remains c rather than
            // c*. This keeps the public complex-coefficient path correct.
            add_transition(&mut row, connected, coefficient * matrix_element.conj());
        }
        row.retain(|(_, coefficient)| !is_exact_zero(*coefficient));
        let mut numerator = Complex64::new(0.0, 0.0);
        for (connected, coefficient) in row {
            let amplitude = if connected == *state {
                reference
            } else {
                amplitudes
                    .iter()
                    .find(|(candidate, _)| *candidate == connected)
                    .map(|(_, value)| *value)
                    .ok_or_else(|| OperatorError::MissingAmplitude {
                        state: connected.clone(),
                    })?
            };
            numerator += coefficient * amplitude;
        }
        Ok(numerator / reference)
    }
}

fn complex_total_cmp(left: &Complex64, right: &Complex64) -> std::cmp::Ordering {
    left.re
        .total_cmp(&right.re)
        .then_with(|| left.im.total_cmp(&right.im))
}

fn deterministic_complex_sum(values: &[Complex64]) -> Complex64 {
    let mut sum = Complex64::new(0.0, 0.0);
    let mut correction = Complex64::new(0.0, 0.0);
    for value in values {
        let next = sum + *value;
        correction.re += if sum.re.abs() >= value.re.abs() {
            (sum.re - next.re) + value.re
        } else {
            (value.re - next.re) + sum.re
        };
        correction.im += if sum.im.abs() >= value.im.abs() {
            (sum.im - next.im) + value.im
        } else {
            (value.im - next.im) + sum.im
        };
        sum = next;
    }
    sum + correction
}

fn is_exact_zero(value: Complex64) -> bool {
    value.re == 0.0 && value.im == 0.0
}

fn add_transition(
    result: &mut Vec<(BasisState, Complex64)>,
    state: BasisState,
    coefficient: Complex64,
) {
    if let Some((_, existing)) = result.iter_mut().find(|(candidate, _)| *candidate == state) {
        *existing += coefficient;
    } else {
        result.push((state, coefficient));
    }
}

fn basis_key(state: &BasisState) -> Vec<u8> {
    state.bits().iter().rev().map(|bit| bit.as_u8()).collect()
}

fn multiply_pauli(left: Pauli, right: Pauli) -> (Pauli, Complex64) {
    let i = Complex64::new(0.0, 1.0);
    match (left, right) {
        (Pauli::I, value) | (value, Pauli::I) => (value, Complex64::new(1.0, 0.0)),
        (Pauli::X, Pauli::X) | (Pauli::Y, Pauli::Y) | (Pauli::Z, Pauli::Z) => {
            (Pauli::I, Complex64::new(1.0, 0.0))
        }
        (Pauli::X, Pauli::Y) => (Pauli::Z, i),
        (Pauli::Y, Pauli::X) => (Pauli::Z, -i),
        (Pauli::Y, Pauli::Z) => (Pauli::X, i),
        (Pauli::Z, Pauli::Y) => (Pauli::X, -i),
        (Pauli::Z, Pauli::X) => (Pauli::Y, i),
        (Pauli::X, Pauli::Z) => (Pauli::Y, -i),
    }
}
