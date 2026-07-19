//! Exact bases, matrices, eigensolvers, thermodynamics, and evolution.
//!
//! The reference backend is deliberately small-system oriented. It consumes
//! the checked basis and Hamiltonian types from `qslib-core`, preserves their
//! canonical state order, and reports residuals rather than hiding numerical
//! failure behind a best-effort result.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

use qslib_core::{
    BasisBit, BasisError, BasisState, Complex64, Hamiltonian, OperatorError, PackedState, SiteCount,
};
use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};

/// Errors returned by exact basis, matrix, and numerical operations.
#[derive(Clone, Debug, PartialEq)]
pub enum ExactError {
    /// The core basis rejected an input.
    Basis(BasisError),
    /// The core Hamiltonian rejected a state or operator action.
    Operator(OperatorError),
    /// A dimension or allocation size overflowed.
    DimensionOverflow,
    /// A matrix was not square or a vector had the wrong length.
    Shape {
        /// Required length or dimension.
        expected: usize,
        /// Supplied length or dimension.
        actual: usize,
    },
    /// A matrix failed the Hermiticity requirement.
    NonHermitian {
        /// Row containing the mismatch.
        row: usize,
        /// Column containing the mismatch.
        column: usize,
        /// Difference from the conjugate-transposed value.
        value: Complex64,
    },
    /// A Hamiltonian action leaves a requested exact sector.
    StateOutsideBasis {
        /// Connected state that is absent from the supplied basis.
        state: BasisState,
    },
    /// A matrix or vector contains a non-finite scalar.
    NonFiniteValue {
        /// Stable context for the invalid value.
        context: &'static str,
        /// Flat element index.
        index: usize,
        /// Invalid scalar.
        value: Complex64,
    },
    /// A numerical iteration did not converge.
    NonConvergent {
        /// Number of iterations attempted.
        iterations: usize,
    },
    /// A scalar parameter was invalid.
    InvalidParameter(&'static str),
}

impl Display for ExactError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Basis(error) => error.fmt(f),
            Self::Operator(error) => error.fmt(f),
            Self::DimensionOverflow => f.write_str("exact dimension overflowed"),
            Self::Shape { expected, actual } => {
                write!(f, "expected length {expected}, received {actual}")
            }
            Self::NonHermitian { row, column, value } => write!(
                f,
                "matrix entry ({row},{column})={value:?} violates Hermiticity"
            ),
            Self::StateOutsideBasis { state } => write!(
                f,
                "Hamiltonian action leaves supplied exact basis ({} sites)",
                state.len()
            ),
            Self::NonFiniteValue {
                context,
                index,
                value,
            } => write!(f, "non-finite {context} scalar at index {index}: {value:?}"),
            Self::NonConvergent { iterations } => write!(
                f,
                "Hermitian eigensolver did not converge in {iterations} iterations"
            ),
            Self::InvalidParameter(name) => write!(f, "invalid exact parameter {name}"),
        }
    }
}

impl std::error::Error for ExactError {}
impl From<BasisError> for ExactError {
    fn from(error: BasisError) -> Self {
        Self::Basis(error)
    }
}
impl From<OperatorError> for ExactError {
    fn from(error: OperatorError) -> Self {
        Self::Operator(error)
    }
}

/// An explicit, ordered exact basis with a packed-state lookup index.
#[derive(Clone, Debug, PartialEq)]
pub struct ExactBasis {
    states: Vec<BasisState>,
    lookup: HashMap<Vec<u64>, usize>,
}

impl ExactBasis {
    /// Enumerate the complete computational basis in increasing packed order.
    pub fn full(site_count: SiteCount) -> Result<Self, ExactError> {
        let states = qslib_core::FullBasis::new(site_count)?
            .map(packed_to_dense)
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_states(states)
    }

    /// Enumerate a fixed-Hamming-weight basis in increasing packed order.
    pub fn fixed_weight(site_count: SiteCount, weight: usize) -> Result<Self, ExactError> {
        let states = qslib_core::SectorBasis::new(site_count, weight)?
            .map(packed_to_dense)
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_states(states)
    }

    fn from_states(states: Vec<BasisState>) -> Result<Self, ExactError> {
        if states.is_empty() {
            return Err(ExactError::DimensionOverflow);
        }
        let mut lookup = HashMap::with_capacity(states.len());
        for (index, state) in states.iter().enumerate() {
            lookup.insert(state.pack()?.words_le().to_vec(), index);
        }
        Ok(Self { states, lookup })
    }

    /// Return states in canonical column and row order.
    pub fn states(&self) -> &[BasisState] {
        &self.states
    }
    /// Return the Hilbert-space dimension represented by this basis.
    pub fn dimension(&self) -> usize {
        self.states.len()
    }
    /// Find a state column index, if it belongs to this basis.
    pub fn index_of(&self, state: &BasisState) -> Option<usize> {
        state
            .pack()
            .ok()
            .and_then(|packed| self.lookup.get(packed.words_le()).copied())
    }
}

fn packed_to_dense(packed: PackedState) -> Result<BasisState, ExactError> {
    let mut bits = Vec::with_capacity(packed.site_count());
    for site in 0..packed.site_count() {
        bits.push(if packed.bit(site)? == BasisBit::One {
            BasisBit::One
        } else {
            BasisBit::Zero
        });
    }
    Ok(BasisState::from_bits(&bits)?)
}

/// A row-major dense complex matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct DenseMatrix {
    dimension: usize,
    data: Vec<Complex64>,
}

impl DenseMatrix {
    /// Construct a matrix from a Hamiltonian and an exact basis.
    pub fn from_hamiltonian(h: &Hamiltonian, basis: &ExactBasis) -> Result<Self, ExactError> {
        if h.site_count().get() != basis.states[0].len() {
            return Err(ExactError::Shape {
                expected: h.site_count().get(),
                actual: basis.states[0].len(),
            });
        }
        let dimension = basis.dimension();
        let length = dimension
            .checked_mul(dimension)
            .ok_or(ExactError::DimensionOverflow)?;
        let mut data = vec![Complex64::new(0.0, 0.0); length];
        for (column, state) in basis.states.iter().enumerate() {
            for (connected, coefficient) in h.apply(state)? {
                let row = basis
                    .index_of(&connected)
                    .ok_or(ExactError::StateOutsideBasis { state: connected })?;
                data[row * dimension + column] += coefficient;
            }
        }
        Ok(Self { dimension, data })
    }
    /// Construct a matrix from row-major entries.
    pub fn new(dimension: usize, data: Vec<Complex64>) -> Result<Self, ExactError> {
        let expected = dimension
            .checked_mul(dimension)
            .ok_or(ExactError::DimensionOverflow)?;
        if data.len() != expected {
            return Err(ExactError::Shape {
                expected,
                actual: data.len(),
            });
        }
        for (index, value) in data.iter().copied().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(ExactError::NonFiniteValue {
                    context: "matrix",
                    index,
                    value,
                });
            }
        }
        Ok(Self { dimension, data })
    }
    /// Return the matrix dimension.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Borrow row-major entries.
    pub fn as_slice(&self) -> &[Complex64] {
        &self.data
    }
    /// Return one matrix element.
    pub fn get(&self, row: usize, column: usize) -> Option<Complex64> {
        (row < self.dimension && column < self.dimension)
            .then(|| self.data[row * self.dimension + column])
    }
    /// Multiply by a vector using the stored row-major matrix.
    pub fn apply(&self, vector: &[Complex64]) -> Result<Vec<Complex64>, ExactError> {
        if vector.len() != self.dimension {
            return Err(ExactError::Shape {
                expected: self.dimension,
                actual: vector.len(),
            });
        }
        for (index, value) in vector.iter().copied().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(ExactError::NonFiniteValue {
                    context: "vector",
                    index,
                    value,
                });
            }
        }
        Ok((0..self.dimension)
            .map(|row| {
                (0..self.dimension)
                    .map(|column| self.data[row * self.dimension + column] * vector[column])
                    .sum()
            })
            .collect())
    }
    /// Validate Hermiticity with a quantity-specific tolerance.
    pub fn validate_hermitian(&self, tolerance: f64) -> Result<(), ExactError> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(ExactError::InvalidParameter("tolerance"));
        }
        for row in 0..self.dimension {
            for column in 0..self.dimension {
                let difference = self.data[row * self.dimension + column]
                    - self.data[column * self.dimension + row].conj();
                if !difference.re.is_finite()
                    || !difference.im.is_finite()
                    || difference.norm() > tolerance
                {
                    return Err(ExactError::NonHermitian {
                        row,
                        column,
                        value: difference,
                    });
                }
            }
        }
        Ok(())
    }
}

/// A compressed sparse row matrix with deterministic column order.
#[derive(Clone, Debug, PartialEq)]
pub struct CsrMatrix {
    dimension: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<Complex64>,
}

impl CsrMatrix {
    /// Build CSR storage from a Hamiltonian and exact basis.
    pub fn from_hamiltonian(h: &Hamiltonian, basis: &ExactBasis) -> Result<Self, ExactError> {
        if h.site_count().get() != basis.states[0].len() {
            return Err(ExactError::Shape {
                expected: h.site_count().get(),
                actual: basis.states[0].len(),
            });
        }
        let mut rows = vec![Vec::<(usize, Complex64)>::new(); basis.dimension()];
        for (column, state) in basis.states.iter().enumerate() {
            for (connected, coefficient) in h.apply(state)? {
                let row = basis
                    .index_of(&connected)
                    .ok_or(ExactError::StateOutsideBasis { state: connected })?;
                rows[row].push((column, coefficient));
            }
        }
        let mut row_offsets = vec![0];
        let mut column_indices = Vec::new();
        let mut values = Vec::new();
        for row in rows {
            for (column, value) in row {
                if value != Complex64::new(0.0, 0.0) {
                    column_indices.push(column);
                    values.push(value);
                }
            }
            row_offsets.push(values.len());
        }
        Ok(Self {
            dimension: basis.dimension(),
            row_offsets,
            column_indices,
            values,
        })
    }
    /// Multiply by a vector.
    pub fn apply(&self, vector: &[Complex64]) -> Result<Vec<Complex64>, ExactError> {
        if vector.len() != self.dimension {
            return Err(ExactError::Shape {
                expected: self.dimension,
                actual: vector.len(),
            });
        }
        for (index, value) in vector.iter().copied().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(ExactError::NonFiniteValue {
                    context: "vector",
                    index,
                    value,
                });
            }
        }
        Ok((0..self.dimension)
            .map(|row| {
                (self.row_offsets[row]..self.row_offsets[row + 1])
                    .map(|index| self.values[index] * vector[self.column_indices[index]])
                    .sum()
            })
            .collect())
    }
    /// Validate CSR Hermiticity without materializing a dense matrix.
    pub fn validate_hermitian(&self, tolerance: f64) -> Result<(), ExactError> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(ExactError::InvalidParameter("tolerance"));
        }
        let mut entries = HashMap::new();
        for row in 0..self.dimension {
            for index in self.row_offsets[row]..self.row_offsets[row + 1] {
                let value = self.values[index];
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(ExactError::NonFiniteValue {
                        context: "CSR matrix",
                        index,
                        value,
                    });
                }
                entries.insert((row, self.column_indices[index]), value);
            }
        }
        for (&(row, column), &value) in &entries {
            let conjugate = entries
                .get(&(column, row))
                .copied()
                .unwrap_or(Complex64::new(0.0, 0.0))
                .conj();
            let difference = value - conjugate;
            if !difference.re.is_finite()
                || !difference.im.is_finite()
                || difference.norm() > tolerance
            {
                return Err(ExactError::NonHermitian {
                    row,
                    column,
                    value: difference,
                });
            }
        }
        Ok(())
    }
    /// Return the matrix dimension.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Borrow CSR row offsets.
    pub fn row_offsets(&self) -> &[usize] {
        &self.row_offsets
    }
    /// Borrow CSR column indices.
    pub fn column_indices(&self) -> &[usize] {
        &self.column_indices
    }
    /// Borrow CSR values.
    pub fn values(&self) -> &[Complex64] {
        &self.values
    }
}

/// A Hermitian eigensystem with eigenvectors stored as columns.
#[derive(Clone, Debug, PartialEq)]
pub struct Eigensystem {
    values: Vec<f64>,
    vectors: Vec<Vec<Complex64>>,
    residuals: Vec<f64>,
}
impl Eigensystem {
    /// Borrow ascending eigenvalues.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    /// Borrow one normalized eigenvector by ascending eigenvalue index.
    pub fn vector(&self, index: usize) -> Option<&[Complex64]> {
        self.vectors.get(index).map(Vec::as_slice)
    }
    /// Borrow residual norms for each eigenpair.
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }
    /// Construct the orthogonal projector onto selected eigenvectors.
    pub fn projector(&self, indices: &[usize]) -> Result<DenseMatrix, ExactError> {
        let dimension = self
            .vectors
            .first()
            .map(Vec::len)
            .ok_or(ExactError::DimensionOverflow)?;
        let mut data = vec![Complex64::new(0.0, 0.0); dimension * dimension];
        for &index in indices {
            let vector = self.vector(index).ok_or(ExactError::Shape {
                expected: self.vectors.len(),
                actual: index,
            })?;
            for row in 0..dimension {
                for column in 0..dimension {
                    data[row * dimension + column] += vector[row] * vector[column].conj();
                }
            }
        }
        DenseMatrix::new(dimension, data)
    }
}

/// Diagonalize a finite Hermitian matrix with a deterministic real-Jacobi reference solver.
pub fn diagonalize_hermitian(matrix: &DenseMatrix) -> Result<Eigensystem, ExactError> {
    matrix.validate_hermitian(1.0e-10)?;
    let n = matrix.dimension;
    if n == 0 {
        return Err(ExactError::DimensionOverflow);
    }
    let real_n = n.checked_mul(2).ok_or(ExactError::DimensionOverflow)?;
    let real_dimension = real_n
        .checked_mul(real_n)
        .ok_or(ExactError::DimensionOverflow)?;
    let mut a = vec![0.0; real_dimension];
    for row in 0..n {
        for column in 0..n {
            let value = matrix.data[row * n + column];
            a[row * real_n + column] = value.re;
            a[row * real_n + column + n] = -value.im;
            a[(row + n) * real_n + column] = value.im;
            a[(row + n) * real_n + column + n] = value.re;
        }
    }
    let mut vectors = (0..real_n)
        .map(|i| {
            let mut v = vec![0.0; real_n];
            v[i] = 1.0;
            v
        })
        .collect::<Vec<_>>();
    let max_iterations = real_dimension
        .checked_mul(100)
        .ok_or(ExactError::DimensionOverflow)?;
    for iteration in 0..max_iterations {
        let (p, q, magnitude) = max_off_diagonal(&a, real_n);
        if magnitude < 1.0e-13 {
            let mut pairs = (0..real_n)
                .map(|i| (a[i * real_n + i], i))
                .collect::<Vec<_>>();
            pairs.sort_by(|left, right| left.0.total_cmp(&right.0));
            let mut values = Vec::with_capacity(n);
            let mut complex_vectors = Vec::with_capacity(n);
            let mut cursor = 0;
            while cursor < pairs.len() && values.len() < n {
                let end = (cursor + 2).min(pairs.len());
                let physical_value = pairs[cursor].0 + (pairs[end - 1].0 - pairs[cursor].0) / 2.0;
                for (_, index) in &pairs[cursor..end] {
                    let candidate = (0..n)
                        .map(|site| {
                            Complex64::new(vectors[*index][site], vectors[*index][site + n])
                        })
                        .collect::<Vec<_>>();
                    if let Some(vector) = orthonormalize(candidate, &complex_vectors) {
                        values.push(physical_value);
                        complex_vectors.push(vector);
                        break;
                    }
                }
                cursor = end;
            }
            if values.len() != n {
                return Err(ExactError::InvalidParameter("degenerate eigenbasis"));
            }
            let residuals = complex_vectors
                .iter()
                .zip(values.iter())
                .map(|(vector, value)| residual_norm(matrix, vector, *value))
                .collect::<Result<Vec<_>, _>>()?;
            if residuals.iter().any(|residual| *residual > 1.0e-8) {
                return Err(ExactError::NonConvergent {
                    iterations: iteration,
                });
            }
            return Ok(Eigensystem {
                values,
                vectors: complex_vectors,
                residuals,
            });
        }
        jacobi_rotate(&mut a, &mut vectors, real_n, p, q);
        if iteration + 1 == max_iterations {
            return Err(ExactError::NonConvergent {
                iterations: max_iterations,
            });
        }
    }
    Err(ExactError::NonConvergent {
        iterations: max_iterations,
    })
}

fn orthonormalize(
    mut candidate: Vec<Complex64>,
    prior: &[Vec<Complex64>],
) -> Option<Vec<Complex64>> {
    for vector in prior {
        let projection = dot(vector, &candidate);
        for (value, prior_value) in candidate.iter_mut().zip(vector) {
            *value -= *prior_value * projection;
        }
    }
    let length = norm(&candidate);
    if length <= 1.0e-12 {
        None
    } else {
        for value in &mut candidate {
            *value /= length;
        }
        Some(candidate)
    }
}

fn max_off_diagonal(a: &[f64], n: usize) -> (usize, usize, f64) {
    let mut best = (0, 1.min(n.saturating_sub(1)), 0.0);
    for row in 0..n {
        for column in row + 1..n {
            let value = a[row * n + column].abs();
            if value > best.2 {
                best = (row, column, value);
            }
        }
    }
    best
}
fn jacobi_rotate(a: &mut [f64], vectors: &mut [Vec<f64>], n: usize, p: usize, q: usize) {
    let app = a[p * n + p];
    let aqq = a[q * n + q];
    let apq = a[p * n + q];
    let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
    let c = phi.cos();
    let s = phi.sin();
    for k in 0..n {
        let akp = a[k * n + p];
        let akq = a[k * n + q];
        a[k * n + p] = c * akp - s * akq;
        a[k * n + q] = s * akp + c * akq;
    }
    for k in 0..n {
        let apk = a[p * n + k];
        let aqk = a[q * n + k];
        a[p * n + k] = c * apk - s * aqk;
        a[q * n + k] = s * apk + c * aqk;
    }
    for component in 0..n {
        let vkp = vectors[p][component];
        let vkq = vectors[q][component];
        vectors[p][component] = c * vkp - s * vkq;
        vectors[q][component] = s * vkp + c * vkq;
    }
}
fn residual_norm(
    matrix: &DenseMatrix,
    vector: &[Complex64],
    value: f64,
) -> Result<f64, ExactError> {
    let applied = matrix.apply(vector)?;
    Ok(applied
        .iter()
        .zip(vector)
        .map(|(left, right)| (*left - *right * value).norm_sqr())
        .sum::<f64>()
        .sqrt())
}

/// The lowest-energy eigenpair and its residual diagnostic.
#[derive(Clone, Debug, PartialEq)]
pub struct GroundState {
    energy: f64,
    vector: Vec<Complex64>,
    residual: f64,
}
impl GroundState {
    /// Select the lowest eigenvalue from a spectrum.
    pub fn from_spectrum(spectrum: &Eigensystem) -> Result<Self, ExactError> {
        let (index, energy) = spectrum
            .values
            .iter()
            .enumerate()
            .min_by(|left, right| left.1.total_cmp(right.1))
            .ok_or(ExactError::DimensionOverflow)?;
        Ok(Self {
            energy: *energy,
            vector: spectrum.vectors[index].clone(),
            residual: spectrum.residuals[index],
        })
    }
    /// Return the ground-state energy.
    pub fn energy(&self) -> f64 {
        self.energy
    }
    /// Borrow the normalized ground-state vector.
    pub fn vector(&self) -> &[Complex64] {
        &self.vector
    }
    /// Return the norm of `H|psi>-E|psi>`.
    pub fn residual(&self) -> f64 {
        self.residual
    }
}

/// Compute the lowest eigenpair of CSR storage with deterministic
/// fully-reorthogonalized Lanczos.
pub fn ground_state_sparse(
    matrix: &CsrMatrix,
    tolerance: f64,
    max_iterations: usize,
) -> Result<GroundState, ExactError> {
    matrix.validate_hermitian(1.0e-10)?;
    let mut best: Option<GroundState> = None;
    let mut last_error = None;
    for start_index in 0..matrix.dimension {
        let mut start = vec![Complex64::new(0.0, 0.0); matrix.dimension];
        start[start_index] = Complex64::new(1.0, 0.0);
        match lanczos_from_start(matrix, tolerance, max_iterations, start) {
            Ok(candidate) => {
                if best
                    .as_ref()
                    .is_none_or(|current| candidate.energy < current.energy)
                {
                    best = Some(candidate);
                }
            }
            Err(error) => last_error = Some(error),
        }
    }
    if let Some(error) = last_error {
        return Err(error);
    }
    best.ok_or(ExactError::NonConvergent {
        iterations: max_iterations,
    })
}

fn lanczos_from_start(
    matrix: &CsrMatrix,
    tolerance: f64,
    max_iterations: usize,
    start: Vec<Complex64>,
) -> Result<GroundState, ExactError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(ExactError::InvalidParameter("tolerance"));
    }
    if max_iterations == 0 {
        return Err(ExactError::InvalidParameter("max_iterations"));
    }
    let n = matrix.dimension;
    let start_norm = norm(&start);
    if start_norm == 0.0 {
        return Err(ExactError::InvalidParameter("start vector"));
    }
    let mut current = start
        .into_iter()
        .map(|value| value / start_norm)
        .collect::<Vec<_>>();
    let mut previous = vec![Complex64::new(0.0, 0.0); n];
    let mut basis = Vec::new();
    let mut alphas = Vec::new();
    let mut betas = Vec::new();
    let iterations = max_iterations.min(n);
    for step in 0..iterations {
        let mut work = matrix.apply(&current)?;
        let alpha = dot(&current, &work).re;
        alphas.push(alpha);
        for (value, vector_value) in work.iter_mut().zip(&current) {
            *value -= *vector_value * alpha;
        }
        if step > 0 {
            for (value, vector_value) in work.iter_mut().zip(&previous) {
                *value -= *vector_value * betas[step - 1];
            }
        }
        basis.push(current.clone());
        for prior in &basis {
            let projection = dot(prior, &work);
            for (value, prior_value) in work.iter_mut().zip(prior) {
                *value -= *prior_value * projection;
            }
        }
        let beta = norm(&work);
        if step + 1 < iterations {
            betas.push(beta);
        }
        previous = current;
        if beta <= tolerance || step + 1 == iterations {
            let dimension = alphas.len();
            let mut tridiagonal = vec![Complex64::new(0.0, 0.0); dimension * dimension];
            for index in 0..dimension {
                tridiagonal[index * dimension + index] = Complex64::new(alphas[index], 0.0);
                if index + 1 < dimension {
                    tridiagonal[index * dimension + index + 1] = Complex64::new(betas[index], 0.0);
                    tridiagonal[(index + 1) * dimension + index] =
                        Complex64::new(betas[index], 0.0);
                }
            }
            let tri = DenseMatrix::new(dimension, tridiagonal)?;
            let tri_spectrum = diagonalize_hermitian(&tri)?;
            let tri_ground = GroundState::from_spectrum(&tri_spectrum)?;
            let mut vector = vec![Complex64::new(0.0, 0.0); n];
            for (coefficient, lanczos_vector) in tri_ground.vector().iter().zip(&basis) {
                for (value, component) in vector.iter_mut().zip(lanczos_vector) {
                    *value += *component * *coefficient;
                }
            }
            let residual = residual_norm_csr(matrix, &vector, tri_ground.energy())?;
            if residual <= tolerance {
                return Ok(GroundState {
                    energy: tri_ground.energy(),
                    vector,
                    residual,
                });
            }
        }
        if beta == 0.0 {
            break;
        }
        current = work.into_iter().map(|value| value / beta).collect();
    }
    Err(ExactError::NonConvergent { iterations })
}

fn dot(left: &[Complex64], right: &[Complex64]) -> Complex64 {
    left.iter().zip(right).map(|(a, b)| a.conj() * *b).sum()
}

fn norm(vector: &[Complex64]) -> f64 {
    dot(vector, vector).re.sqrt()
}

fn residual_norm_csr(
    matrix: &CsrMatrix,
    vector: &[Complex64],
    value: f64,
) -> Result<f64, ExactError> {
    Ok(matrix
        .apply(vector)?
        .iter()
        .zip(vector)
        .map(|(left, right)| (*left - *right * value).norm_sqr())
        .sum::<f64>()
        .sqrt())
}

/// Exact canonical thermal sums evaluated from an eigensystem.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalSummary {
    partition_function: f64,
    log_partition_function: f64,
    partition_function_overflowed: bool,
    energy: f64,
    heat_capacity: f64,
}
impl ThermalSummary {
    /// Evaluate partition function, mean energy, and heat capacity at `beta`.
    pub fn from_spectrum(spectrum: &Eigensystem, beta: f64) -> Result<Self, ExactError> {
        if !beta.is_finite() || beta < 0.0 {
            return Err(ExactError::InvalidParameter("beta"));
        }
        if spectrum.values.is_empty() {
            return Err(ExactError::InvalidParameter("partition function"));
        }
        let minimum = spectrum
            .values
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let weights = spectrum
            .values
            .iter()
            .map(|value| (-beta * (value - minimum)).exp())
            .collect::<Vec<_>>();
        let shifted_z: f64 = weights.iter().sum();
        if !shifted_z.is_finite() || shifted_z == 0.0 {
            return Err(ExactError::InvalidParameter("partition function"));
        }
        let log_partition_function = -beta * minimum + shifted_z.ln();
        let z = log_partition_function.exp();
        let mean_delta = weights
            .iter()
            .zip(&spectrum.values)
            .map(|(weight, value)| weight * (value - minimum))
            .sum::<f64>()
            / shifted_z;
        let variance = weights
            .iter()
            .zip(&spectrum.values)
            .map(|(weight, value)| {
                let delta = value - minimum - mean_delta;
                weight * delta * delta
            })
            .sum::<f64>()
            / shifted_z;
        let energy = minimum + mean_delta;
        Ok(Self {
            partition_function: z,
            log_partition_function,
            partition_function_overflowed: !z.is_finite(),
            energy,
            heat_capacity: beta * beta * variance,
        })
    }
    /// Return `Z(beta)`.
    pub fn partition_function(&self) -> f64 {
        self.partition_function
    }
    /// Return the numerically stable natural logarithm of `Z(beta)`.
    pub fn log_partition_function(&self) -> f64 {
        self.log_partition_function
    }
    /// Return whether `Z` exceeded the finite `f64` range; use `log_partition_function` then.
    pub fn partition_function_overflowed(&self) -> bool {
        self.partition_function_overflowed
    }
    /// Return the canonical mean energy.
    pub fn energy(&self) -> f64 {
        self.energy
    }
    /// Return `beta^2 Var(H)`.
    pub fn heat_capacity(&self) -> f64 {
        self.heat_capacity
    }
}

/// Evolve a state by `exp(-i H t)` or normalized `exp(-H tau)` using spectral data.
pub fn evolve(
    matrix: &DenseMatrix,
    initial: &[Complex64],
    time: f64,
    imaginary: bool,
) -> Result<Vec<Complex64>, ExactError> {
    if !time.is_finite() || (imaginary && time < 0.0) {
        return Err(ExactError::InvalidParameter("time"));
    }
    for (index, value) in initial.iter().copied().enumerate() {
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(ExactError::NonFiniteValue {
                context: "initial vector",
                index,
                value,
            });
        }
    }
    let spectrum = diagonalize_hermitian(matrix)?;
    if initial.len() != matrix.dimension {
        return Err(ExactError::Shape {
            expected: matrix.dimension,
            actual: initial.len(),
        });
    }
    let mut result = vec![Complex64::new(0.0, 0.0); matrix.dimension];
    let energy_shift = spectrum
        .values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    for (index, vector) in spectrum.vectors.iter().enumerate() {
        let overlap: Complex64 = vector
            .iter()
            .zip(initial)
            .map(|(basis, value)| basis.conj() * *value)
            .sum();
        let factor = if imaginary {
            Complex64::new((-(spectrum.values[index] - energy_shift) * time).exp(), 0.0)
        } else {
            Complex64::from_polar(1.0, -spectrum.values[index] * time)
        };
        for (component, value) in vector.iter().enumerate() {
            result[component] += *value * overlap * factor;
        }
    }
    if imaginary {
        let norm = result
            .iter()
            .map(|value| value.norm_sqr())
            .sum::<f64>()
            .sqrt();
        if norm == 0.0 || !norm.is_finite() {
            return Err(ExactError::InvalidParameter("imaginary-time norm"));
        }
        for value in &mut result {
            *value /= norm;
        }
    }
    Ok(result)
}
