//! Caller-supplied TDVP statistics, QGT storage, and regularized solves.
//!
//! The module intentionally does not evaluate neural-network derivatives. A
//! caller supplies weighted local energies and complex log derivatives in
//! row-major sample order. This keeps the numerical contract usable from
//! Rust, Python, Torch, or another autodiff frontend without embedding a
//! framework in the core solver.

use blake3::Hasher;
use qslib_core::Complex64;
use std::fmt::{self, Display, Formatter};

/// Errors returned by TDVP statistics and linear algebra kernels.
#[derive(Clone, Debug, PartialEq)]
pub enum TDVPError {
    /// An input shape or scalar parameter is invalid.
    InvalidParameter(&'static str),
    /// A flat array has an unexpected length.
    Shape {
        /// Required flat length.
        expected: usize,
        /// Supplied flat length.
        actual: usize,
    },
    /// An input or derived scalar is not finite.
    NonFinite(&'static str),
    /// No positive total sample weight was supplied.
    InsufficientSamples,
    /// An iterative solve did not meet its tolerance.
    NonConvergent {
        /// Iterations attempted.
        iterations: usize,
        /// Final residual norm.
        residual: f64,
    },
}
impl Display for TDVPError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter(name) => write!(f, "invalid TDVP parameter {name}"),
            Self::Shape { expected, actual } => {
                write!(f, "TDVP expected length {expected}, received {actual}")
            }
            Self::NonFinite(context) => write!(f, "non-finite TDVP {context}"),
            Self::InsufficientSamples => f.write_str("TDVP requires positive sample weight"),
            Self::NonConvergent {
                iterations,
                residual,
            } => write!(
                f,
                "TDVP solve did not converge after {iterations} iterations (residual {residual})"
            ),
        }
    }
}
impl std::error::Error for TDVPError {}

/// Real-parameter TDVP equation-of-motion convention.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TDVPMode {
    /// Solve `S dot(theta) = Im F`.
    RealTime,
    /// Solve `S dot(theta) = -Re F`.
    ImaginaryTime,
}

/// Matrix representation used by a TDVP solve.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TDVPRepresentation {
    /// A checked dense parameter-space QGT.
    Dense,
}

/// Numerical solver used for the reported direction.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TDVPSolver {
    /// Conjugate gradients on the regularized QGT.
    ConjugateGradient,
    /// Symmetric eigensystem solve or spectral projection.
    Eigenbasis,
}

/// Deterministic parameter metadata used for checkpoint compatibility.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParameterSpec {
    name: String,
    shape: Vec<usize>,
    offset: usize,
    length: usize,
}
impl ParameterSpec {
    /// Return the stable parameter name.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Return the tensor shape in caller order.
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }
    /// Return the flat offset.
    pub fn offset(&self) -> usize {
        self.offset
    }
    /// Return the number of scalar entries.
    pub fn length(&self) -> usize {
        self.length
    }
}

/// Stable flattening order and fingerprint for real model parameters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParameterLayout {
    specs: Vec<ParameterSpec>,
    parameters: usize,
    fingerprint: String,
}
impl ParameterLayout {
    /// Construct a layout from `(name, shape)` entries in flattening order.
    pub fn new(specifications: Vec<(&str, Vec<usize>)>) -> Result<Self, TDVPError> {
        if specifications.is_empty() {
            return Err(TDVPError::InvalidParameter("parameter layout"));
        }
        let mut specs = Vec::with_capacity(specifications.len());
        let mut offset = 0usize;
        for (name, shape) in specifications {
            if name.is_empty() || specs.iter().any(|spec: &ParameterSpec| spec.name == name) {
                return Err(TDVPError::InvalidParameter("parameter names"));
            }
            let length = shape
                .iter()
                .try_fold(1usize, |product, extent| product.checked_mul(*extent))
                .ok_or(TDVPError::InvalidParameter("parameter shape"))?;
            if length == 0 {
                return Err(TDVPError::InvalidParameter("parameter shape"));
            }
            let next = offset
                .checked_add(length)
                .ok_or(TDVPError::InvalidParameter("parameter count"))?;
            specs.push(ParameterSpec {
                name: name.to_owned(),
                shape,
                offset,
                length,
            });
            offset = next;
        }
        let mut hasher = Hasher::new();
        hasher.update(b"qslib-parameter-layout-v1");
        for spec in &specs {
            hasher.update(&(spec.name.len() as u64).to_le_bytes());
            hasher.update(spec.name.as_bytes());
            hasher.update(&(spec.shape.len() as u64).to_le_bytes());
            for extent in &spec.shape {
                hasher.update(&(*extent as u64).to_le_bytes());
            }
        }
        Ok(Self {
            specs,
            parameters: offset,
            fingerprint: format!("blake3-v1-{}", hasher.finalize().to_hex()),
        })
    }
    /// Return all parameter specifications.
    pub fn specs(&self) -> &[ParameterSpec] {
        &self.specs
    }
    /// Return the flat scalar count.
    pub fn parameters(&self) -> usize {
        self.parameters
    }
    /// Return the deterministic layout fingerprint.
    pub fn fingerprint(&self) -> &str {
        &self.fingerprint
    }
}

/// Real and imaginary components of the complex TDVP force covariance.
#[derive(Clone, Debug, PartialEq)]
pub struct ForceStatistics {
    /// Real component of `F`.
    pub real: Vec<f64>,
    /// Imaginary component of `F`.
    pub imag: Vec<f64>,
}

/// Weighted sufficient statistics for real-parameter TDVP.
#[derive(Clone, Debug, PartialEq)]
pub struct TDVPStatistics {
    mode: TDVPMode,
    rhs: Vec<f64>,
    force: ForceStatistics,
    qgt: DenseQgt,
    energy_mean: Complex64,
    energy_variance: f64,
    samples: usize,
}
impl TDVPStatistics {
    /// Return the equation-of-motion mode used to construct `rhs`.
    pub fn mode(&self) -> TDVPMode {
        self.mode
    }
    /// Return the mode-specific right-hand side.
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }
    /// Return the unprojected complex force covariance.
    pub fn force(&self) -> &ForceStatistics {
        &self.force
    }
    /// Return the dense real QGT.
    pub fn qgt(&self) -> &DenseQgt {
        &self.qgt
    }
    /// Return the normalized energy mean.
    pub fn energy_mean(&self) -> Complex64 {
        self.energy_mean
    }
    /// Return `E[|E_loc - E|^2]`.
    pub fn energy_variance(&self) -> f64 {
        self.energy_variance
    }
    /// Return the number of sample rows.
    pub fn samples(&self) -> usize {
        self.samples
    }
    /// Return the number of parameters.
    pub fn parameters(&self) -> usize {
        self.rhs.len()
    }
    /// Solve this statistics object with explicit regularization and clipping.
    pub fn solve(&self, options: TDVPSolveOptions) -> Result<TDVPSolveResult, TDVPError> {
        if !options.tolerance.is_finite() || options.tolerance <= 0.0 || options.max_iterations == 0
        {
            return Err(TDVPError::InvalidParameter("TDVP solve options"));
        }
        let regularization = options.regularization.clone();
        let (spectrum_values, spectrum_vectors, spectrum_iterations) =
            symmetric_eigen(self.qgt.as_slice(), self.parameters())?;
        let spectrum_scale = spectrum_values
            .iter()
            .map(|value| value.abs())
            .fold(0.0, f64::max)
            .max(f64::MIN_POSITIVE);
        if spectrum_values
            .iter()
            .any(|value| *value < -1.0e-10 * spectrum_scale)
        {
            return Err(TDVPError::InvalidParameter("negative QGT eigenvalue"));
        }
        let (mut direction, selected_lambda, iterations, converged, solver) = match regularization {
            Regularization::FixedTikhonov { lambda } => {
                let cg =
                    self.qgt
                        .solve(&self.rhs, lambda, options.tolerance, options.max_iterations)?;
                if !cg.converged {
                    return Err(TDVPError::NonConvergent {
                        iterations: cg.iterations,
                        residual: cg.residual_norm,
                    });
                }
                (
                    cg.direction,
                    Some(lambda),
                    cg.iterations,
                    cg.converged,
                    TDVPSolver::ConjugateGradient,
                )
            }
            Regularization::Gcv { ref lambda_grid } => {
                let projected = matvec_transpose(&spectrum_vectors, self.parameters(), &self.rhs);
                let (lambda, _) = select_gcv_lambda(&spectrum_values, &projected, lambda_grid)?;
                let coefficients = projected
                    .iter()
                    .zip(&spectrum_values)
                    .map(|(value, eigenvalue)| value / (eigenvalue.max(0.0) + lambda))
                    .collect::<Vec<_>>();
                (
                    matvec(&spectrum_vectors, self.parameters(), &coefficients),
                    Some(lambda),
                    spectrum_iterations,
                    true,
                    TDVPSolver::Eigenbasis,
                )
            }
            Regularization::SpectralCutoff { relative_cutoff } => {
                if !relative_cutoff.is_finite() || !(0.0..=1.0).contains(&relative_cutoff) {
                    return Err(TDVPError::InvalidParameter("relative cutoff"));
                }
                let maximum = spectrum_values.iter().copied().fold(0.0, f64::max);
                let projected = matvec_transpose(&spectrum_vectors, self.parameters(), &self.rhs);
                let coefficients = projected
                    .iter()
                    .zip(&spectrum_values)
                    .map(|(value, eigenvalue)| {
                        if *eigenvalue > 0.0 && *eigenvalue >= relative_cutoff * maximum {
                            value / eigenvalue
                        } else {
                            0.0
                        }
                    })
                    .collect::<Vec<_>>();
                (
                    matvec(&spectrum_vectors, self.parameters(), &coefficients),
                    None,
                    spectrum_iterations,
                    true,
                    TDVPSolver::Eigenbasis,
                )
            }
        };
        let preclipping_direction = direction.clone();
        let unclipped_norm = l2_norm(&direction);
        if !unclipped_norm.is_finite() || direction.iter().any(|value| !value.is_finite()) {
            return Err(TDVPError::NonFinite("TDVP direction norm"));
        }
        let clipped = options
            .max_update_norm
            .is_some_and(|limit| limit.is_finite() && limit >= 0.0 && unclipped_norm > limit);
        if let Some(limit) = options.max_update_norm {
            if !limit.is_finite() || limit < 0.0 {
                return Err(TDVPError::InvalidParameter("maximum update norm"));
            }
            if unclipped_norm > limit && unclipped_norm > 0.0 {
                for value in &mut direction {
                    *value *= limit / unclipped_norm;
                }
            }
        }
        let metric_norm_squared = dot(&direction, &self.qgt.matvec(&direction)?)?;
        let linear_residual = self
            .qgt
            .matvec(&direction)?
            .iter()
            .zip(&self.rhs)
            .map(|(left, right)| left - right)
            .collect::<Vec<_>>();
        let linear_residual_squared = dot(&linear_residual, &linear_residual)?;
        let pre_solver_residual_squared =
            if let Regularization::SpectralCutoff { relative_cutoff } =
                options.regularization.clone()
            {
                let projected_direction =
                    matvec_transpose(&spectrum_vectors, self.parameters(), &preclipping_direction);
                let projected_rhs =
                    matvec_transpose(&spectrum_vectors, self.parameters(), &self.rhs);
                let maximum = spectrum_values.iter().copied().fold(0.0, f64::max);
                let retained = spectrum_values
                    .iter()
                    .enumerate()
                    .filter(|(_, eigenvalue)| {
                        **eigenvalue > 0.0 && **eigenvalue >= relative_cutoff * maximum
                    })
                    .map(|(index, eigenvalue)| {
                        eigenvalue * projected_direction[index] - projected_rhs[index]
                    })
                    .collect::<Vec<_>>();
                dot(&retained, &retained)?
            } else {
                let mut applied = self.qgt.matvec(&preclipping_direction)?;
                if let Some(lambda) = selected_lambda {
                    for (value, direction_value) in applied.iter_mut().zip(&preclipping_direction) {
                        *value += lambda * *direction_value;
                    }
                }
                let difference = applied
                    .iter()
                    .zip(&self.rhs)
                    .map(|(left, right)| left - right)
                    .collect::<Vec<_>>();
                dot(&difference, &difference)?
            };
        let force_projection = dot(&direction, &self.rhs)?;
        let (residual_squared, residual_floor_applied, metric_floor_applied) =
            projected_residual_diagnostics(
                self.energy_variance,
                metric_norm_squared,
                force_projection,
            )?;
        let normalized = if self.energy_variance > 0.0 {
            residual_squared / self.energy_variance
        } else {
            residual_squared
        };
        Ok(TDVPSolveResult {
            direction,
            unclipped_norm,
            clipped,
            selected_lambda,
            regularization: options.regularization,
            representation: TDVPRepresentation::Dense,
            solver,
            tolerance: options.tolerance,
            iterations,
            converged,
            linear_residual_squared,
            pre_solver_residual_squared,
            metric_norm: metric_norm_squared.max(0.0).sqrt(),
            residual_squared,
            normalized_residual: normalized,
            zero_variance_policy: "normalized residual uses unnormalized r2 when Var(H)=0",
            residual_floor_applied,
            metric_floor_applied,
            mode: self.mode,
        })
    }
}

/// Estimate weighted QGT, force, energy, and variance statistics.
pub fn estimate_tdvp(
    weights: &[f64],
    local_energies: &[Complex64],
    derivatives: &[Complex64],
    parameter_count: usize,
    mode: TDVPMode,
) -> Result<TDVPStatistics, TDVPError> {
    if weights.is_empty() || parameter_count == 0 {
        return Err(TDVPError::InsufficientSamples);
    }
    if weights.len() != local_energies.len() {
        return Err(TDVPError::Shape {
            expected: weights.len(),
            actual: local_energies.len(),
        });
    }
    let expected = weights
        .len()
        .checked_mul(parameter_count)
        .ok_or(TDVPError::InvalidParameter("derivative shape"))?;
    if derivatives.len() != expected {
        return Err(TDVPError::Shape {
            expected,
            actual: derivatives.len(),
        });
    }
    let mut total_weight = 0.0;
    for (weight, energy) in weights.iter().zip(local_energies) {
        if !weight.is_finite() || *weight < 0.0 {
            return Err(TDVPError::InvalidParameter("sample weight"));
        }
        if !energy.re.is_finite() || !energy.im.is_finite() {
            return Err(TDVPError::NonFinite("local energy"));
        }
        total_weight += weight;
    }
    if !total_weight.is_finite() || total_weight <= 0.0 {
        return Err(TDVPError::InsufficientSamples);
    }
    let normalized = weights
        .iter()
        .map(|weight| weight / total_weight)
        .collect::<Vec<_>>();
    let energy_mean = local_energies
        .iter()
        .zip(&normalized)
        .map(|(energy, weight)| *energy * *weight)
        .sum::<Complex64>();
    let mut derivative_mean = vec![Complex64::new(0.0, 0.0); parameter_count];
    for (sample, weight) in normalized.iter().enumerate() {
        for parameter in 0..parameter_count {
            derivative_mean[parameter] +=
                derivatives[sample * parameter_count + parameter] * *weight;
        }
    }
    let mut qgt = vec![0.0; parameter_count * parameter_count];
    let mut force = ForceStatistics {
        real: vec![0.0; parameter_count],
        imag: vec![0.0; parameter_count],
    };
    let mut energy_variance = 0.0;
    for (sample, weight) in normalized.iter().enumerate() {
        let centered_energy = local_energies[sample] - energy_mean;
        energy_variance += *weight * centered_energy.norm_sqr();
        for left in 0..parameter_count {
            let centered_left =
                derivatives[sample * parameter_count + left] - derivative_mean[left];
            let force_value = centered_left.conj() * centered_energy * *weight;
            force.real[left] += force_value.re;
            force.imag[left] += force_value.im;
            for right in 0..parameter_count {
                let centered_right =
                    derivatives[sample * parameter_count + right] - derivative_mean[right];
                qgt[left * parameter_count + right] +=
                    *weight * (centered_left.conj() * centered_right).re;
            }
        }
    }
    let rhs = match mode {
        TDVPMode::RealTime => force.imag.clone(),
        TDVPMode::ImaginaryTime => force.real.iter().map(|value| -*value).collect(),
    };
    if !energy_variance.is_finite()
        || qgt.iter().any(|value| !value.is_finite())
        || rhs.iter().any(|value| !value.is_finite())
    {
        return Err(TDVPError::NonFinite("derived statistics"));
    }
    Ok(TDVPStatistics {
        mode,
        rhs,
        force,
        qgt: DenseQgt {
            parameters: parameter_count,
            data: qgt,
        },
        energy_mean,
        energy_variance,
        samples: weights.len(),
    })
}

/// Estimate TDVP statistics with one unit weight per sample.
pub fn estimate_tdvp_unweighted(
    local_energies: &[Complex64],
    derivatives: &[Complex64],
    parameter_count: usize,
    mode: TDVPMode,
) -> Result<TDVPStatistics, TDVPError> {
    let weights = vec![1.0; local_energies.len()];
    estimate_tdvp(&weights, local_energies, derivatives, parameter_count, mode)
}

/// Evaluate one local energy from a diagonal value and supplied amplitude ratios.
///
/// Each pair is `(H_{b b'}, psi(b') / psi(b))` in row-oriented local-energy
/// convention. `Hamiltonian::apply` is column-oriented and returns
/// `c * P_ba`; do not conjugate that opaque coefficient wholesale for complex
/// `c`, because the row convention is `c * P_ab`. Prefer a row/local-energy
/// API, or conjugate only the Pauli matrix element when constructing the pair.
pub fn local_energy_from_ratios(
    diagonal: Complex64,
    coefficient_ratios: &[(Complex64, Complex64)],
) -> Result<Complex64, TDVPError> {
    if !diagonal.re.is_finite() || !diagonal.im.is_finite() {
        return Err(TDVPError::NonFinite("local-energy diagonal"));
    }
    let mut energy = diagonal;
    for (coefficient, ratio) in coefficient_ratios {
        if !coefficient.re.is_finite()
            || !coefficient.im.is_finite()
            || !ratio.re.is_finite()
            || !ratio.im.is_finite()
        {
            return Err(TDVPError::NonFinite("local-energy ratio"));
        }
        energy += *coefficient * *ratio;
    }
    if !energy.re.is_finite() || !energy.im.is_finite() {
        return Err(TDVPError::NonFinite("local energy"));
    }
    Ok(energy)
}

/// One chunk in a streamed QGT matrix-vector product.
#[derive(Clone, Copy, Debug)]
pub struct QgtSampleChunk<'a> {
    /// Non-negative weights, one per derivative row.
    pub weights: &'a [f64],
    /// Complex derivative rows in row-major sample order.
    pub derivatives: &'a [Complex64],
}

/// Compute a QGT matrix-vector product from streamed derivative chunks.
///
/// `mean_derivatives` must be the globally weighted derivative mean. Chunks
/// may be produced by a file reader, distributed reducer, or autodiff frontend
/// and are never assembled into a sample-space identity matrix.
pub fn qgt_vector_product_stream<'a, I>(
    chunks: I,
    mean_derivatives: &[Complex64],
    parameter_count: usize,
    vector: &[f64],
) -> Result<Vec<f64>, TDVPError>
where
    I: IntoIterator<Item = QgtSampleChunk<'a>>,
{
    if parameter_count == 0
        || mean_derivatives.len() != parameter_count
        || vector.len() != parameter_count
    {
        return Err(TDVPError::Shape {
            expected: parameter_count,
            actual: vector.len(),
        });
    }
    if mean_derivatives
        .iter()
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
        || vector.iter().any(|value| !value.is_finite())
    {
        return Err(TDVPError::NonFinite("streamed QGT input"));
    }
    let mut result = vec![0.0; parameter_count];
    let mut total_weight = 0.0;
    for chunk in chunks {
        let expected = chunk
            .weights
            .len()
            .checked_mul(parameter_count)
            .ok_or(TDVPError::InvalidParameter("streamed derivative shape"))?;
        if chunk.derivatives.len() != expected {
            return Err(TDVPError::Shape {
                expected,
                actual: chunk.derivatives.len(),
            });
        }
        for (sample, weight) in chunk.weights.iter().enumerate() {
            if !weight.is_finite() || *weight < 0.0 {
                return Err(TDVPError::InvalidParameter("streamed sample weight"));
            }
            total_weight += *weight;
            let mut projection = Complex64::new(0.0, 0.0);
            for parameter in 0..parameter_count {
                let centered = chunk.derivatives[sample * parameter_count + parameter]
                    - mean_derivatives[parameter];
                projection += centered * vector[parameter];
            }
            for parameter in 0..parameter_count {
                let centered = chunk.derivatives[sample * parameter_count + parameter]
                    - mean_derivatives[parameter];
                result[parameter] += *weight * (centered.conj() * projection).re;
            }
        }
    }
    if !total_weight.is_finite() || total_weight <= 0.0 {
        return Err(TDVPError::InsufficientSamples);
    }
    for value in &mut result {
        *value /= total_weight;
    }
    if result.iter().any(|value| !value.is_finite()) {
        return Err(TDVPError::NonFinite("streamed QGT result"));
    }
    Ok(result)
}

/// A symmetric row-major QGT matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct DenseQgt {
    parameters: usize,
    data: Vec<f64>,
}

/// Checked eigensystem metadata for a dense QGT.
#[derive(Clone, Debug, PartialEq)]
pub struct QgtSpectrum {
    /// Ascending eigenvalues.
    pub eigenvalues: Vec<f64>,
    /// Row-major eigenvectors with eigenvectors in columns.
    pub eigenvectors: Vec<f64>,
    /// Largest `||S v - lambda v||` residual.
    pub max_residual: f64,
}
impl DenseQgt {
    /// Construct a checked square matrix.
    pub fn new(parameters: usize, data: Vec<f64>) -> Result<Self, TDVPError> {
        if parameters == 0 {
            return Err(TDVPError::InvalidParameter("QGT dimension"));
        }
        let expected = parameters
            .checked_mul(parameters)
            .ok_or(TDVPError::InvalidParameter("QGT dimension"))?;
        if data.len() != expected {
            return Err(TDVPError::Shape {
                expected,
                actual: data.len(),
            });
        }
        if data.iter().any(|value| !value.is_finite()) {
            return Err(TDVPError::NonFinite("QGT"));
        }
        let scale = data.iter().map(|value| value.abs()).fold(0.0, f64::max);
        let symmetry_tolerance = 1.0e-12 * scale.max(f64::MIN_POSITIVE);
        for row in 0..parameters {
            for column in (row + 1)..parameters {
                if (data[row * parameters + column] - data[column * parameters + row]).abs()
                    > symmetry_tolerance
                {
                    return Err(TDVPError::InvalidParameter("symmetric QGT"));
                }
            }
        }
        Ok(Self { parameters, data })
    }
    /// Return the flat scalar count.
    pub fn parameters(&self) -> usize {
        self.parameters
    }
    /// Return the row-major matrix data.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }
    /// Return the checked symmetric eigensystem and residual diagnostic.
    pub fn spectrum(&self) -> Result<QgtSpectrum, TDVPError> {
        let (eigenvalues, eigenvectors, _) = symmetric_eigen(self.as_slice(), self.parameters)?;
        let mut max_residual: f64 = 0.0;
        for column in 0..self.parameters {
            let vector = (0..self.parameters)
                .map(|row| eigenvectors[row * self.parameters + column])
                .collect::<Vec<_>>();
            let applied = self.matvec(&vector)?;
            let residual = applied
                .iter()
                .zip(&vector)
                .map(|(left, right)| left - eigenvalues[column] * right)
                .fold(0.0, f64::hypot);
            max_residual = max_residual.max(residual);
        }
        if !max_residual.is_finite() {
            return Err(TDVPError::NonFinite("QGT eigen residual"));
        }
        Ok(QgtSpectrum {
            eigenvalues,
            eigenvectors,
            max_residual,
        })
    }
    /// Apply the QGT matrix to a vector.
    pub fn matvec(&self, vector: &[f64]) -> Result<Vec<f64>, TDVPError> {
        if vector.len() != self.parameters {
            return Err(TDVPError::Shape {
                expected: self.parameters,
                actual: vector.len(),
            });
        }
        if vector.iter().any(|value| !value.is_finite()) {
            return Err(TDVPError::NonFinite("QGT vector"));
        }
        let result = (0..self.parameters)
            .map(|row| {
                (0..self.parameters)
                    .map(|column| self.data[row * self.parameters + column] * vector[column])
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        if result.iter().any(|value| !value.is_finite()) {
            return Err(TDVPError::NonFinite("QGT matvec"));
        }
        Ok(result)
    }
    /// Check positive semidefiniteness through the symmetric eigenspectrum.
    pub fn is_positive_semidefinite(&self, tolerance: f64) -> bool {
        symmetric_eigen(&self.data, self.parameters)
            .map(|(values, _, _)| values.iter().all(|value| *value >= -tolerance))
            .unwrap_or(false)
    }
    /// Solve `(S + lambda I)x = rhs` by conjugate gradients.
    pub fn solve(
        &self,
        rhs: &[f64],
        lambda: f64,
        tolerance: f64,
        max_iterations: usize,
    ) -> Result<CGResult, TDVPError> {
        if !lambda.is_finite() || lambda < 0.0 {
            return Err(TDVPError::InvalidParameter("Tikhonov lambda"));
        }
        solve_cg(
            |vector| {
                let mut result = self.matvec(vector)?;
                for (value, component) in result.iter_mut().zip(vector) {
                    *value += lambda * *component;
                }
                Ok(result)
            },
            rhs,
            tolerance,
            max_iterations,
        )
    }
}

/// Conjugate-gradient termination and direction data.
#[derive(Clone, Debug, PartialEq)]
pub struct CGResult {
    direction: Vec<f64>,
    residual_norm: f64,
    iterations: usize,
    converged: bool,
}
impl CGResult {
    /// Return the solved direction.
    pub fn direction(&self) -> &[f64] {
        &self.direction
    }
    /// Return the final residual norm.
    pub fn residual_norm(&self) -> f64 {
        self.residual_norm
    }
    /// Return iterations performed.
    pub fn iterations(&self) -> usize {
        self.iterations
    }
    /// Return whether the requested tolerance was met.
    pub fn converged(&self) -> bool {
        self.converged
    }
}

/// Solve a symmetric positive-semidefinite matrix supplied as a matvec.
pub fn solve_cg<F>(
    mut matvec_fn: F,
    rhs: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Result<CGResult, TDVPError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, TDVPError>,
{
    if rhs.is_empty() || !tolerance.is_finite() || tolerance <= 0.0 || max_iterations == 0 {
        return Err(TDVPError::InvalidParameter("CG options"));
    }
    if rhs.iter().any(|value| !value.is_finite()) {
        return Err(TDVPError::NonFinite("CG right-hand side"));
    }
    let mut x = vec![0.0; rhs.len()];
    let mut residual = rhs.to_vec();
    let mut direction = residual.clone();
    let mut residual_squared = dot(&residual, &residual)?;
    let initial = residual_squared.sqrt();
    if initial == 0.0 {
        return Ok(CGResult {
            direction: x,
            residual_norm: 0.0,
            iterations: 0,
            converged: true,
        });
    }
    let target = tolerance * initial;
    for iteration in 0..max_iterations {
        let applied = matvec_fn(&direction)?;
        let denominator = dot(&direction, &applied)?;
        if !denominator.is_finite() || denominator <= 0.0 {
            return Ok(CGResult {
                direction: x,
                residual_norm: residual_squared.sqrt(),
                iterations: iteration,
                converged: false,
            });
        }
        let alpha = residual_squared / denominator;
        for index in 0..x.len() {
            x[index] += alpha * direction[index];
            residual[index] -= alpha * applied[index];
        }
        let next_squared = dot(&residual, &residual)?;
        if next_squared.sqrt() <= target {
            return Ok(CGResult {
                direction: x,
                residual_norm: next_squared.sqrt(),
                iterations: iteration + 1,
                converged: true,
            });
        }
        let beta = next_squared / residual_squared;
        for index in 0..direction.len() {
            direction[index] = residual[index] + beta * direction[index];
        }
        residual_squared = next_squared;
    }
    Ok(CGResult {
        direction: x,
        residual_norm: residual_squared.sqrt(),
        iterations: max_iterations,
        converged: false,
    })
}

/// Regularization choices for a TDVP direction.
#[derive(Clone, Debug, PartialEq)]
pub enum Regularization {
    /// Add a fixed non-negative Tikhonov shift.
    FixedTikhonov {
        /// Non-negative diagonal shift.
        lambda: f64,
    },
    /// Select a positive shift by generalized cross-validation.
    Gcv {
        /// Positive candidate shifts.
        lambda_grid: Vec<f64>,
    },
    /// Discard eigenmodes below a relative QGT cutoff.
    SpectralCutoff {
        /// Fraction of the largest eigenvalue retained.
        relative_cutoff: f64,
    },
}

/// Options for a TDVP solve.
#[derive(Clone, Debug, PartialEq)]
pub struct TDVPSolveOptions {
    /// Regularization strategy.
    pub regularization: Regularization,
    /// CG stopping tolerance for fixed-shift solves.
    pub tolerance: f64,
    /// Maximum CG iterations.
    pub max_iterations: usize,
    /// Optional Euclidean update-norm cap.
    pub max_update_norm: Option<f64>,
}
impl Default for TDVPSolveOptions {
    fn default() -> Self {
        Self {
            regularization: Regularization::FixedTikhonov { lambda: 1.0e-8 },
            tolerance: 1.0e-10,
            max_iterations: 256,
            max_update_norm: None,
        }
    }
}

/// A solved direction and explicit residual/clipping diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct TDVPSolveResult {
    direction: Vec<f64>,
    unclipped_norm: f64,
    clipped: bool,
    selected_lambda: Option<f64>,
    regularization: Regularization,
    representation: TDVPRepresentation,
    solver: TDVPSolver,
    tolerance: f64,
    iterations: usize,
    converged: bool,
    linear_residual_squared: f64,
    pre_solver_residual_squared: f64,
    metric_norm: f64,
    residual_squared: f64,
    normalized_residual: f64,
    zero_variance_policy: &'static str,
    residual_floor_applied: bool,
    metric_floor_applied: bool,
    mode: TDVPMode,
}
impl TDVPSolveResult {
    /// Return the accepted direction.
    pub fn direction(&self) -> &[f64] {
        &self.direction
    }
    /// Return the norm before clipping.
    pub fn unclipped_norm(&self) -> f64 {
        self.unclipped_norm
    }
    /// Return whether clipping was applied.
    pub fn clipped(&self) -> bool {
        self.clipped
    }
    /// Return the effective Tikhonov shift, including fixed shifts and GCV selections.
    pub fn selected_lambda(&self) -> Option<f64> {
        self.selected_lambda
    }
    /// Return the regularization strategy used for the direction.
    pub fn regularization(&self) -> &Regularization {
        &self.regularization
    }
    /// Return the matrix representation used by this Rust solve surface.
    pub fn representation(&self) -> TDVPRepresentation {
        self.representation
    }
    /// Return the numerical solver used for this direction.
    pub fn solver(&self) -> TDVPSolver {
        self.solver
    }
    /// Return the solve tolerance.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
    /// Return iterations performed by the iterative backend.
    pub fn iterations(&self) -> usize {
        self.iterations
    }
    /// Return whether the numerical solve met its tolerance.
    pub fn converged(&self) -> bool {
        self.converged
    }
    /// Return the post-clipping unregularized linear residual squared.
    pub fn linear_residual_squared(&self) -> f64 {
        self.linear_residual_squared
    }
    /// Return the regularized pre-clipping solver residual squared.
    pub fn pre_solver_residual_squared(&self) -> f64 {
        self.pre_solver_residual_squared
    }
    /// Return the QGT metric norm.
    pub fn metric_norm(&self) -> f64 {
        self.metric_norm
    }
    /// Return the squared projected residual.
    pub fn residual_squared(&self) -> f64 {
        self.residual_squared
    }
    /// Return residual normalized by energy variance when available.
    pub fn normalized_residual(&self) -> f64 {
        self.normalized_residual
    }
    /// Return the zero-energy-variance normalization policy.
    pub fn zero_variance_policy(&self) -> &'static str {
        self.zero_variance_policy
    }
    /// Return the TDVP mode used for the direction.
    pub fn mode(&self) -> TDVPMode {
        self.mode
    }
    /// Return whether bounded roundoff was floored in projected residual `r2`.
    pub fn residual_floor_applied(&self) -> bool {
        self.residual_floor_applied
    }
    /// Return whether a bounded negative metric norm was floored.
    pub fn metric_floor_applied(&self) -> bool {
        self.metric_floor_applied
    }
}

/// Evaluate GCV scores for ridge shifts in an eigenbasis.
pub fn gcv_tikhonov_scores(
    eigenvalues: &[f64],
    projected_rhs: &[f64],
    lambda_grid: &[f64],
) -> Result<Vec<f64>, TDVPError> {
    if eigenvalues.is_empty() || eigenvalues.len() != projected_rhs.len() || lambda_grid.is_empty()
    {
        return Err(TDVPError::InvalidParameter("GCV inputs"));
    }
    if eigenvalues
        .iter()
        .chain(projected_rhs)
        .any(|value| !value.is_finite())
    {
        return Err(TDVPError::NonFinite("GCV input"));
    }
    let eigenvalue_scale = eigenvalues
        .iter()
        .map(|value| value.abs())
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);
    if eigenvalues
        .iter()
        .any(|value| *value < -1.0e-10 * eigenvalue_scale)
    {
        return Err(TDVPError::InvalidParameter("negative GCV eigenvalue"));
    }
    let mut scores = Vec::with_capacity(lambda_grid.len());
    for lambda in lambda_grid {
        if !lambda.is_finite() || *lambda <= 0.0 {
            return Err(TDVPError::InvalidParameter("GCV lambda"));
        }
        let fractions = eigenvalues
            .iter()
            .map(|value| lambda / (value.max(0.0) + lambda))
            .collect::<Vec<_>>();
        let numerator = fractions
            .iter()
            .zip(projected_rhs)
            .map(|(fraction, value)| (fraction * value).powi(2))
            .sum::<f64>()
            / fractions.len() as f64;
        let denominator = (fractions.iter().sum::<f64>() / fractions.len() as f64).powi(2);
        let score = numerator / denominator.max(f64::MIN_POSITIVE);
        if !score.is_finite() {
            return Err(TDVPError::NonFinite("GCV score"));
        }
        scores.push(score);
    }
    Ok(scores)
}

fn select_gcv_lambda(
    eigenvalues: &[f64],
    projected_rhs: &[f64],
    lambda_grid: &[f64],
) -> Result<(f64, f64), TDVPError> {
    let scores = gcv_tikhonov_scores(eigenvalues, projected_rhs, lambda_grid)?;
    let (index, score) = scores
        .iter()
        .enumerate()
        .min_by(|left, right| left.1.total_cmp(right.1))
        .map(|(index, score)| (index, *score))
        .ok_or(TDVPError::InvalidParameter("GCV grid"))?;
    Ok((lambda_grid[index], score))
}

/// Evaluate the convention-defined projected TDVP residual.
pub fn projected_residual_squared(
    energy_variance: f64,
    metric_norm_squared: f64,
    force_projection: f64,
) -> Result<(f64, bool), TDVPError> {
    let (value, residual_floor, metric_floor) =
        projected_residual_diagnostics(energy_variance, metric_norm_squared, force_projection)?;
    Ok((value, residual_floor || metric_floor))
}

fn projected_residual_diagnostics(
    energy_variance: f64,
    metric_norm_squared: f64,
    force_projection: f64,
) -> Result<(f64, bool, bool), TDVPError> {
    if !energy_variance.is_finite()
        || !metric_norm_squared.is_finite()
        || !force_projection.is_finite()
        || energy_variance < 0.0
    {
        return Err(TDVPError::InvalidParameter("projected residual inputs"));
    }
    let metric_floor_applied = metric_norm_squared < 0.0;
    let raw = energy_variance + metric_norm_squared.max(0.0) - 2.0 * force_projection;
    if !raw.is_finite() {
        return Err(TDVPError::NonFinite("projected residual"));
    }
    let scale = energy_variance
        .abs()
        .max(metric_norm_squared.abs())
        .max((2.0 * force_projection).abs())
        .max(f64::MIN_POSITIVE);
    let tolerance = 1.0e-12 * scale;
    if metric_norm_squared < -tolerance {
        return Err(TDVPError::InvalidParameter("negative QGT metric norm"));
    }
    if raw < -tolerance {
        return Err(TDVPError::InvalidParameter("negative projected residual"));
    }
    Ok((raw.max(0.0), raw < 0.0, metric_floor_applied))
}

fn dot(left: &[f64], right: &[f64]) -> Result<f64, TDVPError> {
    if left.len() != right.len() {
        return Err(TDVPError::Shape {
            expected: left.len(),
            actual: right.len(),
        });
    }
    let value = left
        .iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum::<f64>();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TDVPError::NonFinite("dot product"))
    }
}
fn l2_norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |norm, value| norm.hypot(*value))
}
fn matvec(matrix: &[f64], dimension: usize, vector: &[f64]) -> Vec<f64> {
    (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix[row * dimension + column] * vector[column])
                .sum()
        })
        .collect()
}
fn matvec_transpose(matrix: &[f64], dimension: usize, vector: &[f64]) -> Vec<f64> {
    (0..dimension)
        .map(|column| {
            (0..dimension)
                .map(|row| matrix[row * dimension + column] * vector[row])
                .sum()
        })
        .collect()
}

fn symmetric_eigen(
    matrix: &[f64],
    dimension: usize,
) -> Result<(Vec<f64>, Vec<f64>, usize), TDVPError> {
    let expected = dimension
        .checked_mul(dimension)
        .ok_or(TDVPError::InvalidParameter("QGT dimension"))?;
    if dimension == 0 || matrix.len() != expected {
        return Err(TDVPError::Shape {
            expected,
            actual: matrix.len(),
        });
    }
    if matrix.iter().any(|value| !value.is_finite()) {
        return Err(TDVPError::NonFinite("QGT eigenproblem"));
    }
    let scale = matrix.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let symmetry_tolerance = 1.0e-10 * scale.max(f64::MIN_POSITIVE);
    for row in 0..dimension {
        for column in (row + 1)..dimension {
            if (matrix[row * dimension + column] - matrix[column * dimension + row]).abs()
                > symmetry_tolerance
            {
                return Err(TDVPError::InvalidParameter("symmetric QGT"));
            }
        }
    }
    let mut values = matrix.to_vec();
    let mut vectors = vec![0.0; dimension * dimension];
    for index in 0..dimension {
        vectors[index * dimension + index] = 1.0;
    }
    let max_iterations = 100usize
        .checked_mul(dimension)
        .and_then(|value| value.checked_mul(dimension))
        .ok_or(TDVPError::InvalidParameter("QGT eigensolver iterations"))?;
    let convergence_tolerance = 1.0e-14 * scale.max(f64::MIN_POSITIVE);
    let mut converged = false;
    let mut rotations = 0usize;
    for _ in 0..max_iterations {
        let mut p = 0;
        let mut q = 0;
        let mut maximum = 0.0;
        for row in 0..dimension {
            for column in (row + 1)..dimension {
                if values[row * dimension + column].abs() > maximum {
                    maximum = values[row * dimension + column].abs();
                    p = row;
                    q = column;
                }
            }
        }
        if maximum < convergence_tolerance {
            converged = true;
            break;
        }
        let theta = 0.5
            * (2.0 * values[p * dimension + q])
                .atan2(values[q * dimension + q] - values[p * dimension + p]);
        let (sin, cos) = theta.sin_cos();
        rotations += 1;
        for row in 0..dimension {
            let left = values[row * dimension + p];
            let right = values[row * dimension + q];
            values[row * dimension + p] = cos * left - sin * right;
            values[row * dimension + q] = sin * left + cos * right;
        }
        for column in 0..dimension {
            let top = values[p * dimension + column];
            let bottom = values[q * dimension + column];
            values[p * dimension + column] = cos * top - sin * bottom;
            values[q * dimension + column] = sin * top + cos * bottom;
        }
        for row in 0..dimension {
            let left = vectors[row * dimension + p];
            let right = vectors[row * dimension + q];
            vectors[row * dimension + p] = cos * left - sin * right;
            vectors[row * dimension + q] = sin * left + cos * right;
        }
    }
    if !converged {
        return Err(TDVPError::InvalidParameter("QGT eigensolver convergence"));
    }
    let mut eigenvalues = (0..dimension)
        .map(|index| values[index * dimension + index])
        .collect::<Vec<_>>();
    let mut order = (0..dimension).collect::<Vec<_>>();
    order.sort_by(|left, right| eigenvalues[*left].total_cmp(&eigenvalues[*right]));
    let sorted_values = order
        .iter()
        .map(|index| eigenvalues[*index])
        .collect::<Vec<_>>();
    let mut sorted_vectors = vec![0.0; dimension * dimension];
    for (column, source) in order.iter().enumerate() {
        for row in 0..dimension {
            sorted_vectors[row * dimension + column] = vectors[row * dimension + source];
        }
    }
    eigenvalues = sorted_values;
    if eigenvalues.iter().any(|value| !value.is_finite()) {
        return Err(TDVPError::NonFinite("QGT eigenvalue"));
    }
    Ok((eigenvalues, sorted_vectors, rotations))
}
