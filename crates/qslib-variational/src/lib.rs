//! Variational statistics, TDVP linear solves, and integration kernels.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

use qslib_core::Complex64;
use std::fmt::{self, Display, Formatter};

mod evolution;
mod tdvp;

pub use evolution::*;
pub use tdvp::*;

/// Errors from weighted statistics and chain diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub enum StatisticsError {
    /// A weight, sample, or diagnostic parameter was non-finite.
    NonFinite(&'static str),
    /// A weight was negative or the total weight was zero.
    InvalidWeight,
    /// Two fixed-disorder records used the same identifier.
    DuplicateRealization,
    /// An estimator needs more samples or chains.
    InsufficientSamples,
    /// Chains have inconsistent lengths.
    ChainLengthMismatch,
}
impl Display for StatisticsError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(context) => write!(f, "non-finite {context}"),
            Self::InvalidWeight => {
                f.write_str("weights must be finite and non-negative with positive total")
            }
            Self::DuplicateRealization => f.write_str("realization identifiers must be unique"),
            Self::InsufficientSamples => {
                f.write_str("insufficient samples for the requested estimator")
            }
            Self::ChainLengthMismatch => f.write_str("chains must have equal non-empty lengths"),
        }
    }
}
impl std::error::Error for StatisticsError {}

/// Stable weighted online population moments.
#[derive(Clone, Debug, PartialEq)]
pub struct WeightedMoments {
    weight: f64,
    mean: f64,
    m2: f64,
}
impl WeightedMoments {
    /// Construct an empty accumulator.
    pub fn new() -> Result<Self, StatisticsError> {
        Ok(Self {
            weight: 0.0,
            mean: 0.0,
            m2: 0.0,
        })
    }
    /// Add one finite value with a non-negative weight.
    pub fn update(&mut self, value: f64, weight: f64) -> Result<(), StatisticsError> {
        if !value.is_finite() {
            return Err(StatisticsError::NonFinite("sample"));
        }
        validate_weight(weight)?;
        if weight == 0.0 {
            return Ok(());
        }
        if self.weight == 0.0 {
            self.weight = weight;
            self.mean = value;
            return Ok(());
        }
        let total = self.weight + weight;
        if !total.is_finite() {
            return Err(StatisticsError::NonFinite("accumulated weight"));
        }
        let delta = value - self.mean;
        let candidate_mean = self.mean + delta * weight / total;
        let candidate_m2 = self.m2 + delta * delta * self.weight * weight / total;
        if !candidate_m2.is_finite() || !candidate_mean.is_finite() {
            return Err(StatisticsError::NonFinite("accumulated moments"));
        }
        self.mean = candidate_mean;
        self.m2 = candidate_m2;
        self.weight = total;
        Ok(())
    }
    /// Merge another accumulator using the parallel-moments formula.
    pub fn merge(&mut self, other: &Self) -> Result<(), StatisticsError> {
        if other.weight == 0.0 {
            return Ok(());
        }
        if self.weight == 0.0 {
            *self = other.clone();
            return Ok(());
        }
        let total = self.weight + other.weight;
        if !total.is_finite() {
            return Err(StatisticsError::NonFinite("merged weight"));
        }
        let delta = other.mean - self.mean;
        let candidate_m2 = self.m2 + other.m2 + delta * delta * self.weight * other.weight / total;
        let candidate_mean = self.mean + delta * other.weight / total;
        if !candidate_m2.is_finite() || !candidate_mean.is_finite() {
            return Err(StatisticsError::NonFinite("merged moments"));
        }
        self.m2 = candidate_m2;
        self.mean = candidate_mean;
        self.weight = total;
        Ok(())
    }
    /// Return the weighted mean.
    pub fn mean(&self) -> Result<f64, StatisticsError> {
        if self.weight == 0.0 {
            Err(StatisticsError::InsufficientSamples)
        } else {
            Ok(self.mean)
        }
    }
    /// Return the weighted population variance.
    pub fn variance(&self) -> Result<f64, StatisticsError> {
        if self.weight == 0.0 {
            Err(StatisticsError::InsufficientSamples)
        } else {
            Ok(self.m2 / self.weight)
        }
    }
    /// Return the accumulated weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
}

/// Separate real and imaginary online moments for complex estimators.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexWeightedMoments {
    real: WeightedMoments,
    imaginary: WeightedMoments,
}
impl ComplexWeightedMoments {
    /// Construct empty complex moments.
    pub fn new() -> Result<Self, StatisticsError> {
        Ok(Self {
            real: WeightedMoments::new()?,
            imaginary: WeightedMoments::new()?,
        })
    }
    /// Add one complex value with a non-negative weight.
    pub fn update(&mut self, value: Complex64, weight: f64) -> Result<(), StatisticsError> {
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(StatisticsError::NonFinite("sample"));
        }
        validate_weight(weight)?;
        let mut real = self.real.clone();
        let mut imaginary = self.imaginary.clone();
        real.update(value.re, weight)?;
        imaginary.update(value.im, weight)?;
        self.real = real;
        self.imaginary = imaginary;
        Ok(())
    }
    /// Return the complex weighted mean.
    pub fn mean(&self) -> Result<Complex64, StatisticsError> {
        Ok(Complex64::new(self.real.mean()?, self.imaginary.mean()?))
    }
    /// Return the imaginary-component variance.
    pub fn imaginary_variance(&self) -> Result<f64, StatisticsError> {
        self.imaginary.variance()
    }
    /// Return the real-component variance.
    pub fn real_variance(&self) -> Result<f64, StatisticsError> {
        self.real.variance()
    }
}

fn validate_weight(weight: f64) -> Result<(), StatisticsError> {
    if !weight.is_finite() {
        Err(StatisticsError::NonFinite("weight"))
    } else if weight < 0.0 {
        Err(StatisticsError::InvalidWeight)
    } else {
        Ok(())
    }
}

/// Integrated autocorrelation estimate using the initial-positive-sequence rule.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AutocorrelationEstimate {
    integrated_time: f64,
    effective_sample_size: f64,
    algorithm: &'static str,
}
impl AutocorrelationEstimate {
    /// Return the integrated autocorrelation time.
    pub fn integrated_time(&self) -> f64 {
        self.integrated_time
    }
    /// Return `N / tau_int`.
    pub fn effective_sample_size(&self) -> f64 {
        self.effective_sample_size
    }
    /// Return the named autocorrelation algorithm.
    pub fn algorithm(&self) -> &'static str {
        self.algorithm
    }
}

/// Estimate autocorrelation time and effective sample size for one chain.
pub fn autocorrelation(
    samples: &[f64],
    max_lag: usize,
) -> Result<AutocorrelationEstimate, StatisticsError> {
    if samples.len() < 2 || max_lag == 0 {
        return Err(StatisticsError::InsufficientSamples);
    }
    if samples.iter().any(|sample| !sample.is_finite()) {
        return Err(StatisticsError::NonFinite("sample"));
    }
    let sum = samples.iter().sum::<f64>();
    if !sum.is_finite() {
        return Err(StatisticsError::NonFinite("autocorrelation mean"));
    }
    let mean = sum / samples.len() as f64;
    let gamma0 = samples
        .iter()
        .map(|sample| (sample - mean).powi(2))
        .sum::<f64>()
        / samples.len() as f64;
    if !gamma0.is_finite() {
        return Err(StatisticsError::NonFinite("autocorrelation variance"));
    }
    if gamma0 == 0.0 {
        return Ok(AutocorrelationEstimate {
            integrated_time: 1.0,
            effective_sample_size: samples.len() as f64,
            algorithm: "geyer_initial_positive_sequence_common_n_tau_floor_1",
        });
    }
    let mut tau = -1.0;
    let mut lag = 0;
    while lag < max_lag.min(samples.len() - 1) {
        let covariance = samples[..samples.len() - lag]
            .iter()
            .zip(&samples[lag..])
            .map(|(left, right)| (left - mean) * (right - mean))
            .sum::<f64>()
            / samples.len() as f64;
        let next_lag = lag + 1;
        let next_covariance = samples[..samples.len() - next_lag]
            .iter()
            .zip(&samples[next_lag..])
            .map(|(left, right)| (left - mean) * (right - mean))
            .sum::<f64>()
            / samples.len() as f64;
        let pair = covariance + next_covariance;
        if pair <= 0.0 {
            break;
        }
        tau += 2.0 * pair / gamma0;
        if !tau.is_finite() {
            return Err(StatisticsError::NonFinite("autocorrelation time"));
        }
        lag += 2;
    }
    Ok(AutocorrelationEstimate {
        integrated_time: tau.max(1.0),
        effective_sample_size: samples.len() as f64 / tau.max(1.0),
        algorithm: "geyer_initial_positive_sequence_common_n_tau_floor_1",
    })
}

/// Named split-free classic Gelman-Rubin diagnostic for equal-length chains.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RHatDiagnostic {
    value: f64,
    within_chain: f64,
    between_chain: f64,
    algorithm: &'static str,
}
impl RHatDiagnostic {
    /// Return R-hat.
    pub fn value(&self) -> f64 {
        self.value
    }
    /// Return within-chain variance `W`.
    pub fn within_chain_variance(&self) -> f64 {
        self.within_chain
    }
    /// Return between-chain variance `B`.
    pub fn between_chain_variance(&self) -> f64 {
        self.between_chain
    }
    /// Return the named algorithm.
    pub fn algorithm(&self) -> &'static str {
        self.algorithm
    }
}

/// Estimate split-free classic Gelman-Rubin R-hat for equal-length chains.
pub fn r_hat(chains: &[Vec<f64>]) -> Result<RHatDiagnostic, StatisticsError> {
    if chains.len() < 2 || chains.iter().any(|chain| chain.len() < 2) {
        return Err(StatisticsError::InsufficientSamples);
    }
    let length = chains[0].len();
    if chains.iter().any(|chain| chain.len() != length) {
        return Err(StatisticsError::ChainLengthMismatch);
    }
    if chains.iter().flatten().any(|sample| !sample.is_finite()) {
        return Err(StatisticsError::NonFinite("sample"));
    }
    let means = chains
        .iter()
        .map(|chain| chain.iter().sum::<f64>() / length as f64)
        .collect::<Vec<_>>();
    if means.iter().any(|mean| !mean.is_finite()) {
        return Err(StatisticsError::NonFinite("chain mean"));
    }
    let grand = means.iter().sum::<f64>() / means.len() as f64;
    if !grand.is_finite() {
        return Err(StatisticsError::NonFinite("grand mean"));
    }
    let within = chains
        .iter()
        .map(|chain| {
            let mean = chain.iter().sum::<f64>() / length as f64;
            chain
                .iter()
                .map(|sample| (sample - mean).powi(2))
                .sum::<f64>()
                / (length - 1) as f64
        })
        .sum::<f64>()
        / chains.len() as f64;
    let between = length as f64 * means.iter().map(|mean| (mean - grand).powi(2)).sum::<f64>()
        / (means.len() - 1) as f64;
    if !within.is_finite() || !between.is_finite() {
        return Err(StatisticsError::NonFinite("R-hat variance"));
    }
    if within == 0.0 {
        return Ok(RHatDiagnostic {
            value: if between == 0.0 { 1.0 } else { f64::INFINITY },
            within_chain: within,
            between_chain: between,
            algorithm: "gelman_rubin_classic_equal_length",
        });
    }
    Ok(RHatDiagnostic {
        value: (((length - 1) as f64 / length as f64 * within + between / length as f64) / within)
            .sqrt()
            .max(1.0),
        within_chain: within,
        between_chain: between,
        algorithm: "gelman_rubin_classic_equal_length",
    })
}

/// Weighted disorder-realization summary retaining identifiers and ensemble spread.
#[derive(Clone, Debug, PartialEq)]
pub struct DisorderSummary {
    realization_ids: Vec<String>,
    records: Vec<RealizationEstimate>,
    mean: f64,
    between_realization_variance: Option<f64>,
    sampling_variance: Option<f64>,
}
impl DisorderSummary {
    /// Return the weighted realization mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Return weighted between-realization population variance.
    pub fn between_realization_variance(&self) -> Option<f64> {
        self.between_realization_variance
    }
    /// Return realization identifiers in input order.
    pub fn realizations(&self) -> &[String] {
        &self.realization_ids
    }
    /// Return the retained fixed-realization estimates.
    pub fn records(&self) -> &[RealizationEstimate] {
        &self.records
    }
    /// Return weighted within-realization sampling variance when supplied.
    pub fn sampling_variance(&self) -> Option<f64> {
        self.sampling_variance
    }
}

/// One fixed-disorder estimate with optional within-realization sampling variance.
#[derive(Clone, Debug, PartialEq)]
pub struct RealizationEstimate {
    id: String,
    mean: f64,
    weight: f64,
    sampling_variance: Option<f64>,
}
impl RealizationEstimate {
    /// Construct a fixed-realization estimate.
    pub fn new(
        id: impl Into<String>,
        mean: f64,
        weight: f64,
        sampling_variance: Option<f64>,
    ) -> Result<Self, StatisticsError> {
        let id = id.into();
        if id.is_empty() {
            return Err(StatisticsError::InsufficientSamples);
        }
        if !mean.is_finite()
            || sampling_variance.is_some_and(|value| !value.is_finite() || value < 0.0)
        {
            return Err(StatisticsError::NonFinite("realization estimate"));
        }
        validate_weight(weight)?;
        Ok(Self {
            id,
            mean,
            weight,
            sampling_variance,
        })
    }
    /// Return realization ID.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Return fixed-realization mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Return aggregation weight.
    pub fn weight(&self) -> f64 {
        self.weight
    }
    /// Return optional within-realization sampling variance.
    pub fn sampling_variance(&self) -> Option<f64> {
        self.sampling_variance
    }
}

/// Aggregate fixed-realization estimates while retaining ensemble variation.
pub fn disorder_average(
    realizations: &[(&str, f64, f64)],
) -> Result<DisorderSummary, StatisticsError> {
    if realizations.is_empty() {
        return Err(StatisticsError::InsufficientSamples);
    }
    let records = realizations
        .iter()
        .map(|(id, mean, weight)| RealizationEstimate::new(*id, *mean, *weight, None))
        .collect::<Result<Vec<_>, _>>()?;
    disorder_average_with_uncertainty(&records)
}

/// Aggregate retained realization records with separate within-realization uncertainty.
pub fn disorder_average_with_uncertainty(
    records: &[RealizationEstimate],
) -> Result<DisorderSummary, StatisticsError> {
    if records.is_empty() {
        return Err(StatisticsError::InsufficientSamples);
    }
    for (index, record) in records.iter().enumerate() {
        if records[..index]
            .iter()
            .any(|previous| previous.id == record.id)
        {
            return Err(StatisticsError::DuplicateRealization);
        }
    }
    let mut moments = WeightedMoments::new()?;
    for record in records {
        moments.update(record.mean, record.weight)?;
    }
    if moments.weight() == 0.0 {
        return Err(StatisticsError::InvalidWeight);
    }
    let mean = moments.mean()?;
    let positive_records = records.iter().filter(|record| record.weight > 0.0);
    let positive: Vec<_> = positive_records.collect();
    let variance = (positive.len() >= 2)
        .then(|| moments.variance())
        .transpose()?;
    let sampling_variance = (!positive.is_empty()
        && positive
            .iter()
            .all(|record| record.sampling_variance.is_some()))
    .then(|| {
        let total_weight = positive.iter().map(|record| record.weight).sum::<f64>();
        positive
            .iter()
            .map(|record| {
                let normalized = record.weight / total_weight;
                normalized * normalized * record.sampling_variance.unwrap_or(0.0)
            })
            .sum::<f64>()
    });
    Ok(DisorderSummary {
        realization_ids: records.iter().map(|record| record.id.clone()).collect(),
        records: records.to_vec(),
        mean,
        between_realization_variance: variance,
        sampling_variance,
    })
}
