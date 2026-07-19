//! Variational statistics, TDVP linear solves, and integration kernels.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

use qslib_core::Complex64;
use std::fmt::{self, Display, Formatter};

/// Errors from weighted statistics and chain diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub enum StatisticsError {
    /// A weight, sample, or diagnostic parameter was non-finite.
    NonFinite(&'static str),
    /// A weight was negative or the total weight was zero.
    InvalidWeight,
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
        let delta = value - self.mean;
        self.mean += delta * weight / total;
        self.m2 += delta * delta * self.weight * weight / total;
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
        let delta = other.mean - self.mean;
        self.m2 += other.m2 + delta * delta * self.weight * other.weight / total;
        self.mean += delta * other.weight / total;
        self.weight = total;
        Ok(())
    }
    /// Return the weighted mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Return the weighted population variance.
    pub fn variance(&self) -> f64 {
        if self.weight == 0.0 {
            f64::NAN
        } else {
            self.m2 / self.weight
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
        self.real.update(value.re, weight)?;
        self.imaginary.update(value.im, weight)
    }
    /// Return the complex weighted mean.
    pub fn mean(&self) -> Complex64 {
        Complex64::new(self.real.mean(), self.imaginary.mean())
    }
    /// Return the imaginary-component variance.
    pub fn imaginary_variance(&self) -> f64 {
        self.imaginary.variance()
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
    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    let gamma0 = samples
        .iter()
        .map(|sample| (sample - mean).powi(2))
        .sum::<f64>()
        / samples.len() as f64;
    if gamma0 == 0.0 {
        return Ok(AutocorrelationEstimate {
            integrated_time: 1.0,
            effective_sample_size: samples.len() as f64,
        });
    }
    let mut tau = 1.0;
    for lag in 1..=max_lag.min(samples.len() - 1) {
        let covariance = samples[..samples.len() - lag]
            .iter()
            .zip(&samples[lag..])
            .map(|(left, right)| (left - mean) * (right - mean))
            .sum::<f64>()
            / (samples.len() - lag) as f64;
        let rho = covariance / gamma0;
        if rho <= 0.0 {
            break;
        }
        tau += 2.0 * rho;
    }
    Ok(AutocorrelationEstimate {
        integrated_time: tau,
        effective_sample_size: samples.len() as f64 / tau,
    })
}

/// Estimate split-free classic Gelman-Rubin R-hat for equal-length chains.
pub fn r_hat(chains: &[Vec<f64>]) -> Result<f64, StatisticsError> {
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
    let grand = means.iter().sum::<f64>() / means.len() as f64;
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
    if within == 0.0 {
        return Ok(if between == 0.0 { 1.0 } else { f64::INFINITY });
    }
    Ok(
        (((length - 1) as f64 / length as f64 * within + between / length as f64) / within)
            .sqrt()
            .max(1.0),
    )
}

/// Weighted disorder-realization summary retaining identifiers and ensemble spread.
#[derive(Clone, Debug, PartialEq)]
pub struct DisorderSummary {
    realization_ids: Vec<String>,
    mean: f64,
    between_realization_variance: f64,
}
impl DisorderSummary {
    /// Return the weighted realization mean.
    pub fn mean(&self) -> f64 {
        self.mean
    }
    /// Return weighted between-realization population variance.
    pub fn between_realization_variance(&self) -> f64 {
        self.between_realization_variance
    }
    /// Return realization identifiers in input order.
    pub fn realizations(&self) -> &[String] {
        &self.realization_ids
    }
}

/// Aggregate fixed-realization estimates while retaining ensemble variation.
pub fn disorder_average(
    realizations: &[(&str, f64, f64)],
) -> Result<DisorderSummary, StatisticsError> {
    if realizations.is_empty() {
        return Err(StatisticsError::InsufficientSamples);
    }
    let mut moments = WeightedMoments::new()?;
    for (_, value, weight) in realizations {
        moments.update(*value, *weight)?;
    }
    if moments.weight() == 0.0 {
        return Err(StatisticsError::InvalidWeight);
    }
    let variance = realizations
        .iter()
        .map(|(_, value, weight)| weight * (*value - moments.mean()).powi(2))
        .sum::<f64>()
        / moments.weight();
    Ok(DisorderSummary {
        realization_ids: realizations
            .iter()
            .map(|(id, _, _)| (*id).to_string())
            .collect(),
        mean: moments.mean(),
        between_realization_variance: variance,
    })
}
