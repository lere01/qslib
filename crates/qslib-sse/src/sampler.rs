//! Sign-safe diagonal insertion/removal sampling.

#![allow(deprecated)]

use crate::{
    BasisSseState, Operator, SseModel, SseModelError, ThermodynamicAccumulator,
    ThermodynamicResults,
};
use rand_core::RngCore;

/// Sampler construction or update failure.
#[derive(Clone, Debug, PartialEq)]
pub enum SamplerError {
    /// Invalid inverse temperature.
    InvalidBeta(f64),
    /// Model has no diagonal terms.
    NoDiagonalTerms,
    /// Invalid simulation configuration.
    InvalidConfig(&'static str),
    /// Model/state failure.
    Model(SseModelError),
}
impl From<SseModelError> for SamplerError {
    fn from(error: SseModelError) -> Self {
        Self::Model(error)
    }
}
impl std::fmt::Display for SamplerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidBeta(value) => write!(f, "beta must be positive and finite, got {value}"),
            Self::NoDiagonalTerms => f.write_str("SSE model has no diagonal terms"),
            Self::InvalidConfig(name) => write!(f, "invalid sampler configuration: {name}"),
            Self::Model(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for SamplerError {}

/// Insertion/removal acceptance counts from one diagonal sweep.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct DiagonalSweepStats {
    /// Proposed insertions.
    pub insertions_proposed: usize,
    /// Accepted insertions.
    pub insertions_accepted: usize,
    /// Proposed removals.
    pub removals_proposed: usize,
    /// Accepted removals.
    pub removals_accepted: usize,
}

/// Pair insertion/removal counts for off-diagonal vertices.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct OffDiagonalSweepStats {
    /// Pair proposals considered.
    pub proposals: usize,
    /// Pair proposals accepted.
    pub accepted: usize,
}

/// Thermalization and measurement controls.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SimulationConfig {
    /// Sweeps discarded before measurement.
    pub thermalization_sweeps: usize,
    /// Number of measured sweeps.
    pub measurement_sweeps: usize,
    /// Sweeps between measurements.
    pub sweeps_per_measurement: usize,
}
impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            thermalization_sweeps: 1_000,
            measurement_sweeps: 10_000,
            sweeps_per_measurement: 1,
        }
    }
}

/// Aggregate result from one sampled chain.
#[derive(Clone, Debug, PartialEq)]
pub struct SimulationResults {
    /// Expansion-order thermodynamics.
    pub thermodynamics: ThermodynamicResults,
    /// Aggregate diagonal update statistics.
    pub diagonal: DiagonalSweepStats,
    /// Aggregate off-diagonal pair statistics.
    pub off_diagonal: OffDiagonalSweepStats,
}

/// Fixed-length SSE sampler with sign-safe diagonal updates.
pub struct SseSampler<M, R> {
    model: M,
    state: BasisSseState,
    beta: f64,
    rng: R,
}
impl<M: SseModel, R: RngCore> SseSampler<M, R> {
    /// Construct a sampler after validating beta, diagonal support, and trace closure.
    pub fn new(model: M, state: BasisSseState, beta: f64, rng: R) -> Result<Self, SamplerError> {
        if !beta.is_finite() || beta <= 0.0 {
            return Err(SamplerError::InvalidBeta(beta));
        }
        if model.diagonal_term_indices().is_empty() {
            return Err(SamplerError::NoDiagonalTerms);
        }
        state.validate_trace(&model)?;
        Ok(Self {
            model,
            state,
            beta,
            rng,
        })
    }
    /// Borrow the model.
    pub fn model(&self) -> &M {
        &self.model
    }
    /// Borrow the current state.
    pub fn state(&self) -> &BasisSseState {
        &self.state
    }
    /// Return inverse temperature.
    pub fn beta(&self) -> f64 {
        self.beta
    }
    /// Return instantaneous energy estimator.
    pub fn energy_estimator(&self) -> f64 {
        self.model.energy_shift() - self.state.expansion_order() as f64 / self.beta
    }
    /// Perform one diagonal insertion/removal sweep.
    pub fn diagonal_sweep(&mut self) -> Result<DiagonalSweepStats, SamplerError> {
        let diagonal = self.model.diagonal_term_indices().to_vec();
        let mut stats = DiagonalSweepStats::default();
        let cutoff = self.state.operator_string().len();
        let initial_bits = self.state.bits().to_vec();
        let mut working_bits = initial_bits.clone();
        for position in 0..cutoff {
            let current = self.state.operator_string()[position];
            match current {
                Operator::Identity => {
                    stats.insertions_proposed += 1;
                    let choice = (self.rng.next_u64() as usize) % diagonal.len();
                    let term_index = diagonal[choice] as usize;
                    let weight = self.model.matrix_element(term_index, &working_bits)?;
                    let order = self.state.expansion_order();
                    let denominator = cutoff.saturating_sub(order).max(1) as f64;
                    let probability =
                        (self.beta * diagonal.len() as f64 * weight / denominator).min(1.0);
                    if random_unit(&mut self.rng) < probability {
                        self.state.operator_string_mut()[position] =
                            Operator::diagonal(diagonal[choice]);
                        stats.insertions_accepted += 1;
                    }
                }
                Operator::Diagonal(term_index) => {
                    stats.removals_proposed += 1;
                    let weight = self
                        .model
                        .matrix_element(term_index as usize, &working_bits)?;
                    if weight <= 0.0 {
                        self.state.operator_string_mut()[position] = Operator::identity();
                        stats.removals_accepted += 1;
                        continue;
                    }
                    let order = self.state.expansion_order();
                    let probability = ((cutoff.saturating_sub(order) + 1) as f64
                        / (self.beta * diagonal.len() as f64 * weight))
                        .min(1.0);
                    if random_unit(&mut self.rng) < probability {
                        self.state.operator_string_mut()[position] = Operator::identity();
                        stats.removals_accepted += 1;
                    }
                }
                Operator::OffDiagonal(term_index) => {
                    let value = self
                        .model
                        .matrix_element(term_index as usize, &working_bits)?;
                    if value <= 0.0 {
                        return Err(SamplerError::Model(
                            SseModelError::NonPositiveMatrixElement {
                                term_index: term_index as usize,
                                value,
                            },
                        ));
                    }
                    self.model
                        .apply_off_diagonal(term_index as usize, &mut working_bits)?;
                }
            }
        }
        if working_bits != initial_bits {
            return Err(SamplerError::Model(SseModelError::TraceNotClosed));
        }
        Ok(stats)
    }
    /// Propose symmetric pair insertions/removals of identical off-diagonal vertices.
    pub fn off_diagonal_pair_sweep(&mut self) -> Result<OffDiagonalSweepStats, SamplerError> {
        let off_diagonal = (0..self.model.num_terms())
            .filter(|&index| {
                self.model.operator_kind(index).ok() == Some(crate::OperatorKind::OffDiagonal)
            })
            .map(|index| index as u32)
            .collect::<Vec<_>>();
        let cutoff = self.state.operator_string().len();
        let mut stats = OffDiagonalSweepStats::default();
        if cutoff < 2 || off_diagonal.is_empty() {
            return Ok(stats);
        }
        for _ in 0..cutoff {
            let mut left = (self.rng.next_u64() as usize) % cutoff;
            let mut right = (self.rng.next_u64() as usize) % cutoff;
            if left == right {
                right = (right + 1) % cutoff;
            }
            if left > right {
                std::mem::swap(&mut left, &mut right);
            }
            let term_index = off_diagonal[(self.rng.next_u64() as usize) % off_diagonal.len()];
            let mut candidate = self.state.clone();
            let pair = (
                candidate.operator_string()[left],
                candidate.operator_string()[right],
            );
            let toggled = match pair {
                (Operator::Identity, Operator::Identity) => {
                    candidate.operator_string_mut()[left] = Operator::off_diagonal(term_index);
                    candidate.operator_string_mut()[right] = Operator::off_diagonal(term_index);
                    true
                }
                (Operator::OffDiagonal(a), Operator::OffDiagonal(b))
                    if a == term_index && b == term_index =>
                {
                    candidate.operator_string_mut()[left] = Operator::identity();
                    candidate.operator_string_mut()[right] = Operator::identity();
                    true
                }
                _ => false,
            };
            if !toggled {
                continue;
            }
            stats.proposals += 1;
            let old_log = configuration_log_weight(&self.state, &self.model, self.beta)?;
            let new_log = match configuration_log_weight(&candidate, &self.model, self.beta) {
                Ok(value) => value,
                Err(SamplerError::Model(SseModelError::TraceNotClosed))
                | Err(SamplerError::Model(SseModelError::NonPositiveMatrixElement { .. })) => {
                    continue;
                }
                Err(error) => return Err(error),
            };
            if random_unit(&mut self.rng) < (new_log - old_log).exp().min(1.0) {
                self.state = candidate;
                stats.accepted += 1;
            }
        }
        Ok(stats)
    }
    /// Run thermalization and measurement sweeps.
    pub fn run(&mut self, config: SimulationConfig) -> Result<SimulationResults, SamplerError> {
        if config.sweeps_per_measurement == 0 {
            return Err(SamplerError::InvalidConfig("sweeps_per_measurement"));
        }
        let mut diagonal = DiagonalSweepStats::default();
        let mut off_diagonal = OffDiagonalSweepStats::default();
        for _ in 0..config.thermalization_sweeps {
            let stats = self.diagonal_sweep()?;
            diagonal = add_stats(diagonal, stats);
            off_diagonal = add_off_stats(off_diagonal, self.off_diagonal_pair_sweep()?);
        }
        let mut accumulator = ThermodynamicAccumulator::default();
        for _ in 0..config.measurement_sweeps {
            for _ in 0..config.sweeps_per_measurement {
                let stats = self.diagonal_sweep()?;
                diagonal = add_stats(diagonal, stats);
                off_diagonal = add_off_stats(off_diagonal, self.off_diagonal_pair_sweep()?);
            }
            accumulator.record(self.state.expansion_order());
        }
        let thermodynamics = accumulator
            .results(self.beta, self.model.energy_shift(), self.model.num_sites())
            .ok_or(SamplerError::InvalidConfig("no measurements"))?;
        Ok(SimulationResults {
            thermodynamics,
            diagonal,
            off_diagonal,
        })
    }
}
fn add_stats(left: DiagonalSweepStats, right: DiagonalSweepStats) -> DiagonalSweepStats {
    DiagonalSweepStats {
        insertions_proposed: left.insertions_proposed + right.insertions_proposed,
        insertions_accepted: left.insertions_accepted + right.insertions_accepted,
        removals_proposed: left.removals_proposed + right.removals_proposed,
        removals_accepted: left.removals_accepted + right.removals_accepted,
    }
}
fn random_unit<R: RngCore>(rng: &mut R) -> f64 {
    (rng.next_u64() as f64) / (u64::MAX as f64)
}
fn add_off_stats(
    left: OffDiagonalSweepStats,
    right: OffDiagonalSweepStats,
) -> OffDiagonalSweepStats {
    OffDiagonalSweepStats {
        proposals: left.proposals + right.proposals,
        accepted: left.accepted + right.accepted,
    }
}
fn configuration_log_weight<M: SseModel>(
    state: &BasisSseState,
    model: &M,
    beta: f64,
) -> Result<f64, SamplerError> {
    let propagation = state.propagate(model)?;
    let order = state.expansion_order();
    let cutoff = state.operator_string().len();
    Ok(
        order as f64 * beta.ln() + factorial_log(cutoff - order) - factorial_log(cutoff)
            + propagation.log_weight,
    )
}
fn factorial_log(value: usize) -> f64 {
    (1..=value).map(|item| (item as f64).ln()).sum()
}
