//! Transactional Euler and Heun integration for checked flat parameter state.

use crate::DenseQgt;
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display, Formatter};

const SEED_ALGORITHM_VERSION: u32 = 1;

/// Errors produced by the evolution driver.
#[derive(Clone, Debug, PartialEq)]
pub enum EvolutionError {
    /// A configuration or state value is invalid.
    InvalidParameter(&'static str),
    /// A callback returned a domain-specific failure.
    Callback(String),
    /// A callback returned a vector with the wrong length.
    Shape {
        /// Expected vector length.
        expected: usize,
        /// Received vector length.
        actual: usize,
    },
    /// A non-finite state, velocity, or diagnostic was encountered.
    NonFinite(&'static str),
    /// Metadata could not be encoded or decoded.
    Serialization(String),
    /// A checkpoint and flat state use different layout fingerprints.
    LayoutMismatch,
    /// A checkpoint uses trajectory-changing integration controls that differ.
    ConfigurationMismatch,
    /// The adaptive step reached its lower bound before meeting tolerance.
    StepSizeUnderflow {
        /// Lower-bounded attempted step size.
        dt: f64,
        /// Error norm at that step size.
        error: f64,
    },
}
impl Display for EvolutionError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter(context) => write!(f, "invalid evolution {context}"),
            Self::Callback(message) => write!(f, "evolution callback failed: {message}"),
            Self::Shape { expected, actual } => {
                write!(
                    f,
                    "evolution vector shape mismatch: expected {expected}, got {actual}"
                )
            }
            Self::NonFinite(context) => write!(f, "non-finite evolution {context}"),
            Self::Serialization(message) => {
                write!(f, "evolution metadata serialization failed: {message}")
            }
            Self::LayoutMismatch => {
                f.write_str("evolution layout fingerprint does not match flat state")
            }
            Self::ConfigurationMismatch => {
                f.write_str("evolution checkpoint controls do not match configuration")
            }
            Self::StepSizeUnderflow { dt, error } => {
                write!(f, "adaptive step reached dt={dt} with error={error}")
            }
        }
    }
}
impl std::error::Error for EvolutionError {}

/// Reference one-step integration methods.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum IntegrationMethod {
    /// First-order forward Euler.
    Euler,
    /// Second-order explicit trapezoidal (Heun) integration.
    Heun,
}

/// Metric used to normalize an adaptive step error.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum ErrorMetric {
    /// Euclidean norm of the predictor-corrector difference.
    Euclidean,
    /// QGT norm supplied by the first velocity evaluation.
    Qgt,
}

/// Configuration for transactional evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvolutionConfig {
    /// Integration method.
    pub method: IntegrationMethod,
    /// Whether Heun attempts are accepted only when their error meets tolerance.
    pub adaptive: bool,
    /// Initial or fixed step size.
    pub dt: f64,
    /// Adaptive error tolerance in the configured metric.
    pub step_tolerance: f64,
    /// Lower adaptive step bound.
    pub dt_min: f64,
    /// Upper adaptive step bound.
    pub dt_max: f64,
    /// Multiplicative safety factor for accepted and rejected adaptive steps.
    pub safety_factor: f64,
    /// Adaptive error metric.
    pub error_metric: ErrorMetric,
    /// Master seed used to derive deterministic stage seeds.
    pub seed: u64,
}
impl Default for EvolutionConfig {
    fn default() -> Self {
        Self {
            method: IntegrationMethod::Heun,
            adaptive: false,
            dt: 0.01,
            step_tolerance: 1.0e-6,
            dt_min: 1.0e-8,
            dt_max: 1.0,
            safety_factor: 0.9,
            error_metric: ErrorMetric::Euclidean,
            seed: 0,
        }
    }
}
impl EvolutionConfig {
    fn validate(&self) -> Result<(), EvolutionError> {
        if self.adaptive && self.method == IntegrationMethod::Euler {
            return Err(EvolutionError::InvalidParameter(
                "adaptive Euler requires an error estimator",
            ));
        }
        if !self.dt.is_finite() || self.dt <= 0.0 {
            return Err(EvolutionError::InvalidParameter("step size"));
        }
        if !self.step_tolerance.is_finite() || self.step_tolerance <= 0.0 {
            return Err(EvolutionError::InvalidParameter("step tolerance"));
        }
        if !self.dt_min.is_finite() || self.dt_min <= 0.0 || self.dt_min > self.dt_max {
            return Err(EvolutionError::InvalidParameter("step bounds"));
        }
        if !self.dt_max.is_finite()
            || self.dt_max <= 0.0
            || !self.safety_factor.is_finite()
            || self.safety_factor <= 0.0
            || self.safety_factor > 1.0
        {
            return Err(EvolutionError::InvalidParameter("adaptive controls"));
        }
        Ok(())
    }
}

/// Checked flat parameter storage owned by the generic integrator.
#[derive(Clone, Debug, PartialEq)]
pub struct FlatState {
    fingerprint: String,
    values: Vec<f64>,
}
impl FlatState {
    /// Construct a finite non-empty flat state with a layout fingerprint.
    pub fn new(fingerprint: impl Into<String>, values: Vec<f64>) -> Result<Self, EvolutionError> {
        let fingerprint = fingerprint.into();
        if fingerprint.is_empty() || values.is_empty() {
            return Err(EvolutionError::InvalidParameter("flat state"));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(EvolutionError::NonFinite("flat state"));
        }
        Ok(Self {
            fingerprint,
            values,
        })
    }
    /// Return the layout fingerprint.
    pub fn fingerprint(&self) -> &str {
        &self.fingerprint
    }
    /// Return flat parameter values.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

/// Deterministic seed assigned to one velocity stage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct StageSeed {
    master_seed: u64,
    accepted_step: u64,
    rejected_attempt: u64,
    stage: u32,
    derived: u64,
}
impl StageSeed {
    /// Return the configured master seed.
    pub fn master_seed(self) -> u64 {
        self.master_seed
    }
    /// Return the accepted-step index at which this stage was evaluated.
    pub fn accepted_step(self) -> u64 {
        self.accepted_step
    }
    /// Return the number of rejected attempts already recorded.
    pub fn rejected_attempt(self) -> u64 {
        self.rejected_attempt
    }
    /// Return the integration stage number.
    pub fn stage(self) -> u32 {
        self.stage
    }
    /// Return the derived deterministic seed value.
    pub fn value(self) -> u64 {
        self.derived
    }
}

/// Velocity and optional QGT metric returned by a callback.
#[derive(Clone, Debug, PartialEq)]
pub struct Velocity {
    direction: Vec<f64>,
    qgt: Option<DenseQgt>,
}
impl Velocity {
    /// Construct a Euclidean velocity.
    pub fn new(direction: Vec<f64>) -> Self {
        Self {
            direction,
            qgt: None,
        }
    }
    /// Construct a velocity with its real QGT metric.
    pub fn with_qgt(direction: Vec<f64>, qgt: DenseQgt) -> Self {
        Self {
            direction,
            qgt: Some(qgt),
        }
    }
    /// Return the direction vector.
    pub fn direction(&self) -> &[f64] {
        &self.direction
    }
    /// Return the optional QGT metric.
    pub fn qgt(&self) -> Option<&DenseQgt> {
        self.qgt.as_ref()
    }
}

/// Serializable accepted-boundary evolution metadata.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct EvolutionMetadata {
    schema_version: u32,
    time: f64,
    proposed_dt: f64,
    accepted_steps: u64,
    rejected_steps: u64,
    seed: u64,
    layout_fingerprint: String,
    method: IntegrationMethod,
    error_metric: ErrorMetric,
    adaptive: bool,
    step_tolerance: f64,
    dt_min: f64,
    dt_max: f64,
    safety_factor: f64,
    seed_algorithm_version: u32,
}
impl EvolutionMetadata {
    /// Reconstruct accepted-boundary metadata from a durable checkpoint.
    ///
    /// `EvolutionDriver::from_parts` performs the complete configuration and
    /// numerical validity check after this value is assembled.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        time: f64,
        proposed_dt: f64,
        accepted_steps: u64,
        rejected_steps: u64,
        seed: u64,
        layout_fingerprint: impl Into<String>,
        method: IntegrationMethod,
        error_metric: ErrorMetric,
        adaptive: bool,
        step_tolerance: f64,
        dt_min: f64,
        dt_max: f64,
        safety_factor: f64,
        seed_algorithm_version: u32,
    ) -> Self {
        Self {
            schema_version: 1,
            time,
            proposed_dt,
            accepted_steps,
            rejected_steps,
            seed,
            layout_fingerprint: layout_fingerprint.into(),
            method,
            error_metric,
            adaptive,
            step_tolerance,
            dt_min,
            dt_max,
            safety_factor,
            seed_algorithm_version,
        }
    }
    /// Encode metadata as stable JSON for a later schema envelope.
    pub fn to_json(&self) -> Result<String, EvolutionError> {
        serde_json::to_string(self)
            .map_err(|error| EvolutionError::Serialization(error.to_string()))
    }
    /// Decode metadata and reject unknown fields.
    pub fn from_json(json: &str) -> Result<Self, EvolutionError> {
        serde_json::from_str(json).map_err(|error| EvolutionError::Serialization(error.to_string()))
    }
    /// Return physical time at the accepted boundary.
    pub fn time(&self) -> f64 {
        self.time
    }
    /// Return the proposed next step.
    pub fn proposed_dt(&self) -> f64 {
        self.proposed_dt
    }
    /// Return accepted-step count.
    pub fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }
    /// Return rejected-attempt count.
    pub fn rejected_steps(&self) -> u64 {
        self.rejected_steps
    }
    /// Return layout fingerprint.
    pub fn layout_fingerprint(&self) -> &str {
        &self.layout_fingerprint
    }
    /// Return the integration method.
    pub fn method(&self) -> IntegrationMethod {
        self.method
    }
    /// Return whether adaptive acceptance is enabled.
    pub fn adaptive(&self) -> bool {
        self.adaptive
    }
    /// Return the adaptive error metric.
    pub fn error_metric(&self) -> ErrorMetric {
        self.error_metric
    }
    /// Return the evolution master seed.
    pub fn seed(&self) -> u64 {
        self.seed
    }
    /// Return the adaptive tolerance.
    pub fn step_tolerance(&self) -> f64 {
        self.step_tolerance
    }
    /// Return the lower step bound.
    pub fn dt_min(&self) -> f64 {
        self.dt_min
    }
    /// Return the upper step bound.
    pub fn dt_max(&self) -> f64 {
        self.dt_max
    }
    /// Return the adaptive safety factor.
    pub fn safety_factor(&self) -> f64 {
        self.safety_factor
    }
    /// Return the seed derivation algorithm version.
    pub fn seed_algorithm_version(&self) -> u32 {
        self.seed_algorithm_version
    }
}

/// Accepted-boundary state of an evolution trajectory.
#[derive(Clone, Debug, PartialEq)]
pub struct EvolutionState {
    flat_state: FlatState,
    metadata: EvolutionMetadata,
}
impl EvolutionState {
    /// Return the checked flat state.
    pub fn flat_state(&self) -> &FlatState {
        &self.flat_state
    }
    /// Return accepted-boundary metadata.
    pub fn metadata(&self) -> &EvolutionMetadata {
        &self.metadata
    }
}

/// Result of one accepted integration step.
#[derive(Clone, Debug, PartialEq)]
pub struct StepOutcome {
    dt: f64,
    error_norm: f64,
    rejected_attempts: u64,
}
impl StepOutcome {
    /// Return the accepted step size.
    pub fn dt(&self) -> f64 {
        self.dt
    }
    /// Return the final attempted error norm.
    pub fn error_norm(&self) -> f64 {
        self.error_norm
    }
    /// Return rejected attempts incurred before acceptance.
    pub fn rejected_attempts(&self) -> u64 {
        self.rejected_attempts
    }
}

/// Generic transactional integrator over a checked flat state.
#[derive(Clone, Debug, PartialEq)]
pub struct EvolutionDriver {
    state: EvolutionState,
    config: EvolutionConfig,
}
impl EvolutionDriver {
    /// Construct a trajectory at time zero.
    pub fn new(flat_state: FlatState, config: EvolutionConfig) -> Result<Self, EvolutionError> {
        config.validate()?;
        let metadata = EvolutionMetadata {
            schema_version: 1,
            time: 0.0,
            proposed_dt: config.dt.clamp(config.dt_min, config.dt_max),
            accepted_steps: 0,
            rejected_steps: 0,
            seed: config.seed,
            layout_fingerprint: flat_state.fingerprint.clone(),
            method: config.method,
            error_metric: config.error_metric,
            adaptive: config.adaptive,
            step_tolerance: config.step_tolerance,
            dt_min: config.dt_min,
            dt_max: config.dt_max,
            safety_factor: config.safety_factor,
            seed_algorithm_version: SEED_ALGORITHM_VERSION,
        };
        Ok(Self {
            state: EvolutionState {
                flat_state,
                metadata,
            },
            config,
        })
    }
    /// Restore an accepted-boundary state after checking layout compatibility.
    pub fn from_parts(
        flat_state: FlatState,
        metadata: EvolutionMetadata,
        config: EvolutionConfig,
    ) -> Result<Self, EvolutionError> {
        config.validate()?;
        if flat_state.fingerprint != metadata.layout_fingerprint {
            return Err(EvolutionError::LayoutMismatch);
        }
        if config.seed != metadata.seed
            || config.method != metadata.method
            || config.error_metric != metadata.error_metric
            || config.adaptive != metadata.adaptive
            || config.step_tolerance != metadata.step_tolerance
            || config.dt_min != metadata.dt_min
            || config.dt_max != metadata.dt_max
            || config.safety_factor != metadata.safety_factor
            || metadata.seed_algorithm_version != SEED_ALGORITHM_VERSION
        {
            return Err(EvolutionError::ConfigurationMismatch);
        }
        if metadata.schema_version != 1
            || !metadata.time.is_finite()
            || !metadata.proposed_dt.is_finite()
            || metadata.proposed_dt <= 0.0
            || metadata.proposed_dt < config.dt_min
            || metadata.proposed_dt > config.dt_max
        {
            return Err(EvolutionError::NonFinite("evolution metadata"));
        }
        Ok(Self {
            state: EvolutionState {
                flat_state,
                metadata,
            },
            config,
        })
    }
    /// Return accepted-boundary state.
    pub fn state(&self) -> &EvolutionState {
        &self.state
    }

    /// Advance through one accepted step, retrying rejected adaptive attempts transactionally.
    pub fn advance<F>(&mut self, mut callback: F) -> Result<StepOutcome, EvolutionError>
    where
        F: FnMut(&[f64], f64, StageSeed) -> Result<Velocity, EvolutionError>,
    {
        let mut dt = self.state.metadata.proposed_dt;
        let rejected_before = self.state.metadata.rejected_steps;
        loop {
            let attempt = self.attempt(dt, self.state.metadata.rejected_steps, &mut callback)?;
            if !self.config.adaptive || attempt.error_norm <= self.config.step_tolerance {
                let next_time = self.state.metadata.time + dt;
                if !next_time.is_finite() {
                    return Err(EvolutionError::NonFinite("accepted time"));
                }
                let next_accepted = self
                    .state
                    .metadata
                    .accepted_steps
                    .checked_add(1)
                    .ok_or(EvolutionError::InvalidParameter("accepted-step count"))?;
                let next_dt = if self.config.adaptive {
                    adaptive_next_dt(
                        dt,
                        attempt.error_norm,
                        self.config.step_tolerance,
                        self.config.dt_min,
                        self.config.dt_max,
                        self.config.safety_factor,
                    )
                } else {
                    dt
                };
                let rejected_attempts = self
                    .state
                    .metadata
                    .rejected_steps
                    .checked_sub(rejected_before)
                    .ok_or(EvolutionError::InvalidParameter("rejected-step count"))?;
                self.state.flat_state.values = attempt.candidate;
                self.state.metadata.time = next_time;
                self.state.metadata.accepted_steps = next_accepted;
                self.state.metadata.proposed_dt = next_dt;
                return Ok(StepOutcome {
                    dt,
                    error_norm: attempt.error_norm,
                    rejected_attempts,
                });
            }
            let next_rejected = self
                .state
                .metadata
                .rejected_steps
                .checked_add(1)
                .ok_or(EvolutionError::InvalidParameter("rejected-step count"))?;
            let next = adaptive_next_dt(
                dt,
                attempt.error_norm,
                self.config.step_tolerance,
                self.config.dt_min,
                self.config.dt_max,
                self.config.safety_factor,
            );
            if next >= dt && dt <= self.config.dt_min {
                return Err(EvolutionError::StepSizeUnderflow {
                    dt,
                    error: attempt.error_norm,
                });
            }
            self.state.metadata.rejected_steps = next_rejected;
            dt = next;
            self.state.metadata.proposed_dt = dt;
        }
    }

    /// Advance a fixed number of accepted steps and invoke `observer` only after each acceptance.
    pub fn run<F, O>(
        &mut self,
        accepted_steps: u64,
        mut observer: O,
        mut callback: F,
    ) -> Result<Vec<StepOutcome>, EvolutionError>
    where
        F: FnMut(&[f64], f64, StageSeed) -> Result<Velocity, EvolutionError>,
        O: FnMut(&EvolutionState),
    {
        let mut outcomes = Vec::with_capacity(accepted_steps as usize);
        for _ in 0..accepted_steps {
            let outcome = self.advance(&mut callback)?;
            observer(&self.state);
            outcomes.push(outcome);
        }
        Ok(outcomes)
    }

    fn attempt<F>(
        &self,
        dt: f64,
        rejected_attempt: u64,
        callback: &mut F,
    ) -> Result<Attempt, EvolutionError>
    where
        F: FnMut(&[f64], f64, StageSeed) -> Result<Velocity, EvolutionError>,
    {
        let base = &self.state.flat_state.values;
        let time = self.state.metadata.time;
        let first = callback(
            base,
            time,
            stage_seed(
                self.config.seed,
                self.state.metadata.accepted_steps,
                rejected_attempt,
                0,
            ),
        )?;
        validate_velocity(&first, base.len())?;
        let mut euler = base.clone();
        add_scaled(&mut euler, first.direction(), dt)?;
        if self.config.method == IntegrationMethod::Euler {
            return Ok(Attempt {
                candidate: euler,
                error_norm: 0.0,
            });
        }
        let second = callback(
            &euler,
            time + dt,
            stage_seed(
                self.config.seed,
                self.state.metadata.accepted_steps,
                rejected_attempt,
                1,
            ),
        )?;
        validate_velocity(&second, base.len())?;
        let mut corrected = base.clone();
        for (index, value) in corrected.iter_mut().enumerate() {
            *value += 0.5 * dt * (first.direction()[index] + second.direction()[index]);
            if !value.is_finite() {
                return Err(EvolutionError::NonFinite("corrected state"));
            }
        }
        let difference = corrected
            .iter()
            .zip(&euler)
            .map(|(left, right)| left - right)
            .collect::<Vec<_>>();
        let error_norm = match self.config.error_metric {
            ErrorMetric::Euclidean => l2_norm(&difference),
            ErrorMetric::Qgt => {
                let qgt = first.qgt().ok_or(EvolutionError::InvalidParameter(
                    "QGT error metric without QGT",
                ))?;
                if !qgt_is_positive_semidefinite_scale_relative(qgt)? {
                    return Err(EvolutionError::InvalidParameter("indefinite QGT metric"));
                }
                let applied = qgt
                    .matvec(&difference)
                    .map_err(|_| EvolutionError::InvalidParameter("QGT error metric"))?;
                qgt_norm(&difference, &applied)?.sqrt()
            }
        };
        if !error_norm.is_finite() {
            return Err(EvolutionError::NonFinite("step error"));
        }
        Ok(Attempt {
            candidate: corrected,
            error_norm,
        })
    }
}

struct Attempt {
    candidate: Vec<f64>,
    error_norm: f64,
}

fn validate_velocity(velocity: &Velocity, expected: usize) -> Result<(), EvolutionError> {
    if velocity.direction.len() != expected {
        return Err(EvolutionError::Shape {
            expected,
            actual: velocity.direction.len(),
        });
    }
    if velocity.direction.iter().any(|value| !value.is_finite()) {
        return Err(EvolutionError::NonFinite("velocity"));
    }
    if let Some(qgt) = &velocity.qgt {
        if qgt.parameters() != expected {
            return Err(EvolutionError::Shape {
                expected,
                actual: qgt.parameters(),
            });
        }
    }
    Ok(())
}
fn add_scaled(values: &mut [f64], direction: &[f64], scale: f64) -> Result<(), EvolutionError> {
    for (value, increment) in values.iter_mut().zip(direction) {
        *value += scale * *increment;
        if !value.is_finite() {
            return Err(EvolutionError::NonFinite("updated state"));
        }
    }
    Ok(())
}
fn l2_norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |norm, value| norm.hypot(*value))
}
fn adaptive_next_dt(
    dt: f64,
    error: f64,
    tolerance: f64,
    dt_min: f64,
    dt_max: f64,
    safety_factor: f64,
) -> f64 {
    let growth = if error == 0.0 {
        2.0
    } else {
        (tolerance / error).sqrt()
    };
    (dt * safety_factor * growth).clamp(dt_min, dt_max)
}
fn qgt_norm(left: &[f64], right: &[f64]) -> Result<f64, EvolutionError> {
    let value = left.iter().zip(right).map(|(a, b)| a * b).sum::<f64>();
    let scale = left
        .iter()
        .zip(right)
        .map(|(a, b)| (a * b).abs())
        .sum::<f64>()
        .max(f64::MIN_POSITIVE);
    if !value.is_finite() || value < -1.0e-12 * scale {
        return Err(EvolutionError::NonFinite("metric norm"));
    }
    Ok(value.max(0.0))
}
fn qgt_is_positive_semidefinite_scale_relative(qgt: &DenseQgt) -> Result<bool, EvolutionError> {
    let spectrum = qgt
        .spectrum()
        .map_err(|_| EvolutionError::InvalidParameter("QGT spectrum"))?;
    let scale = spectrum
        .eigenvalues
        .iter()
        .map(|value| value.abs())
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);
    Ok(spectrum
        .eigenvalues
        .iter()
        .all(|value| *value >= -1.0e-10 * scale))
}
fn stage_seed(
    master_seed: u64,
    accepted_step: u64,
    rejected_attempt: u64,
    stage: u32,
) -> StageSeed {
    let mut value = master_seed ^ 0x9e37_79b9_7f4a_7c15;
    for component in [accepted_step, u64::from(stage)] {
        value ^= component
            .wrapping_add(0x9e37_79b9_7f4a_7c15)
            .rotate_left(17);
        value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value ^= value >> 31;
    }
    StageSeed {
        master_seed,
        accepted_step,
        rejected_attempt,
        stage,
        derived: value,
    }
}
