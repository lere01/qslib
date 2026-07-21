//! Sign-safe diagonal insertion/removal sampling.

use crate::{
    BasisSseState, Operator, SseModel, SseModelError, ThermodynamicAccumulator,
    ThermodynamicResults,
};
use rand_chacha::ChaCha20Rng;
use rand_core::{Rng, SeedableRng};
use std::sync::{Arc, mpsc};

/// Sampler construction or update failure.
#[derive(Clone, Debug, PartialEq)]
pub enum SamplerError {
    /// Invalid inverse temperature.
    InvalidBeta(f64),
    /// Model has no diagonal terms.
    NoDiagonalTerms,
    /// A TFIM cluster update was requested for an incompatible model.
    UnsupportedClusterUpdate,
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
            Self::UnsupportedClusterUpdate => {
                f.write_str("the model does not support the TFIM cluster update")
            }
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

/// Structural and acceptance counts from cluster updates.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ClusterSweepStats {
    /// Connected components identified by the linked-cluster breakup.
    pub clusters: usize,
    /// Components selected for flipping.
    pub flipped_clusters: usize,
    /// Diagonal/off-diagonal partner vertices toggled.
    pub vertices_toggled: usize,
    /// Cluster proposals attempted.
    pub proposals: usize,
    /// Cluster proposals accepted.
    pub proposals_accepted: usize,
}

/// Boundary basis-state flip counts from one basis sweep.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BasisSweepStats {
    /// Site-flip proposals considered.
    pub proposals: usize,
    /// Accepted site flips.
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
    /// Aggregate boundary basis-state update statistics.
    pub basis: BasisSweepStats,
}

/// Adaptive-cutoff SSE sampler with sign-safe local updates.
pub struct SseSampler<M, R> {
    model: M,
    state: BasisSseState,
    beta: f64,
    rng: R,
}
impl<M: SseModel, R: Rng> SseSampler<M, R> {
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
    /// Grow the cutoff when the identity headroom is below 25 percent.
    pub fn ensure_operator_headroom(&mut self) -> bool {
        let cutoff = self.state.operator_string().len();
        if !self.identity_headroom_low() {
            return false;
        }
        self.state
            .grow_operator_string(cutoff.saturating_mul(2).max(cutoff.saturating_add(1)));
        true
    }
    fn identity_headroom_low(&self) -> bool {
        let cutoff = self.state.operator_string().len();
        let order = self.state.expansion_order();
        let minimum_empty = (cutoff / 4).max(16);
        cutoff.saturating_sub(order) < minimum_empty
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
    /// Propose independent symmetric flips of the imaginary-time boundary state.
    ///
    /// This update is necessary for trace sampling: holding the boundary bits
    /// fixed would sample only one diagonal sector, which is exact only for
    /// special symmetric initial states.
    pub fn basis_state_sweep(&mut self) -> Result<BasisSweepStats, SamplerError> {
        let mut stats = BasisSweepStats::default();
        let old_log = configuration_log_weight(&self.state, &self.model, self.beta)?;
        let site = (self.rng.next_u64() as usize) % self.state.bits().len();
        stats.proposals = 1;
        let mut candidate = self.state.clone();
        let bit = candidate
            .bits_mut()
            .get_mut(site)
            .ok_or(SamplerError::Model(SseModelError::InvalidLength {
                expected: self.model.num_sites(),
                actual: self.state.bits().len(),
            }))?;
        *bit = match *bit {
            qslib_core::BasisBit::Zero => qslib_core::BasisBit::One,
            qslib_core::BasisBit::One => qslib_core::BasisBit::Zero,
        };
        let new_log = match configuration_log_weight(&candidate, &self.model, self.beta) {
            Ok(value) => value,
            Err(SamplerError::Model(SseModelError::TraceNotClosed))
            | Err(SamplerError::Model(SseModelError::NonPositiveMatrixElement { .. })) => {
                return Ok(stats);
            }
            Err(error) => return Err(error),
        };
        if random_unit(&mut self.rng) < (new_log - old_log).exp().min(1.0) {
            self.state = candidate;
            stats.accepted += 1;
        }
        Ok(stats)
    }
    /// Perform the deterministic transverse-field Ising linked-cluster update.
    ///
    /// Every vertex is broken up into world-line legs, legs are joined through
    /// bond vertices and around the imaginary-time boundary, and each connected
    /// component is flipped with probability one half. Flipping a component
    /// toggles the transverse `SiteConstant`/`SpinFlip` partner vertices it
    /// crosses, which is weight-neutral because partners share one matrix
    /// element. Models must advertise
    /// [`SseModel::supports_tfim_cluster_update`].
    pub fn tfim_cluster_sweep(&mut self) -> Result<ClusterSweepStats, SamplerError> {
        if !self.model.supports_tfim_cluster_update() {
            return Err(SamplerError::UnsupportedClusterUpdate);
        }
        self.cluster_sweep(false)
    }
    /// Perform a trace-preserving cluster proposal followed by a Metropolis
    /// correction for occupation-dependent Rydberg diagonal weights.
    ///
    /// This global proposal is retained primarily as a correctness reference;
    /// its acceptance may be poor for large or strongly interacting systems.
    pub fn rydberg_global_cluster_sweep(&mut self) -> Result<ClusterSweepStats, SamplerError> {
        self.cluster_sweep(true)
    }
    fn cluster_sweep(
        &mut self,
        metropolis_correction: bool,
    ) -> Result<ClusterSweepStats, SamplerError> {
        let original = self.state.clone();
        let original_log_weight = if metropolis_correction {
            self.state.propagate(&self.model)?.log_weight
        } else {
            0.0
        };

        let num_sites = self.model.num_sites();
        let mut first_leg = vec![None; num_sites];
        let mut last_leg = vec![None; num_sites];
        let mut union_find = UnionFind::default();
        let mut site_vertices = Vec::new();
        for (position, operator) in self.state.operator_string().iter().enumerate() {
            let Some(term_index) = operator.term_index() else {
                continue;
            };
            let term = *self
                .model
                .term(term_index)
                .ok_or(SseModelError::InvalidTermIndex {
                    term_index,
                    num_terms: self.model.num_terms(),
                })?;
            match term.sites() {
                (site, None) => {
                    let incoming = union_find.add();
                    let outgoing = union_find.add();
                    link_world_line(
                        site as usize,
                        incoming,
                        outgoing,
                        &mut first_leg,
                        &mut last_leg,
                        &mut union_find,
                    );
                    if self.model.transverse_partner(term_index).is_some() {
                        site_vertices.push(SiteVertex {
                            position,
                            term_index,
                            incoming,
                            outgoing,
                        });
                    } else {
                        union_find.union(incoming, outgoing);
                    }
                }
                (site_i, Some(site_j)) => {
                    let i_in = union_find.add();
                    let i_out = union_find.add();
                    let j_in = union_find.add();
                    let j_out = union_find.add();
                    link_world_line(
                        site_i as usize,
                        i_in,
                        i_out,
                        &mut first_leg,
                        &mut last_leg,
                        &mut union_find,
                    );
                    link_world_line(
                        site_j as usize,
                        j_in,
                        j_out,
                        &mut first_leg,
                        &mut last_leg,
                        &mut union_find,
                    );
                    union_find.union(i_in, i_out);
                    union_find.union(i_in, j_in);
                    union_find.union(i_in, j_out);
                }
            }
        }
        for site in 0..num_sites {
            if let (Some(first), Some(last)) = (first_leg[site], last_leg[site]) {
                union_find.union(first, last);
            }
        }

        let mut flip = vec![false; union_find.len()];
        let mut assigned = vec![false; union_find.len()];
        let mut stats = ClusterSweepStats::default();
        for leg in 0..union_find.len() {
            let root = union_find.find(leg);
            if !assigned[root] {
                assigned[root] = true;
                flip[root] = random_bool(&mut self.rng);
                stats.clusters += 1;
                stats.flipped_clusters += usize::from(flip[root]);
            }
        }
        for (site, &first) in first_leg.iter().enumerate() {
            match first {
                Some(leg) if flip[union_find.find(leg)] => {
                    flip_bit(&mut self.state.bits_mut()[site]);
                }
                None => {
                    stats.clusters += 1;
                    if random_bool(&mut self.rng) {
                        flip_bit(&mut self.state.bits_mut()[site]);
                        stats.flipped_clusters += 1;
                    }
                }
                _ => {}
            }
        }
        for vertex in site_vertices {
            if flip[union_find.find(vertex.incoming)] == flip[union_find.find(vertex.outgoing)] {
                continue;
            }
            let partner = self
                .model
                .transverse_partner(vertex.term_index)
                .ok_or(SamplerError::UnsupportedClusterUpdate)?;
            self.state.operator_string_mut()[vertex.position] =
                match self.model.operator_kind(partner as usize)? {
                    crate::OperatorKind::Diagonal => Operator::diagonal(partner),
                    crate::OperatorKind::OffDiagonal => Operator::off_diagonal(partner),
                    crate::OperatorKind::Identity => {
                        return Err(SamplerError::Model(SseModelError::InvalidOperatorKind));
                    }
                };
            stats.vertices_toggled += 1;
        }

        stats.proposals = 1;
        if metropolis_correction {
            let proposed = match self.state.propagate(&self.model) {
                Ok(value) => value,
                Err(SseModelError::NonPositiveMatrixElement { .. }) => {
                    self.state = original;
                    stats.flipped_clusters = 0;
                    stats.vertices_toggled = 0;
                    return Ok(stats);
                }
                Err(error) => return Err(SamplerError::Model(error)),
            };
            if !proposed.trace_closed {
                return Err(SamplerError::Model(SseModelError::TraceNotClosed));
            }
            let log_acceptance = proposed.log_weight - original_log_weight;
            if random_unit(&mut self.rng) < log_acceptance.min(0.0).exp() {
                stats.proposals_accepted = 1;
            } else {
                self.state = original;
                stats.flipped_clusters = 0;
                stats.vertices_toggled = 0;
            }
        } else {
            self.state.validate_trace(&self.model)?;
            stats.proposals_accepted = 1;
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
        let mut basis = BasisSweepStats::default();
        for _ in 0..config.thermalization_sweeps {
            self.ensure_operator_headroom();
            basis = add_basis_stats(basis, self.basis_state_sweep()?);
            let stats = self.diagonal_sweep()?;
            diagonal = add_stats(diagonal, stats);
            off_diagonal = add_off_stats(off_diagonal, self.off_diagonal_pair_sweep()?);
        }
        let mut accumulator = ThermodynamicAccumulator::default();
        let mut measured = 0;
        while measured < config.measurement_sweeps {
            if self.identity_headroom_low() {
                return Err(SamplerError::InvalidConfig("operator_string_cutoff"));
            }
            for _ in 0..config.sweeps_per_measurement {
                basis = add_basis_stats(basis, self.basis_state_sweep()?);
                let stats = self.diagonal_sweep()?;
                diagonal = add_stats(diagonal, stats);
                off_diagonal = add_off_stats(off_diagonal, self.off_diagonal_pair_sweep()?);
            }
            accumulator.record(self.state.expansion_order());
            measured += 1;
        }
        let thermodynamics = accumulator
            .results(self.beta, self.model.energy_shift(), self.model.num_sites())
            .ok_or(SamplerError::InvalidConfig("no measurements"))?;
        Ok(SimulationResults {
            thermodynamics,
            diagonal,
            off_diagonal,
            basis,
        })
    }
}

/// Run independent logical chains with deterministic scheduling-independent
/// streams. Returned results are sorted by logical chain index.
pub fn run_parallel_chains<M>(
    model: M,
    state: BasisSseState,
    beta: f64,
    config: SimulationConfig,
    master_seed: u64,
    chain_count: usize,
    worker_count: usize,
) -> Result<Vec<SimulationResults>, SamplerError>
where
    M: SseModel + Clone + Send + Sync + 'static,
{
    if chain_count == 0 {
        return Err(SamplerError::InvalidConfig("chain_count"));
    }
    if worker_count == 0 {
        return Err(SamplerError::InvalidConfig("worker_count"));
    }
    let workers = worker_count.min(chain_count);
    let model = Arc::new(model);
    let state = Arc::new(state);
    let (sender, receiver) = mpsc::channel();
    std::thread::scope(|scope| {
        for worker in 0..workers {
            let model = Arc::clone(&model);
            let state = Arc::clone(&state);
            let sender = sender.clone();
            scope.spawn(move || {
                for chain_index in (worker..chain_count).step_by(workers) {
                    let result = SseSampler::new(
                        (*model).clone(),
                        (*state).clone(),
                        beta,
                        ChaCha20Rng::from_seed(crate::derive_chain_seed(
                            master_seed,
                            chain_index as u64,
                        )),
                    )
                    .and_then(|mut sampler| sampler.run(config));
                    if sender.send((chain_index, result)).is_err() {
                        return;
                    }
                }
            });
        }
        drop(sender);
        let mut results: Vec<Option<SimulationResults>> = (0..chain_count).map(|_| None).collect();
        for _ in 0..chain_count {
            let (chain_index, result) = receiver
                .recv()
                .map_err(|_| SamplerError::InvalidConfig("chain worker"))?;
            results[chain_index] = Some(result?);
        }
        Ok(results.into_iter().flatten().collect())
    })
}
fn add_stats(left: DiagonalSweepStats, right: DiagonalSweepStats) -> DiagonalSweepStats {
    DiagonalSweepStats {
        insertions_proposed: left.insertions_proposed + right.insertions_proposed,
        insertions_accepted: left.insertions_accepted + right.insertions_accepted,
        removals_proposed: left.removals_proposed + right.removals_proposed,
        removals_accepted: left.removals_accepted + right.removals_accepted,
    }
}
fn random_unit<R: Rng>(rng: &mut R) -> f64 {
    (rng.next_u64() as f64) / (u64::MAX as f64)
}
fn random_bool<R: Rng>(rng: &mut R) -> bool {
    rng.next_u64() & 1 == 1
}
fn flip_bit(bit: &mut qslib_core::BasisBit) {
    *bit = match *bit {
        qslib_core::BasisBit::Zero => qslib_core::BasisBit::One,
        qslib_core::BasisBit::One => qslib_core::BasisBit::Zero,
    };
}

/// A toggleable transverse vertex recorded during the linked-cluster breakup.
#[derive(Clone, Copy, Debug)]
struct SiteVertex {
    position: usize,
    term_index: usize,
    incoming: usize,
    outgoing: usize,
}

#[derive(Default)]
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}
impl UnionFind {
    fn add(&mut self) -> usize {
        let index = self.parent.len();
        self.parent.push(index);
        self.rank.push(0);
        index
    }
    fn len(&self) -> usize {
        self.parent.len()
    }
    fn find(&mut self, index: usize) -> usize {
        if self.parent[index] != index {
            self.parent[index] = self.find(self.parent[index]);
        }
        self.parent[index]
    }
    fn union(&mut self, left: usize, right: usize) {
        let mut left_root = self.find(left);
        let mut right_root = self.find(right);
        if left_root == right_root {
            return;
        }
        if self.rank[left_root] < self.rank[right_root] {
            std::mem::swap(&mut left_root, &mut right_root);
        }
        self.parent[right_root] = left_root;
        if self.rank[left_root] == self.rank[right_root] {
            self.rank[left_root] += 1;
        }
    }
}

fn link_world_line(
    site: usize,
    incoming: usize,
    outgoing: usize,
    first_leg: &mut [Option<usize>],
    last_leg: &mut [Option<usize>],
    union_find: &mut UnionFind,
) {
    if let Some(previous) = last_leg[site] {
        union_find.union(previous, incoming);
    } else {
        first_leg[site] = Some(incoming);
    }
    last_leg[site] = Some(outgoing);
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
fn add_basis_stats(left: BasisSweepStats, right: BasisSweepStats) -> BasisSweepStats {
    BasisSweepStats {
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
