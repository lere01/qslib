//! Thermodynamic estimators from SSE expansion-order moments.

/// Constant-memory expansion-order accumulator.
#[derive(Clone, Copy, Debug, Default)]
pub struct ThermodynamicAccumulator {
    samples: u64,
    sum: f64,
    sum_squared: f64,
}
impl ThermodynamicAccumulator {
    /// Record one expansion order.
    pub fn record(&mut self, expansion_order: usize) {
        let value = expansion_order as f64;
        self.samples += 1;
        self.sum += value;
        self.sum_squared += value * value;
    }
    /// Return sample count.
    pub fn samples(&self) -> u64 {
        self.samples
    }
    /// Convert moments to physical thermodynamic estimates.
    pub fn results(
        &self,
        beta: f64,
        energy_shift: f64,
        num_sites: usize,
    ) -> Option<ThermodynamicResults> {
        if self.samples == 0 || !beta.is_finite() || beta <= 0.0 || num_sites == 0 {
            return None;
        }
        let count = self.samples as f64;
        let mean = self.sum / count;
        let second = self.sum_squared / count;
        let energy = energy_shift - mean / beta;
        let variance = if self.samples > 1 {
            ((second - mean * mean) * count / (count - 1.0)).max(0.0)
        } else {
            0.0
        };
        let heat_capacity = second - mean * mean - mean;
        Some(ThermodynamicResults {
            samples: self.samples,
            mean_expansion_order: mean,
            energy,
            energy_standard_error: variance.sqrt() / (count.sqrt() * beta),
            energy_per_site: energy / num_sites as f64,
            heat_capacity,
            heat_capacity_per_site: heat_capacity / num_sites as f64,
        })
    }
}
/// Thermodynamic results derived from expansion-order moments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermodynamicResults {
    /// Number of samples.
    pub samples: u64,
    /// Mean expansion order.
    pub mean_expansion_order: f64,
    /// Total physical energy.
    pub energy: f64,
    /// Naive independent-sample standard error of the energy estimator.
    ///
    /// This does not correct for autocorrelation; independent logical chains
    /// should be combined when a confidence interval is required.
    pub energy_standard_error: f64,
    /// Energy per site.
    pub energy_per_site: f64,
    /// Heat capacity estimator.
    pub heat_capacity: f64,
    /// Heat capacity per site.
    pub heat_capacity_per_site: f64,
}
