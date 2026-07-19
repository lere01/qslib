//! Exercise the converted four-site standalone SSE parity configuration.

use qslib_sse::{
    BasisSseState, LegacyModelKind, LegacySpin, LocalSseModel, Operator, SimulationConfig,
    convert_legacy_bits, run_parallel_chains,
};

fn main() {
    let model = LocalSseModel::tfim(4, &[(0, 1), (1, 2), (2, 3), (3, 0)], 1.0, 0.5)
        .expect("converted four-site TFIM");
    let state = BasisSseState::new(
        convert_legacy_bits(
            LegacyModelKind::Tfim,
            &[
                LegacySpin::Up,
                LegacySpin::Down,
                LegacySpin::Up,
                LegacySpin::Down,
            ],
        ),
        vec![Operator::identity(); 128],
    )
    .expect("valid state");
    let results = run_parallel_chains(
        model,
        state,
        4.0,
        SimulationConfig {
            thermalization_sweeps: 1_000,
            measurement_sweeps: 10_000,
            sweeps_per_measurement: 1,
        },
        24_301,
        4,
        4,
    )
    .expect("parity-path run");
    let energy_per_site = results
        .iter()
        .map(|result| result.thermodynamics.energy_per_site)
        .sum::<f64>()
        / results.len() as f64;
    println!("qslib energy/site = {energy_per_site}; standalone reference = -1.0630828125");
}
