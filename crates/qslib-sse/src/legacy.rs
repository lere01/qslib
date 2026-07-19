//! Explicit adapters for legacy spin labels and deterministic chain seeds.

use qslib_core::{BasisBit, derive_seed, expand_master_seed};

/// Historical two-state label used by older input formats.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LegacySpin {
    /// Conventional spin-up label.
    Up,
    /// Conventional spin-down label.
    Down,
}

/// Model family whose physical meaning determines the legacy-to-canonical map.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LegacyModelKind {
    /// TFIM convention: up is the positive Pauli-Z eigenstate.
    Tfim,
    /// Rydberg convention: up denotes an occupied excitation.
    Rydberg,
}

/// Convert legacy labels into qslib's canonical bit basis.
pub fn convert_legacy_bits(kind: LegacyModelKind, spins: &[LegacySpin]) -> Vec<BasisBit> {
    spins
        .iter()
        .map(|spin| match (kind, spin) {
            (LegacyModelKind::Tfim, LegacySpin::Up)
            | (LegacyModelKind::Rydberg, LegacySpin::Down) => BasisBit::Zero,
            (LegacyModelKind::Tfim, LegacySpin::Down)
            | (LegacyModelKind::Rydberg, LegacySpin::Up) => BasisBit::One,
        })
        .collect()
}

/// Derive an order-independent canonical 32-byte SSE chain seed.
pub fn derive_chain_seed(master_seed: u64, chain_index: u64) -> [u8; 32] {
    let master = expand_master_seed(master_seed);
    derive_seed(&master, "sse_chain", &[chain_index])
}

/// Reproduce the standalone SSE SplitMix-style seed for legacy artifact parity.
///
/// New qslib chains must use [`derive_chain_seed`]; this function exists only
/// when comparing historical standalone runs.
pub fn derive_legacy_chain_seed(master_seed: u64, chain_index: u64) -> u64 {
    let mut value = master_seed ^ chain_index.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

/// Return deterministic seeds for a logical chain range.
pub fn logical_chain_seeds(master_seed: u64, chain_count: usize) -> Vec<[u8; 32]> {
    (0..chain_count)
        .map(|index| derive_chain_seed(master_seed, index as u64))
        .collect()
}
