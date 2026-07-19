//! Versioned, domain-separated random-stream derivation.

use blake3::Hasher;

/// Stable seed-scheme identifier persisted in qslib artifacts.
pub const QSLIB_SEED_SCHEME: &str = "qslib-seed-v1";

/// Expand a convenient `u64` into the canonical 32-byte master seed.
pub fn expand_master_seed(seed: u64) -> [u8; 32] {
    blake3::derive_key("qslib master seed v1", &seed.to_le_bytes())
}

/// Derive a 32-byte child stream using the ADR-0003 canonical framing.
pub fn derive_seed(master_seed: &[u8; 32], domain: &str, indices: &[u64]) -> [u8; 32] {
    let mut hasher = Hasher::new_keyed(master_seed);
    hasher.update(QSLIB_SEED_SCHEME.as_bytes());
    hasher.update(&[0]);
    hasher.update(&(domain.len() as u32).to_le_bytes());
    hasher.update(domain.as_bytes());
    hasher.update(&(indices.len() as u32).to_le_bytes());
    for index in indices {
        hasher.update(&index.to_le_bytes());
    }
    *hasher.finalize().as_bytes()
}
