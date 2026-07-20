#![no_main]

use libfuzzer_sys::fuzz_target;
use qslib_core::{PackedState, WordWidth};

const WIDTHS: [WordWidth; 4] = [WordWidth::U8, WordWidth::U16, WordWidth::U32, WordWidth::U64];

fuzz_target!(|input: &[u8]| {
    let site_count = 1 + usize::from(input.first().copied().unwrap_or(0)) % 257;
    let payload = input.get(1..).unwrap_or_default();

    for width in WIDTHS {
        let serialized_words = site_count.div_ceil(width.bits());
        let byte_count = serialized_words * width.bytes();
        let mut bytes = vec![0_u8; byte_count];
        for (index, byte) in bytes.iter_mut().enumerate() {
            *byte = payload.get(index).copied().unwrap_or(0);
        }

        if let Ok(state) = PackedState::from_bytes(site_count, width, &bytes) {
            let canonical = state
                .to_bytes(width)
                .expect("accepted packed state must serialize");
            let restored = PackedState::from_bytes(site_count, width, &canonical)
                .expect("canonical packed state must deserialize");
            assert_eq!(restored, state);
        }
    }
});
