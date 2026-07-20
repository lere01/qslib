use qslib_core::{
    BasisBit, BasisError, BasisState, FullBasis, PhysicalAxis, SectorBasis, SimulationBasis,
    SiteCount, SiteId, WordWidth,
};
use rand_chacha::ChaCha8Rng;
use rand_core::{Rng, SeedableRng};
use std::str::FromStr;

#[test]
fn checked_identifiers_and_axes_reject_ambiguous_inputs() {
    assert!(matches!(SiteCount::new(0), Err(BasisError::EmptySystem)));
    let sites = SiteCount::new(5).expect("positive site count");
    assert_eq!(sites.get(), 5);
    assert_eq!(SiteId::try_from_usize(4).expect("u32 site id").get(), 4);
    if let Some(value) = (u32::MAX as usize).checked_add(1) {
        assert!(matches!(
            SiteId::try_from_usize(value),
            Err(BasisError::IdentifierOverflow { .. })
        ));
    }
    assert!(matches!(
        sites.validate(SiteId::new(5)),
        Err(BasisError::SiteOutOfRange { .. })
    ));
    assert_eq!(PhysicalAxis::X.as_str(), "x");
    assert_eq!(SimulationBasis::Z.as_str(), "z");
    assert_ne!(PhysicalAxis::X, PhysicalAxis::Y);
    assert_eq!(PhysicalAxis::from_str("y").expect("axis"), PhysicalAxis::Y);
    assert_eq!(
        SimulationBasis::from_str("z").expect("basis"),
        SimulationBasis::Z
    );
    assert!(PhysicalAxis::from_str("Z").is_err());
}

#[test]
fn dense_bits_have_explicit_zero_one_meaning_and_checked_values() {
    let state = BasisState::from_raw_bits(&[1, 0, 1, 1]).expect("binary state");
    assert_eq!(state.len(), 4);
    assert_eq!(
        state.bits(),
        &[BasisBit::One, BasisBit::Zero, BasisBit::One, BasisBit::One]
    );
    assert_eq!(state.hamming_weight(), 3);
    assert_eq!(state.pauli_eigenvalues(), &[-1, 1, -1, -1]);
    let view = state.as_view();
    assert_eq!(view.len(), 4);
    assert_eq!(view.hamming_weight(), 3);
    assert_eq!(view.bits(), state.bits());
    assert!(matches!(
        BasisState::from_raw_bits(&[0, 2]),
        Err(BasisError::InvalidBit { .. })
    ));
}

#[test]
fn little_endian_packed_states_round_trip_and_track_word_width() {
    let state = BasisState::from_raw_bits(&[1, 0, 1, 1]).expect("binary state");
    let packed = state.pack().expect("packed state");
    assert_eq!(packed.words_le(), &[13]);
    assert_eq!(packed.hamming_weight(), 3);
    assert_eq!(packed.bit(0).expect("site 0"), BasisBit::One);
    assert_eq!(packed.bit(1).expect("site 1"), BasisBit::Zero);
    for width in [
        WordWidth::U8,
        WordWidth::U16,
        WordWidth::U32,
        WordWidth::U64,
    ] {
        let bytes = packed.to_bytes(width).expect("serialized");
        let restored =
            qslib_core::PackedState::from_bytes(4, width, &bytes).expect("packed state round trip");
        assert_eq!(restored, packed);
    }

    let mut wide = vec![0; 65];
    wide[64] = 1;
    let wide_state = BasisState::from_raw_bits(&wide).expect("wide binary state");
    let wide_packed = wide_state.pack().expect("wide packed state");
    assert_eq!(wide_packed.words_le(), &[0, 1]);
    let bytes = wide_packed
        .to_bytes(WordWidth::U64)
        .expect("wide serialization");
    assert_eq!(bytes.len(), 16);
    assert_eq!(
        qslib_core::PackedState::from_bytes(65, WordWidth::U64, &bytes).expect("wide round trip"),
        wide_packed
    );
}

#[test]
fn full_basis_is_canonical_little_endian_integer_order() {
    let basis = FullBasis::new(SiteCount::new(3).expect("sites")).expect("full basis");
    let masks: Vec<_> = basis.map(|state| state.words_le()[0]).collect();
    assert_eq!(masks, (0_u64..8).collect::<Vec<_>>());
    assert!(matches!(
        FullBasis::new(SiteCount::new(usize::BITS as usize).expect("boundary sites")),
        Err(BasisError::DimensionOverflow { .. })
    ));
}

#[test]
fn reference_scalars_reject_non_finite_values() {
    let value: qslib_core::Real = 1.25;
    assert_eq!(qslib_core::ensure_finite(value).expect("finite"), value);
    assert!(matches!(
        qslib_core::ensure_finite(f64::NAN),
        Err(BasisError::NonFiniteScalar { .. })
    ));
    let complex: qslib_core::Complex64 = qslib_core::Complex64::new(1.0, -0.5);
    assert_eq!(complex.re, 1.0);
    assert_eq!(complex.im, -0.5);
}

#[test]
fn fixed_weight_basis_is_canonical_and_checked() {
    let basis = SectorBasis::new(SiteCount::new(4).expect("sites"), 2).expect("sector");
    let masks: Vec<_> = basis.map(|state| state.words_le()[0]).collect();
    assert_eq!(masks, vec![3, 5, 6, 9, 10, 12]);
    assert!(matches!(
        SectorBasis::new(SiteCount::new(3).expect("sites"), 4),
        Err(BasisError::WeightOutOfRange { .. })
    ));
    assert_eq!(
        SectorBasis::new(SiteCount::new(4).expect("sites"), 0)
            .expect("zero sector")
            .map(|state| state.words_le()[0])
            .collect::<Vec<_>>(),
        vec![0]
    );
    assert_eq!(
        SectorBasis::new(SiteCount::new(4).expect("sites"), 4)
            .expect("full sector")
            .map(|state| state.words_le()[0])
            .collect::<Vec<_>>(),
        vec![15]
    );
}

#[test]
fn packed_state_rejects_noncanonical_high_bits_and_wrong_serialized_width() {
    assert!(matches!(
        qslib_core::PackedState::from_words(65, &[0, 2]),
        Err(BasisError::NonCanonicalHighBits { .. })
    ));
    assert!(matches!(
        qslib_core::PackedState::from_bytes(9, WordWidth::U8, &[0]),
        Err(BasisError::SerializedLength { .. })
    ));
    assert!(matches!(
        qslib_core::PackedState::from_bytes(9, WordWidth::U8, &[0, 2]),
        Err(BasisError::NonCanonicalHighBits { .. })
    ));
}

#[test]
fn generated_packed_state_properties_are_bounded_and_reproducible() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x5153_5441_5445_0001);
    let widths = [
        WordWidth::U8,
        WordWidth::U16,
        WordWidth::U32,
        WordWidth::U64,
    ];

    let cases = if cfg!(miri) { 16 } else { 512 };
    for _case in 0..cases {
        let site_count: usize = 1 + (rng.next_u32() as usize % 257);
        let raw_bits = (0..site_count)
            .map(|_| (rng.next_u32() & 1) as u8)
            .collect::<Vec<_>>();
        let dense = BasisState::from_raw_bits(&raw_bits).expect("generated binary state");
        let packed = dense.pack().expect("generated packed state");

        assert_eq!(packed.site_count(), site_count);
        assert_eq!(packed.hamming_weight(), dense.hamming_weight());
        for (site, expected) in raw_bits.iter().copied().enumerate() {
            assert_eq!(packed.bit(site).expect("generated site").as_u8(), expected);
        }

        for width in widths {
            let bytes = packed.to_bytes(width).expect("generated serialization");
            let restored = qslib_core::PackedState::from_bytes(site_count, width, &bytes)
                .expect("generated deserialization");
            assert_eq!(restored, packed);
            assert_eq!(
                restored.to_bytes(width).expect("stable serialization"),
                bytes
            );
        }
    }
}

#[test]
fn bounded_state_conversion_fuzz_inputs_never_panic() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x5153_5441_5445_0002);
    let widths = [
        WordWidth::U8,
        WordWidth::U16,
        WordWidth::U32,
        WordWidth::U64,
    ];

    let cases = if cfg!(miri) { 8 } else { 1_024 };
    for _case in 0..cases {
        let site_count: usize = 1 + (rng.next_u32() as usize % 257);
        for width in widths {
            let serialized_words = site_count.div_ceil(width.bits());
            let byte_count = serialized_words * width.bytes();
            let mut bytes = vec![0_u8; byte_count];
            for byte in &mut bytes {
                *byte = rng.next_u32() as u8;
            }

            if let Ok(state) = qslib_core::PackedState::from_bytes(site_count, width, &bytes) {
                let canonical = state.to_bytes(width).expect("accepted state serializes");
                let restored = qslib_core::PackedState::from_bytes(site_count, width, &canonical)
                    .expect("canonical state deserializes");
                assert_eq!(restored, state);
            }
        }
    }
}
