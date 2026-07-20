use qslib_core::{
    BasisState, Boundary, InteractionChannel, InteractionIdentity, InteractionTable, Pauli,
    PauliString, Permutation, RectangularGeometry, SimulationBasis, SiteCount, SiteId,
    SpinInversion, SymmetryCharacter,
};
use rand_chacha::ChaCha8Rng;
use rand_core::{Rng, SeedableRng};

#[test]
fn gather_permutations_compose_invert_and_act_on_bits() {
    let permutation = Permutation::new(
        SiteCount::new(3).unwrap(),
        vec![SiteId::new(1), SiteId::new(2), SiteId::new(0)],
    )
    .unwrap();
    let state = BasisState::from_raw_bits(&[0, 1, 1]).unwrap();
    assert_eq!(
        permutation.apply_state(&state).unwrap(),
        BasisState::from_raw_bits(&[1, 1, 0]).unwrap()
    );
    assert_eq!(
        permutation
            .compose(&permutation.inverse().unwrap())
            .unwrap(),
        Permutation::identity(SiteCount::new(3).unwrap()).unwrap()
    );
    assert!(
        Permutation::new(
            SiteCount::new(3).unwrap(),
            vec![SiteId::new(0), SiteId::new(0), SiteId::new(2)]
        )
        .is_err()
    );
}

#[test]
fn translations_and_point_groups_have_group_laws() {
    let geometry = RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Periodic).unwrap();
    let tx = qslib_core::translation(&geometry, 1, 0).unwrap();
    let ty = qslib_core::translation(&geometry, 0, 1).unwrap();
    assert_eq!(
        tx.apply_state(&BasisState::from_raw_bits(&[1, 0, 0, 0]).unwrap())
            .unwrap(),
        BasisState::from_raw_bits(&[0, 1, 0, 0]).unwrap()
    );
    assert_eq!(tx.compose(&ty).unwrap(), ty.compose(&tx).unwrap());
    assert_eq!(qslib_core::translation_group(&geometry).unwrap().order(), 4);
    assert_eq!(
        qslib_core::square_point_group(&geometry).unwrap().order(),
        8
    );
    assert_eq!(
        qslib_core::rectangle_point_group(&geometry)
            .unwrap()
            .order(),
        4
    );
    assert!(
        qslib_core::square_point_group(
            &RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Open).unwrap()
        )
        .is_err()
    );
    assert!(
        qslib_core::square_point_group(
            &RectangularGeometry::with_kind(
                2,
                2,
                Boundary::Periodic,
                Boundary::Periodic,
                qslib_core::LatticeKind::Triangular,
            )
            .unwrap()
        )
        .is_err()
    );
    assert!(qslib_core::translation(&geometry, isize::MIN, 0).is_err());
}

#[test]
fn spin_inversion_and_orbits_are_explicit() {
    let inversion = SpinInversion::new(SiteCount::new(3).unwrap());
    let state = BasisState::from_raw_bits(&[0, 1, 0]).unwrap();
    assert_eq!(
        inversion.apply(&state).unwrap(),
        BasisState::from_raw_bits(&[1, 0, 1]).unwrap()
    );
    assert!(
        inversion
            .apply(&BasisState::from_raw_bits(&[0, 1]).unwrap())
            .is_err()
    );
    let translation_group = qslib_core::translation_group(
        &RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Periodic).unwrap(),
    )
    .unwrap();
    let orbit = translation_group
        .orbit(&BasisState::from_raw_bits(&[1, 0, 0, 0]).unwrap())
        .unwrap();
    assert_eq!(orbit.len(), 4);
    let expected = orbit
        .iter()
        .min_by_key(|state| {
            state
                .bits()
                .iter()
                .enumerate()
                .map(|(site, bit)| (bit.as_u8() as usize) << site)
                .sum::<usize>()
        })
        .unwrap()
        .clone();
    assert_eq!(
        translation_group
            .canonical_representative(&orbit[0])
            .unwrap(),
        expected
    );
}

#[test]
fn orbit_ordering_supports_states_wider_than_one_machine_word() {
    let sites = SiteCount::new(65).unwrap();
    let mut swap_sources = (0..65).map(SiteId::new).collect::<Vec<_>>();
    swap_sources[0] = SiteId::new(64);
    swap_sources[64] = SiteId::new(0);
    let group = qslib_core::FiniteGroup::new(vec![
        Permutation::identity(sites).unwrap(),
        Permutation::new(sites, swap_sources).unwrap(),
    ])
    .unwrap();
    let mut bits = vec![0_u8; 65];
    bits[0] = 1;
    let state = BasisState::from_raw_bits(&bits).unwrap();
    let representative = group.canonical_representative(&state).unwrap();
    assert_eq!(representative.bits()[0].as_u8(), 1);
}

#[test]
fn trivial_character_projection_is_idempotent_and_normalized() {
    let group = qslib_core::translation_group(
        &RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Periodic).unwrap(),
    )
    .unwrap();
    let character = SymmetryCharacter::trivial(group.order());
    let amplitudes = vec![qslib_core::Complex64::new(1.0, 0.0); 16];
    let projected = group.project_amplitudes(&amplitudes, &character).unwrap();
    assert!(
        projected
            .iter()
            .all(|value| *value == qslib_core::Complex64::new(1.0, 0.0))
    );
    let nonuniform = (0..16)
        .map(|value| qslib_core::Complex64::new(value as f64, 0.0))
        .collect::<Vec<_>>();
    let once = group.project_amplitudes(&nonuniform, &character).unwrap();
    let twice = group.project_amplitudes(&once, &character).unwrap();
    assert_eq!(once, twice);
}

#[test]
fn group_aligned_characters_validate_homomorphism_and_project_nontrivially() {
    let swap = Permutation::new(
        SiteCount::new(2).unwrap(),
        vec![SiteId::new(1), SiteId::new(0)],
    )
    .unwrap();
    let group = qslib_core::FiniteGroup::new(vec![
        Permutation::identity(SiteCount::new(2).unwrap()).unwrap(),
        swap,
    ])
    .unwrap();
    let character = SymmetryCharacter::new_for_group(
        &group,
        vec![
            qslib_core::Complex64::new(1.0, 0.0),
            qslib_core::Complex64::new(-1.0, 0.0),
        ],
    )
    .unwrap();
    let amplitudes = vec![
        qslib_core::Complex64::new(0.0, 0.0),
        qslib_core::Complex64::new(1.0, 0.0),
        qslib_core::Complex64::new(0.0, 0.0),
        qslib_core::Complex64::new(0.0, 0.0),
    ];
    let projected = group.project_amplitudes(&amplitudes, &character).unwrap();
    assert_eq!(projected[1], qslib_core::Complex64::new(0.5, 0.0));
    assert_eq!(projected[2], qslib_core::Complex64::new(-0.5, 0.0));
    assert!(
        SymmetryCharacter::new_for_group(
            &group,
            vec![
                qslib_core::Complex64::new(1.0, 0.0),
                qslib_core::Complex64::new(0.0, 1.0),
            ],
        )
        .is_err()
    );
    let raw = SymmetryCharacter::new(vec![
        qslib_core::Complex64::new(1.0, 0.0),
        qslib_core::Complex64::new(0.0, 1.0),
    ])
    .unwrap();
    assert!(group.project_amplitudes(&amplitudes, &raw).is_err());
    let group_three = qslib_core::FiniteGroup::new(vec![
        Permutation::identity(SiteCount::new(3).unwrap()).unwrap(),
        Permutation::new(
            SiteCount::new(3).unwrap(),
            vec![SiteId::new(1), SiteId::new(0), SiteId::new(2)],
        )
        .unwrap(),
    ])
    .unwrap();
    let other_group = qslib_core::FiniteGroup::new(vec![
        Permutation::identity(SiteCount::new(3).unwrap()).unwrap(),
        Permutation::new(
            SiteCount::new(3).unwrap(),
            vec![SiteId::new(0), SiteId::new(2), SiteId::new(1)],
        )
        .unwrap(),
    ])
    .unwrap();
    let other_character = SymmetryCharacter::new_for_group(
        &other_group,
        vec![
            qslib_core::Complex64::new(1.0, 0.0),
            qslib_core::Complex64::new(-1.0, 0.0),
        ],
    )
    .unwrap();
    let amplitudes_three = vec![qslib_core::Complex64::new(0.0, 0.0); 8];
    assert!(
        group_three
            .project_amplitudes(&amplitudes_three, &other_character)
            .is_err()
    );
    assert!(SymmetryCharacter::new(vec![qslib_core::Complex64::new(f64::NAN, 0.0)]).is_err());
}

#[test]
fn interaction_symmetry_rejects_broken_bonds_and_accepts_translation_invariant_terms() {
    let geometry = RectangularGeometry::new(3, 1, Boundary::Periodic, Boundary::Periodic).unwrap();
    let bonds = geometry
        .bonds(qslib_core::BondMultiplicity::Simple)
        .unwrap();
    let table = InteractionTable::new_with_identities(
        SiteCount::new(3).unwrap(),
        bonds
            .iter()
            .copied()
            .map(|bond| {
                (
                    InteractionIdentity::new(bond, InteractionChannel::IsingZZ),
                    1.0,
                )
            })
            .collect(),
    )
    .unwrap();
    let translation = qslib_core::translation(&geometry, 1, 0).unwrap();
    assert!(qslib_core::is_interaction_symmetry(&table, &translation));
    let bond = bonds[0];
    let broken = InteractionTable::new(
        SiteCount::new(3).unwrap(),
        vec![(bond, InteractionChannel::IsingZZ, 2.0)],
    )
    .unwrap();
    assert!(!qslib_core::is_interaction_symmetry(&broken, &translation));
}

#[test]
fn full_model_validation_checks_onsite_fields_and_spin_inversion() {
    let geometry = RectangularGeometry::new(3, 1, Boundary::Periodic, Boundary::Periodic).unwrap();
    let bonds = geometry
        .bonds(qslib_core::BondMultiplicity::Simple)
        .unwrap();
    let table = InteractionTable::new(
        SiteCount::new(3).unwrap(),
        bonds
            .iter()
            .copied()
            .map(|bond| (bond, InteractionChannel::IsingZZ, 1.0))
            .collect(),
    )
    .unwrap();
    let translation = qslib_core::translation(&geometry, 1, 0).unwrap();
    let uniform =
        qslib_core::tfim(&table, &[0.5, 0.5, 0.5], qslib_core::SimulationBasis::Z).unwrap();
    assert!(qslib_core::validate_model_symmetry(&uniform, &translation).is_ok());
    assert!(qslib_core::validate_spin_inversion(&uniform).is_ok());
    let nonuniform =
        qslib_core::tfim(&table, &[0.5, 0.6, 0.5], qslib_core::SimulationBasis::Z).unwrap();
    let error = qslib_core::validate_model_symmetry(&nonuniform, &translation).unwrap_err();
    assert!(error.to_string().contains("mapped operator"));
}

#[test]
fn pauli_support_can_be_mapped_by_a_permutation() {
    let permutation = Permutation::new(
        SiteCount::new(2).unwrap(),
        vec![SiteId::new(1), SiteId::new(0)],
    )
    .unwrap();
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::X)]).unwrap();
    assert_eq!(
        permutation.map_pauli_string(&operator).unwrap().factors(),
        &[(SiteId::new(1), Pauli::X)]
    );
    let _ = SimulationBasis::Z;
}

#[test]
fn diagonal_and_sublattice_gauges_apply_phases_with_sign_only_classification() {
    let state = BasisState::from_raw_bits(&[1, 1]).unwrap();
    let gauge = qslib_core::DiagonalGauge::new(vec![0.5, 1.0]).unwrap();
    let phase = gauge.phase(&state).unwrap();
    assert!((phase.re + 0.0).abs() < 1.0e-12);
    assert!((phase.im + 1.0).abs() < 1.0e-12);
    assert!(!gauge.is_sign_only());
    let checkerboard = qslib_core::sublattice_gauge(
        &RectangularGeometry::new(2, 2, Boundary::Open, Boundary::Open).unwrap(),
    )
    .unwrap();
    assert!(checkerboard.is_sign_only());
    let checker_state = BasisState::from_raw_bits(&[1, 0, 0, 0]).unwrap();
    assert_eq!(
        checkerboard
            .apply_amplitude(&checker_state, qslib_core::Complex64::new(1.0, 0.0))
            .unwrap(),
        qslib_core::Complex64::new(1.0, 0.0)
    );
    let odd_checker_state = BasisState::from_raw_bits(&[0, 1, 0, 0]).unwrap();
    assert_eq!(
        checkerboard.phase(&odd_checker_state).unwrap(),
        qslib_core::Complex64::new(-1.0, 0.0)
    );
    let mixed = qslib_core::DiagonalGauge::new(vec![1.0e15 + 1.0, 0.5]).unwrap();
    assert!(
        (mixed
            .phase(&BasisState::from_raw_bits(&[1, 1]).unwrap())
            .unwrap()
            .im
            + 1.0)
            .abs()
            < 1.0e-12
    );
    assert!(qslib_core::DiagonalGauge::new(vec![f64::NAN]).is_err());
}

#[test]
fn legacy_x_major_symmetry_adapter_is_explicit() {
    let adapter = qslib_core::legacy_x_major_permutation(2, 3).unwrap();
    let legacy_order_state = BasisState::from_raw_bits(&[0, 1, 0, 0, 0, 0]).unwrap();
    let canonical = adapter.apply_state(&legacy_order_state).unwrap();
    assert_eq!(
        canonical,
        BasisState::from_raw_bits(&[0, 0, 1, 0, 0, 0]).unwrap()
    );
}

#[test]
fn generated_permutations_preserve_inverse_and_composition() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x5359_4d4d_4554_0001);

    let cases = if cfg!(miri) { 32 } else { 256 };
    for _case in 0..cases {
        let site_count: usize = 1 + (rng.next_u32() as usize % 64);
        let sites = SiteCount::new(site_count).expect("generated site count");
        let mut source_indices = (0..site_count).collect::<Vec<_>>();
        for index in (1..site_count).rev() {
            let swap = (rng.next_u32() as usize) % (index + 1);
            source_indices.swap(index, swap);
        }
        let permutation = Permutation::new(
            sites,
            source_indices
                .iter()
                .copied()
                .map(|index| SiteId::new(index as u32))
                .collect(),
        )
        .expect("generated permutation");
        let inverse = permutation.inverse().expect("generated inverse");
        let identity = Permutation::identity(sites).expect("generated identity");

        assert_eq!(
            permutation.compose(&inverse).expect("right identity"),
            identity
        );
        assert_eq!(
            inverse.compose(&permutation).expect("left identity"),
            identity
        );

        let state = BasisState::from_raw_bits(
            &(0..site_count)
                .map(|_| (rng.next_u32() & 1) as u8)
                .collect::<Vec<_>>(),
        )
        .expect("generated state");
        let transformed = permutation.apply_state(&state).expect("action");
        assert_eq!(
            inverse.apply_state(&transformed).expect("inverse action"),
            state
        );
    }
}
