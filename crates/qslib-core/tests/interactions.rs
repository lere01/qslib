use qslib_core::{
    Bond, BondMultiplicity, Boundary, DenseCouplings, InteractionChannel, InteractionError,
    InteractionIdentity, InteractionTable, RectangularGeometry, SiteCount, SiteId, SparseCouplings,
};

fn bond(i: u32, j: u32) -> Bond {
    Bond::new(SiteId::new(i), SiteId::new(j)).unwrap()
}

#[test]
fn weighted_interactions_keep_pair_dependent_signed_and_zero_coefficients() {
    let table = InteractionTable::new(
        SiteCount::new(3).unwrap(),
        vec![
            (bond(2, 0), InteractionChannel::HeisenbergExchange, -1.25),
            (bond(0, 1), InteractionChannel::HeisenbergExchange, 0.0),
        ],
    )
    .unwrap();
    assert_eq!(table.interactions()[0].bond().first(), SiteId::new(0));
    assert_eq!(table.interactions()[0].coefficient(), 0.0);
    assert_eq!(table.interactions()[1].coefficient(), -1.25);
    assert_eq!(table.active_interactions().len(), 1);
    assert_eq!(table.declared_interactions().len(), 2);
}

#[test]
fn named_terms_and_channel_validation_preserve_full_identity() {
    let named_a =
        InteractionIdentity::named(bond(0, 1), InteractionChannel::HeisenbergExchange, "j1")
            .unwrap();
    let named_b =
        InteractionIdentity::named(bond(0, 1), InteractionChannel::HeisenbergExchange, "j2")
            .unwrap();
    let table = InteractionTable::new_with_identities(
        SiteCount::new(2).unwrap(),
        vec![(named_a, 1.0), (named_b, -1.0)],
    )
    .unwrap();
    assert_eq!(table.interactions().len(), 2);
    assert!(InteractionChannel::generic(" ").is_err());
    assert!(InteractionChannel::generic("bad name").is_err());
}

#[test]
fn duplicate_identity_is_rejected_but_distinct_channels_are_allowed() {
    let duplicate = InteractionTable::new(
        SiteCount::new(2).unwrap(),
        vec![
            (bond(0, 1), InteractionChannel::Generic("x".into()), 1.0),
            (bond(1, 0), InteractionChannel::Generic("x".into()), 2.0),
        ],
    );
    assert!(matches!(
        duplicate,
        Err(InteractionError::DuplicateIdentity { .. })
    ));

    let channels = InteractionTable::new(
        SiteCount::new(2).unwrap(),
        vec![
            (bond(0, 1), InteractionChannel::Generic("x".into()), 1.0),
            (bond(0, 1), InteractionChannel::Generic("z".into()), 2.0),
        ],
    )
    .unwrap();
    assert_eq!(channels.interactions().len(), 2);
}

#[test]
fn dense_and_sparse_couplings_validate_symmetry_and_never_double_count() {
    let dense = DenseCouplings::new(
        SiteCount::new(3).unwrap(),
        vec![0.0, 2.0, -1.0, 2.0, 0.0, 0.5, -1.0, 0.5, 0.0],
    )
    .unwrap();
    let interactions = dense
        .to_interactions(InteractionChannel::HeisenbergExchange)
        .unwrap();
    assert_eq!(interactions.len(), 3);
    assert_eq!(interactions[0].coefficient(), 2.0);
    assert!(
        DenseCouplings::new_with_tolerance(
            SiteCount::new(2).unwrap(),
            vec![0.0, 1.0, 2.0, 0.0],
            1.0e-12
        )
        .is_err()
    );
    assert!(
        DenseCouplings::new_with_tolerance(
            SiteCount::new(2).unwrap(),
            vec![0.0, 1.0, 1.0 + 1.0e-13, 0.0],
            1.0e-12
        )
        .is_ok()
    );

    let sparse = SparseCouplings::new(
        SiteCount::new(3).unwrap(),
        vec![
            (SiteId::new(1), SiteId::new(0), 2.0),
            (SiteId::new(2), SiteId::new(0), -1.0),
        ],
    )
    .unwrap();
    assert_eq!(
        sparse
            .to_interactions(InteractionChannel::HeisenbergExchange)
            .unwrap()
            .len(),
        2
    );
    assert!(
        SparseCouplings::new(
            SiteCount::new(2).unwrap(),
            vec![(SiteId::new(0), SiteId::new(0), 1.0)]
        )
        .is_err()
    );
}

#[test]
fn disorder_realization_is_order_and_schedule_independent() {
    let table = InteractionTable::new(
        SiteCount::new(3).unwrap(),
        vec![
            (bond(0, 2), InteractionChannel::HeisenbergExchange, 0.0),
            (bond(0, 1), InteractionChannel::HeisenbergExchange, 0.0),
        ],
    )
    .unwrap();
    let a = table.realize_uniform_disorder([7; 32], -1.0, 1.0).unwrap();
    let b = InteractionTable::new(
        SiteCount::new(3).unwrap(),
        vec![
            (bond(0, 1), InteractionChannel::HeisenbergExchange, 0.0),
            (bond(2, 0), InteractionChannel::HeisenbergExchange, 0.0),
        ],
    )
    .unwrap()
    .realize_uniform_disorder([7; 32], -1.0, 1.0)
    .unwrap();
    assert_eq!(a, b);
    assert_eq!(a.provenance().seed(), [7; 32]);
    assert_eq!(a.provenance().rng_algorithm(), "chacha20");
    assert_eq!(a.provenance().seed_scheme(), "qslib-seed-v1");
    assert_eq!(a.provenance().domain(), "disorder");
    assert_eq!(a.provenance().semantics(), "replacement");
    assert_ne!(
        a,
        table
            .realize_uniform_disorder_at([7; 32], 1, -1.0, 1.0)
            .unwrap()
    );
    assert!(
        a.interactions()
            .iter()
            .all(|term| (-1.0..=1.0).contains(&term.coefficient()))
    );
}

#[test]
fn disorder_seed_scheme_has_a_pinned_reference_vector() {
    let table = InteractionTable::new(
        SiteCount::new(2).unwrap(),
        vec![(bond(0, 1), InteractionChannel::HeisenbergExchange, 0.0)],
    )
    .unwrap();
    let realization = table
        .realize_uniform_disorder_at([7; 32], 3, -1.0, 1.0)
        .unwrap();
    assert_eq!(
        realization.interactions()[0].coefficient(),
        -0.579_846_429_556_808_4
    );
}

#[test]
fn periodic_image_identity_changes_the_disorder_stream() {
    let geometry = RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Periodic).unwrap();
    let images = geometry.bonds(BondMultiplicity::PeriodicImages).unwrap();
    let first = images[0];
    let second = images
        .iter()
        .copied()
        .find(|bond| {
            bond.first() == first.first()
                && bond.second() == first.second()
                && bond.source() != first.source()
        })
        .unwrap();
    let table = InteractionTable::new(
        SiteCount::new(4).unwrap(),
        vec![
            (first, InteractionChannel::HeisenbergExchange, 0.0),
            (second, InteractionChannel::HeisenbergExchange, 0.0),
        ],
    )
    .unwrap();
    let realization = table.realize_uniform_disorder([7; 32], -1.0, 1.0).unwrap();
    assert_ne!(
        realization.interactions()[0].coefficient(),
        realization.interactions()[1].coefficient()
    );
}
