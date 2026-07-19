use qslib_core::{BondMultiplicity, Boundary, RectangularGeometry, SiteId, XMajorAdapter};

#[test]
fn modern_row_major_and_sse_simple_geometry_have_the_same_resolved_pairs() {
    let geometry = RectangularGeometry::new(3, 2, Boundary::Open, Boundary::Open).unwrap();
    let qslib_pairs: Vec<_> = geometry
        .bonds(BondMultiplicity::Simple)
        .unwrap()
        .iter()
        .map(|bond| (bond.first().get(), bond.second().get()))
        .collect();
    // Neutral expected data mirrors the independently specified SSE nearest-neighbour list.
    let sse_pairs = vec![(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)];
    assert_eq!(qslib_pairs, sse_pairs);
}

#[test]
fn legacy_x_major_pairs_fail_without_and_pass_with_named_adapter() {
    let adapter = XMajorAdapter::new(3, 2).unwrap();
    let legacy_pair = (SiteId::new(1), SiteId::new(3));
    let canonical_pair = (
        adapter.to_canonical(legacy_pair.0).unwrap(),
        adapter.to_canonical(legacy_pair.1).unwrap(),
    );
    assert_ne!(legacy_pair, (canonical_pair.0, canonical_pair.1));
    assert_eq!((canonical_pair.0.get(), canonical_pair.1.get()), (3, 4));
}
