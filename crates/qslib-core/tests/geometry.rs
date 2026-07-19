use qslib_core::{
    BondMultiplicity, Boundary, CustomGeometry, LatticeKind, RectangularGeometry, ShellTolerance,
    SiteId, XMajorAdapter,
};

#[test]
fn rectangular_geometry_is_row_major_and_supports_mixed_boundaries() {
    let geometry = RectangularGeometry::new(3, 2, Boundary::Open, Boundary::Periodic).unwrap();
    assert_eq!(geometry.site_count().get(), 6);
    assert_eq!(geometry.site_id(2, 1).unwrap(), SiteId::new(5));
    assert_eq!(
        geometry.coordinate(SiteId::new(4)).unwrap().as_tuple(),
        (1.0, 1.0)
    );

    let bonds = geometry.bonds(BondMultiplicity::Simple).unwrap();
    assert_eq!(bonds.len(), 7);
    assert!(
        bonds
            .iter()
            .all(|bond| { *bond == qslib_core::Bond::new(bond.first(), bond.second()).unwrap() })
    );
    assert_eq!(
        bonds
            .iter()
            .map(|bond| (bond.first().get(), bond.second().get()))
            .collect::<Vec<_>>(),
        vec![(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)]
    );
    assert!(bonds.windows(2).all(|pair| pair[0] <= pair[1]));
    assert!(bonds.iter().all(|bond| bond.first() < bond.second()));
}

#[test]
fn tiny_periodic_lattice_preserves_image_multiplicity() {
    let geometry = RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Periodic).unwrap();
    let simple = geometry.bonds(BondMultiplicity::Simple).unwrap();
    let images = geometry.bonds(BondMultiplicity::PeriodicImages).unwrap();
    assert_eq!(simple.len(), 4);
    assert_eq!(images.len(), 8);
    assert!(images.windows(2).all(|pair| pair[0] <= pair[1]));
    assert_eq!(
        images
            .iter()
            .map(|bond| (bond.source().get(), bond.direction()))
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        8
    );
}

#[test]
fn simple_bonds_are_endpoint_only_and_extent_one_skips_self_images() {
    let geometry = RectangularGeometry::new(1, 3, Boundary::Periodic, Boundary::Open).unwrap();
    let bonds = geometry.bonds(BondMultiplicity::Simple).unwrap();
    assert_eq!(
        bonds
            .iter()
            .map(|bond| bond.image_translation())
            .collect::<Vec<_>>(),
        vec![(0, 0), (0, 0)]
    );
    assert_eq!(
        geometry
            .bonds(BondMultiplicity::PeriodicImages)
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn triangular_embedding_has_canonical_coordinates() {
    let geometry = RectangularGeometry::with_kind(
        2,
        2,
        Boundary::Open,
        Boundary::Open,
        LatticeKind::Triangular,
    )
    .unwrap();
    let coordinate = geometry.coordinate(SiteId::new(3)).unwrap();
    assert!((coordinate.x() - 1.5).abs() < 1.0e-15);
    assert!((coordinate.y() - 3.0_f64.sqrt() / 2.0).abs() < 1.0e-15);
}

#[test]
fn triangular_minimum_image_handles_anisotropic_cells() {
    let geometry = RectangularGeometry::with_kind(
        1,
        100,
        Boundary::Periodic,
        Boundary::Periodic,
        LatticeKind::Triangular,
    )
    .unwrap();
    let displacement = geometry
        .minimum_image_displacement(SiteId::new(0), SiteId::new(50))
        .unwrap();
    assert!(displacement.x().abs() < 1.0e-12);
    assert!((displacement.y() + 50.0 * 3.0_f64.sqrt() / 2.0).abs() < 1.0e-12);
}

#[test]
fn triangular_minimum_image_handles_long_x_short_y_cells() {
    let geometry = RectangularGeometry::with_kind(
        100,
        1,
        Boundary::Periodic,
        Boundary::Periodic,
        LatticeKind::Triangular,
    )
    .unwrap();
    let displacement = geometry
        .minimum_image_displacement(SiteId::new(0), SiteId::new(50))
        .unwrap();
    assert!((displacement.x() - 37.5).abs() < 1.0e-12);
    assert!((displacement.y() + 25.0 * 3.0_f64.sqrt() / 2.0).abs() < 1.0e-12);
}

#[test]
fn custom_geometry_preserves_input_order_and_builds_complete_pairs() {
    let geometry = CustomGeometry::new(vec![
        qslib_core::Coordinate::new(2.0, 0.0).unwrap(),
        qslib_core::Coordinate::new(0.0, 1.0).unwrap(),
        qslib_core::Coordinate::new(-1.0, 0.0).unwrap(),
    ])
    .unwrap();
    assert_eq!(
        geometry.coordinate(SiteId::new(1)).unwrap().as_tuple(),
        (0.0, 1.0)
    );
    assert_eq!(geometry.complete_bonds().unwrap().len(), 3);
}

#[test]
fn invalid_geometry_is_rejected_without_wrapping() {
    assert!(RectangularGeometry::new(0, 2, Boundary::Open, Boundary::Open).is_err());
    let geometry = RectangularGeometry::new(3, 2, Boundary::Open, Boundary::Open).unwrap();
    assert!(geometry.site_id(3, 0).is_err());
    assert!(geometry.coordinate(SiteId::new(6)).is_err());
}

#[test]
fn minimum_image_and_shell_selection_are_explicit_and_deterministic() {
    let geometry = RectangularGeometry::new(4, 1, Boundary::Periodic, Boundary::Open).unwrap();
    assert_eq!(
        geometry
            .minimum_image_displacement(SiteId::new(0), SiteId::new(2))
            .unwrap()
            .as_tuple(),
        (-2.0, 0.0)
    );
    assert_eq!(
        geometry
            .pairs_at_squared_distance(1.0, ShellTolerance::Absolute(0.0))
            .unwrap()
            .len(),
        4
    );
    assert_eq!(
        geometry
            .pairs_at_squared_distance(1.0, ShellTolerance::Relative(1.0e-12))
            .unwrap()
            .len(),
        4
    );
}

#[test]
fn legacy_x_major_order_requires_an_explicit_bijection() {
    let adapter = XMajorAdapter::new(3, 2).unwrap();
    assert_eq!(
        adapter.to_canonical(SiteId::new(2)).unwrap(),
        SiteId::new(1)
    );
    assert_eq!(
        adapter.from_canonical(SiteId::new(1)).unwrap(),
        SiteId::new(2)
    );
    assert!(adapter.to_canonical(SiteId::new(6)).is_err());
}
