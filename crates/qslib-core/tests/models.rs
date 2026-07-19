use qslib_core::{
    BasisState, Boundary, DenseCouplings, InteractionChannel, InteractionTable, ModelError,
    RectangularGeometry, SimulationBasis, SiteCount, SiteId,
};
use qslib_test_support::{ComparisonPolicy, FixtureKind, load_conformance_fixtures};
use serde_json::Value;

fn bond_table(channel: InteractionChannel, coefficient: f64) -> InteractionTable {
    let geometry = RectangularGeometry::new(2, 1, Boundary::Open, Boundary::Open).unwrap();
    let bond = geometry
        .bonds(qslib_core::BondMultiplicity::Simple)
        .unwrap()[0];
    InteractionTable::new(
        SiteCount::new(2).unwrap(),
        vec![(bond, channel, coefficient)],
    )
    .unwrap()
}

fn three_site_table() -> InteractionTable {
    let b01 = qslib_core::Bond::new(SiteId::new(0), SiteId::new(1)).unwrap();
    let b12 = qslib_core::Bond::new(SiteId::new(1), SiteId::new(2)).unwrap();
    InteractionTable::new(
        SiteCount::new(3).unwrap(),
        vec![
            (b01, InteractionChannel::HeisenbergExchange, 2.0),
            (b12, InteractionChannel::HeisenbergExchange, -3.0),
        ],
    )
    .unwrap()
}

fn three_site_ising_table() -> InteractionTable {
    let b01 = qslib_core::Bond::new(SiteId::new(0), SiteId::new(1)).unwrap();
    let b12 = qslib_core::Bond::new(SiteId::new(1), SiteId::new(2)).unwrap();
    InteractionTable::new(
        SiteCount::new(3).unwrap(),
        vec![
            (b01, InteractionChannel::IsingZZ, 2.0),
            (b12, InteractionChannel::IsingZZ, -3.0),
        ],
    )
    .unwrap()
}

fn matrix(model: &qslib_core::Hamiltonian, n: usize) -> Vec<Vec<qslib_core::Complex64>> {
    let dimension = 1usize << n;
    let mut columns = vec![vec![qslib_core::Complex64::new(0.0, 0.0); dimension]; dimension];
    for (mask, column) in columns.iter_mut().enumerate() {
        let bits = (0..n)
            .map(|site| ((mask >> site) & 1) as u8)
            .collect::<Vec<_>>();
        let state = BasisState::from_raw_bits(&bits).unwrap();
        for (connected, coefficient) in model.apply(&state).unwrap() {
            let connected_mask = connected
                .bits()
                .iter()
                .enumerate()
                .map(|(site, bit)| (bit.as_u8() as usize) << site)
                .sum::<usize>();
            column[connected_mask] += coefficient;
        }
    }
    let mut result = vec![vec![qslib_core::Complex64::new(0.0, 0.0); dimension]; dimension];
    for (row, result_row) in result.iter_mut().enumerate() {
        for (value, column) in result_row.iter_mut().zip(columns.iter()) {
            *value = column[row];
        }
    }
    result
}

fn fixture(kind: FixtureKind) -> qslib_test_support::Fixture {
    load_conformance_fixtures()
        .unwrap()
        .into_iter()
        .find(|candidate| candidate.kind == kind)
        .unwrap()
}

fn fixture_matrix(data: &Value, key: &str) -> Vec<Vec<qslib_core::Complex64>> {
    let matrix = &data[key];
    let shape = matrix["shape"].as_array().unwrap();
    let rows = shape[0].as_u64().unwrap() as usize;
    let columns = shape[1].as_u64().unwrap() as usize;
    let entries = matrix["entries"].as_array().unwrap();
    assert_eq!(entries.len(), rows * columns);
    entries
        .chunks(columns)
        .map(|row| {
            row.iter()
                .map(|entry| {
                    qslib_core::Complex64::new(
                        entry["re"].as_f64().unwrap(),
                        entry["im"].as_f64().unwrap(),
                    )
                })
                .collect()
        })
        .collect()
}

fn assert_matrix_matches_fixture(
    actual: &[Vec<qslib_core::Complex64>],
    expected: &[Vec<qslib_core::Complex64>],
    policy: &ComparisonPolicy,
) {
    assert_eq!(actual.len(), expected.len());
    for (actual_row, expected_row) in actual.iter().zip(expected) {
        assert_eq!(actual_row.len(), expected_row.len());
        for (actual_value, expected_value) in actual_row.iter().zip(expected_row) {
            if policy.is_exact() {
                assert_eq!(*actual_value, *expected_value);
            } else {
                let difference = (*actual_value - *expected_value).norm();
                let scale = actual_value.norm().max(expected_value.norm());
                assert!(
                    difference <= policy.absolute_tolerance() + policy.relative_tolerance() * scale
                );
            }
        }
    }
}

#[test]
fn model_matrices_match_independent_neutral_fixtures() {
    let tfim_fixture = fixture(FixtureKind::TfimOneBond);
    let tfim_model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 0.5],
        SimulationBasis::Z,
    )
    .unwrap();
    assert_matrix_matches_fixture(
        &matrix(&tfim_model, 2),
        &fixture_matrix(&tfim_fixture.data, "matrix"),
        tfim_fixture.comparison(),
    );

    let heisenberg_fixture = fixture(FixtureKind::HeisenbergOneBond);
    let heisenberg_model = qslib_core::heisenberg(
        &bond_table(InteractionChannel::HeisenbergExchange, 2.0),
        SimulationBasis::Z,
    )
    .unwrap();
    assert_matrix_matches_fixture(
        &matrix(&heisenberg_model, 2),
        &fixture_matrix(&heisenberg_fixture.data, "matrix"),
        heisenberg_fixture.comparison(),
    );

    let rydberg_fixture = fixture(FixtureKind::RydbergTwoSite);
    let couplings =
        DenseCouplings::new(SiteCount::new(2).unwrap(), vec![0.0, 5.0, 5.0, 0.0]).unwrap();
    let rydberg_model =
        qslib_core::rydberg(&couplings, &[2.0, 2.0], &[3.0, 3.0], SimulationBasis::Z).unwrap();
    assert_matrix_matches_fixture(
        &matrix(&rydberg_model, 2),
        &fixture_matrix(&rydberg_fixture.data, "matrix"),
        rydberg_fixture.comparison(),
    );
}

#[test]
fn basis_rotation_matrices_match_the_independent_hadamard_fixture() {
    let basis_fixture = fixture(FixtureKind::BasisRotationSpectrum);
    let z_model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 1.0),
        &[0.5, 0.5],
        SimulationBasis::Z,
    )
    .unwrap();
    let x_model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 1.0),
        &[0.5, 0.5],
        SimulationBasis::X,
    )
    .unwrap();
    assert_matrix_matches_fixture(
        &matrix(&z_model, 2),
        &fixture_matrix(&basis_fixture.data, "z_basis_matrix"),
        basis_fixture.comparison(),
    );
    assert_matrix_matches_fixture(
        &matrix(&x_model, 2),
        &fixture_matrix(&basis_fixture.data, "x_basis_matrix"),
        basis_fixture.comparison(),
    );
}

#[test]
fn heterogeneous_heisenberg_action_matches_the_independent_fixture() {
    let fixture = fixture(FixtureKind::HeisenbergHeterogeneous);
    let input = BasisState::from_raw_bits(&[0, 1, 0]).unwrap();
    let model = qslib_core::heisenberg(&three_site_table(), SimulationBasis::Z).unwrap();
    let action = model.apply(&input).unwrap();
    assert_eq!(
        action.iter().find(|(state, _)| *state == input).unwrap().1,
        qslib_core::Complex64::new(fixture.data["diagonal_energy"].as_f64().unwrap(), 0.0)
    );
    for connection in fixture.data["connections"].as_array().unwrap() {
        let bits = connection["output_bits"]
            .as_array()
            .unwrap()
            .iter()
            .map(|bit| bit.as_u64().unwrap() as u8)
            .collect::<Vec<_>>();
        let output = BasisState::from_raw_bits(&bits).unwrap();
        assert_eq!(
            action.iter().find(|(state, _)| *state == output).unwrap().1,
            qslib_core::Complex64::new(connection["matrix_element"].as_f64().unwrap(), 0.0)
        );
    }
}

#[test]
fn tfim_local_energy_matches_the_fixture_row_definition() {
    let fixture = fixture(FixtureKind::TfimOneBond);
    let model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 0.5],
        SimulationBasis::Z,
    )
    .unwrap();
    let states = (0..4)
        .map(|mask| {
            BasisState::from_raw_bits(&[(mask & 1) as u8, ((mask >> 1) & 1) as u8]).unwrap()
        })
        .collect::<Vec<_>>();
    let amplitudes = states
        .iter()
        .enumerate()
        .map(|(index, state)| {
            (
                state.clone(),
                qslib_core::Complex64::new((index + 1) as f64, 0.0),
            )
        })
        .collect::<Vec<_>>();
    let expected_row = fixture_matrix(&fixture.data, "matrix")[0]
        .iter()
        .zip(amplitudes.iter())
        .map(|(matrix_element, (_, amplitude))| *matrix_element * *amplitude)
        .fold(qslib_core::Complex64::new(0.0, 0.0), |sum, value| {
            sum + value
        });
    let actual = model.local_energy(&states[0], &amplitudes).unwrap();
    assert_eq!(actual, expected_row);
}

#[test]
fn tfim_z_basis_matches_one_bond_matrix_elements() {
    let model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 1.0],
        SimulationBasis::Z,
    )
    .unwrap();
    let action = model
        .apply(&BasisState::from_raw_bits(&[0, 0]).unwrap())
        .unwrap();
    assert_eq!(action.len(), 3);
    assert_eq!(action[0].0, BasisState::from_raw_bits(&[0, 0]).unwrap());
    assert_eq!(action[0].1, qslib_core::Complex64::new(-2.0, 0.0));
    assert_eq!(action[1].1, qslib_core::Complex64::new(-0.5, 0.0));
    assert_eq!(action[2].1, qslib_core::Complex64::new(-1.0, 0.0));
}

#[test]
fn heisenberg_z_basis_matches_exchange_and_flip_matrix_elements() {
    let model = qslib_core::heisenberg(
        &bond_table(InteractionChannel::HeisenbergExchange, 4.0),
        SimulationBasis::Z,
    )
    .unwrap();
    let action = model
        .apply(&BasisState::from_raw_bits(&[0, 1]).unwrap())
        .unwrap();
    assert_eq!(action.len(), 2);
    assert_eq!(
        action
            .iter()
            .find(|(state, _)| *state == BasisState::from_raw_bits(&[0, 1]).unwrap())
            .unwrap()
            .1,
        qslib_core::Complex64::new(-1.0, 0.0)
    );
    assert_eq!(
        action
            .iter()
            .find(|(state, _)| *state == BasisState::from_raw_bits(&[1, 0]).unwrap())
            .unwrap()
            .1,
        qslib_core::Complex64::new(2.0, 0.0)
    );
}

#[test]
fn rydberg_occupation_and_pair_energy_use_canonical_bit_meaning() {
    let couplings =
        DenseCouplings::new(SiteCount::new(2).unwrap(), vec![0.0, 5.0, 5.0, 0.0]).unwrap();
    let model =
        qslib_core::rydberg(&couplings, &[2.0, 2.0], &[3.0, 3.0], SimulationBasis::Z).unwrap();
    let action = model
        .apply(&BasisState::from_raw_bits(&[1, 1]).unwrap())
        .unwrap();
    assert_eq!(
        action
            .iter()
            .find(|(state, _)| *state == BasisState::from_raw_bits(&[1, 1]).unwrap())
            .unwrap()
            .1,
        qslib_core::Complex64::new(-1.0, 0.0)
    );
    assert!(action.iter().any(|(state, coefficient)| *state
        == BasisState::from_raw_bits(&[0, 1]).unwrap()
        && *coefficient == qslib_core::Complex64::new(-1.0, 0.0)));
}

#[test]
fn models_reject_wrong_field_lengths_and_unsupported_basis() {
    assert!(matches!(
        qslib_core::tfim(
            &bond_table(InteractionChannel::IsingZZ, 1.0),
            &[1.0],
            SimulationBasis::Z
        ),
        Err(ModelError::FieldLength { .. })
    ));
    assert!(
        qslib_core::heisenberg(
            &bond_table(InteractionChannel::HeisenbergExchange, 1.0),
            SimulationBasis::X
        )
        .is_ok()
    );
    assert!(matches!(
        qslib_core::tfim(
            &bond_table(InteractionChannel::HeisenbergExchange, 1.0),
            &[1.0, 1.0],
            SimulationBasis::Z
        ),
        Err(ModelError::UnexpectedChannel { .. })
    ));
    let _ = SiteId::new(0);
}

#[test]
fn tfim_x_basis_rotates_physical_axes_without_changing_spectrum() {
    let z = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 1.0],
        SimulationBasis::Z,
    )
    .unwrap();
    let x = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 1.0],
        SimulationBasis::X,
    )
    .unwrap();
    let state = BasisState::from_raw_bits(&[0, 0]).unwrap();
    assert_eq!(
        z.apply(&state).unwrap()[0].1,
        qslib_core::Complex64::new(-2.0, 0.0)
    );
    assert_eq!(
        x.apply(&state).unwrap()[0].1,
        qslib_core::Complex64::new(-1.5, 0.0)
    );
    let z_matrix = matrix(&z, 2);
    let x_matrix = matrix(&x, 2);
    let trace = |values: &[Vec<qslib_core::Complex64>]| {
        (0..values.len())
            .map(|index| values[index][index])
            .fold(qslib_core::Complex64::new(0.0, 0.0), |sum, value| {
                sum + value
            })
    };
    let trace_square = |values: &[Vec<qslib_core::Complex64>]| {
        (0..values.len())
            .flat_map(|row| {
                (0..values.len()).map(move |column| values[row][column] * values[column][row])
            })
            .fold(qslib_core::Complex64::new(0.0, 0.0), |sum, value| {
                sum + value
            })
    };
    assert_eq!(trace(&z_matrix), trace(&x_matrix));
    assert_eq!(trace_square(&z_matrix), trace_square(&x_matrix));
}

#[test]
fn heterogeneous_tfim_preserves_each_pair_coefficient_and_provenance() {
    let model = qslib_core::tfim(
        &three_site_ising_table(),
        &[0.0, 0.0, 0.0],
        SimulationBasis::Z,
    )
    .unwrap();
    let state = BasisState::from_raw_bits(&[0, 1, 0]).unwrap();
    let action = model.apply(&state).unwrap();
    assert_eq!(
        action
            .iter()
            .find(|(candidate, _)| *candidate == state)
            .unwrap()
            .1,
        qslib_core::Complex64::new(-1.0, 0.0)
    );
    assert_eq!(model.family(), "tfim");
    assert_eq!(model.basis(), SimulationBasis::Z);
    assert_eq!(model.interactions().len(), 2);
    assert_eq!(model.interactions()[0].coefficient(), 2.0);
    assert_eq!(model.interactions()[1].coefficient(), -3.0);
    assert!(matches!(
        model.specification(),
        qslib_core::ModelSpecification::Tfim { fields } if fields == &[0.0, 0.0, 0.0]
    ));
}

#[test]
fn explicit_tfim_matrix_is_hermitian_and_preserves_heterogeneous_fields() {
    let model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 1.0],
        SimulationBasis::Z,
    )
    .unwrap();
    let matrix = matrix(&model, 2);
    assert_eq!(matrix[0][0], qslib_core::Complex64::new(-2.0, 0.0));
    assert_eq!(matrix[1][0], qslib_core::Complex64::new(-0.5, 0.0));
    assert_eq!(matrix[2][0], qslib_core::Complex64::new(-1.0, 0.0));
    assert!(matrix.iter().enumerate().all(|(row, values)| {
        values
            .iter()
            .enumerate()
            .all(|(column, value)| *value == matrix[column][row].conj())
    }));
}

#[test]
fn explicit_heisenberg_matrix_has_exchange_block_and_is_hermitian() {
    let model = qslib_core::heisenberg(
        &bond_table(InteractionChannel::HeisenbergExchange, 4.0),
        SimulationBasis::Z,
    )
    .unwrap();
    let z_matrix = matrix(&model, 2);
    assert_eq!(z_matrix[0][0], qslib_core::Complex64::new(1.0, 0.0));
    assert_eq!(z_matrix[1][1], qslib_core::Complex64::new(-1.0, 0.0));
    assert_eq!(z_matrix[1][2], qslib_core::Complex64::new(2.0, 0.0));
    assert_eq!(z_matrix[2][1], qslib_core::Complex64::new(2.0, 0.0));
    assert!(z_matrix.iter().enumerate().all(|(row, values)| {
        values
            .iter()
            .enumerate()
            .all(|(column, value)| *value == z_matrix[column][row].conj())
    }));
    let x_matrix = matrix(
        &qslib_core::heisenberg(
            &bond_table(InteractionChannel::HeisenbergExchange, 4.0),
            SimulationBasis::X,
        )
        .unwrap(),
        2,
    );
    assert_eq!(z_matrix, x_matrix);
}

#[test]
fn explicit_rydberg_matrix_has_occupation_diagonal_and_drive_connections() {
    let couplings =
        DenseCouplings::new(SiteCount::new(2).unwrap(), vec![0.0, 5.0, 5.0, 0.0]).unwrap();
    let model =
        qslib_core::rydberg(&couplings, &[2.0, 2.0], &[3.0, 3.0], SimulationBasis::Z).unwrap();
    let matrix = matrix(&model, 2);
    assert_eq!(matrix[0][0], qslib_core::Complex64::new(0.0, 0.0));
    assert_eq!(matrix[1][1], qslib_core::Complex64::new(-3.0, 0.0));
    assert_eq!(matrix[2][2], qslib_core::Complex64::new(-3.0, 0.0));
    assert_eq!(matrix[3][3], qslib_core::Complex64::new(-1.0, 0.0));
    assert_eq!(matrix[1][0], qslib_core::Complex64::new(-1.0, 0.0));
    assert_eq!(matrix[2][0], qslib_core::Complex64::new(-1.0, 0.0));
}

#[test]
fn heterogeneous_rydberg_pairs_enter_the_occupation_energy_independently() {
    let couplings = DenseCouplings::new(
        SiteCount::new(3).unwrap(),
        vec![0.0, 2.0, -3.0, 2.0, 0.0, 0.0, -3.0, 0.0, 0.0],
    )
    .unwrap();
    let model = qslib_core::rydberg(
        &couplings,
        &[0.0, 0.0, 0.0],
        &[0.0, 0.0, 0.0],
        SimulationBasis::Z,
    )
    .unwrap();
    let matrix = matrix(&model, 3);
    assert_eq!(matrix[7][7], qslib_core::Complex64::new(-1.0, 0.0));
    assert_eq!(model.interactions().len(), 3);
    assert!(matches!(
        model.specification(),
        qslib_core::ModelSpecification::Rydberg { omega, detuning }
            if omega == &[0.0, 0.0, 0.0] && detuning == &[0.0, 0.0, 0.0]
    ));
}

#[test]
fn heterogeneous_three_site_heisenberg_keeps_signed_pair_coefficients() {
    let model = qslib_core::heisenberg(&three_site_table(), SimulationBasis::Z).unwrap();
    let action = model
        .apply(&BasisState::from_raw_bits(&[0, 1, 0]).unwrap())
        .unwrap();
    assert_eq!(
        action
            .iter()
            .find(|(state, _)| *state == BasisState::from_raw_bits(&[0, 1, 0]).unwrap())
            .unwrap()
            .1,
        qslib_core::Complex64::new(0.25, 0.0)
    );
    assert!(action.iter().any(|(state, coefficient)| *state
        == BasisState::from_raw_bits(&[1, 0, 0]).unwrap()
        && *coefficient == qslib_core::Complex64::new(1.0, 0.0)));
    assert!(action.iter().any(|(state, coefficient)| *state
        == BasisState::from_raw_bits(&[0, 0, 1]).unwrap()
        && *coefficient == qslib_core::Complex64::new(-1.5, 0.0)));
}

#[test]
fn j1j2_shorthand_resolves_axial_and_diagonal_shells() {
    let geometry = RectangularGeometry::new(3, 3, Boundary::Open, Boundary::Open).unwrap();
    let model = qslib_core::j1j2(&geometry, 2.0, -3.0, SimulationBasis::Z).unwrap();
    let j1_bonds = geometry
        .bonds(qslib_core::BondMultiplicity::Simple)
        .unwrap();
    let j2_bonds = geometry
        .pairs_at_squared_distance(2.0, qslib_core::ShellTolerance::Absolute(0.0))
        .unwrap();
    assert_eq!(model.interactions().len(), j1_bonds.len() + j2_bonds.len());
    for bond in &j1_bonds {
        let term = model
            .interactions()
            .iter()
            .find(|term| term.identity().name() == Some("j1") && term.bond() == *bond)
            .unwrap();
        assert_eq!(term.coefficient(), 2.0);
    }
    for bond in &j2_bonds {
        let term = model
            .interactions()
            .iter()
            .find(|term| term.identity().name() == Some("j2") && term.bond() == *bond)
            .unwrap();
        assert_eq!(term.coefficient(), -3.0);
    }
    assert_eq!(model.family(), "j1j2");
    assert!(matches!(
        model.specification(),
        qslib_core::ModelSpecification::J1J2 { geometry: saved, .. } if saved == &geometry
    ));
}

#[test]
fn j1j2_shorthand_rejects_triangular_geometry() {
    let geometry = RectangularGeometry::with_kind(
        3,
        3,
        Boundary::Open,
        Boundary::Open,
        qslib_core::LatticeKind::Triangular,
    )
    .unwrap();
    assert!(qslib_core::j1j2(&geometry, 1.0, 0.5, SimulationBasis::Z).is_err());
}

#[test]
fn disordered_j1j2_preserves_shell_identity_signed_zero_and_periodic_additivity() {
    let geometry = RectangularGeometry::new(2, 2, Boundary::Periodic, Boundary::Periodic).unwrap();
    let j1_bonds = geometry
        .bonds(qslib_core::BondMultiplicity::Simple)
        .unwrap();
    let j2_bonds = geometry
        .pairs_at_squared_distance(2.0, qslib_core::ShellTolerance::Absolute(0.0))
        .unwrap();
    let j1 = vec![2.0, -1.0, 0.0, 3.0][..j1_bonds.len()].to_vec();
    let j2 = vec![4.0, -5.0][..j2_bonds.len()].to_vec();
    let model = qslib_core::j1j2_disordered(&geometry, &j1, &j2, SimulationBasis::Z).unwrap();
    assert_eq!(model.interactions().len(), j1.len() + j2.len());
    assert_eq!(
        model
            .interactions()
            .iter()
            .filter(|term| term.identity().name() == Some("j1"))
            .count(),
        j1.len()
    );
    assert_eq!(
        model
            .interactions()
            .iter()
            .filter(|term| term.identity().name() == Some("j2"))
            .count(),
        j2.len()
    );
    assert!(
        model
            .interactions()
            .iter()
            .any(|term| term.coefficient() == 0.0)
    );
    assert!(matches!(
        qslib_core::j1j2_disordered(&geometry, &j1[..j1.len() - 1], &j2, SimulationBasis::Z),
        Err(ModelError::ShellLength { shell: "j1", .. })
    ));
    let bond = j1_bonds[0];
    let additive_table = InteractionTable::new_with_identities(
        geometry.site_count(),
        vec![
            (
                qslib_core::InteractionIdentity::named(
                    bond,
                    InteractionChannel::HeisenbergExchange,
                    "j1",
                )
                .unwrap(),
                2.0,
            ),
            (
                qslib_core::InteractionIdentity::named(
                    bond,
                    InteractionChannel::HeisenbergExchange,
                    "j2",
                )
                .unwrap(),
                3.0,
            ),
        ],
    )
    .unwrap();
    let additive = qslib_core::heisenberg(&additive_table, SimulationBasis::Z).unwrap();
    assert_eq!(additive.terms().len(), 3);
    assert!(
        additive
            .terms()
            .iter()
            .all(|(coefficient, _)| *coefficient == qslib_core::Complex64::new(1.25, 0.0))
    );
}

#[test]
fn tfim_local_energy_uses_wavefunction_ratios() {
    let model = qslib_core::tfim(
        &bond_table(InteractionChannel::IsingZZ, 2.0),
        &[0.5, 1.0],
        SimulationBasis::Z,
    )
    .unwrap();
    let amplitudes = vec![
        (
            BasisState::from_raw_bits(&[0, 0]).unwrap(),
            qslib_core::Complex64::new(1.0, 0.0),
        ),
        (
            BasisState::from_raw_bits(&[1, 0]).unwrap(),
            qslib_core::Complex64::new(2.0, 0.0),
        ),
        (
            BasisState::from_raw_bits(&[0, 1]).unwrap(),
            qslib_core::Complex64::new(3.0, 0.0),
        ),
    ];
    assert_eq!(
        model.local_energy(&amplitudes[0].0, &amplitudes).unwrap(),
        qslib_core::Complex64::new(-6.0, 0.0)
    );
}
