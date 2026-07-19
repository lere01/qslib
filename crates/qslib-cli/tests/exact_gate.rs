use qslib::core::{Complex64, Hamiltonian, Pauli, PauliString, SiteCount, SiteId};
use qslib::exact::{ExactBasis, GroundState, diagonalize_hermitian};

#[test]
fn cli_dependency_graph_can_run_the_heterogeneous_ground_state_gate() {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(4).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![
            (-1.0, vec![(0, Pauli::Z), (1, Pauli::Z)]),
            (-2.0, vec![(1, Pauli::Z), (2, Pauli::Z)]),
            (0.5, vec![(2, Pauli::Z), (3, Pauli::Z)]),
        ]
        .into_iter()
        .map(|(coefficient, support)| {
            (
                Complex64::new(coefficient, 0.0),
                PauliString::new(
                    support
                        .into_iter()
                        .map(|(site, pauli)| (SiteId::new(site), pauli))
                        .collect(),
                )
                .unwrap(),
            )
        })
        .collect(),
    )
    .unwrap();
    let basis = ExactBasis::full(SiteCount::new(4).unwrap()).unwrap();
    let matrix = qslib::exact::DenseMatrix::from_hamiltonian(&h, &basis).unwrap();
    let ground = GroundState::from_spectrum(&diagonalize_hermitian(&matrix).unwrap()).unwrap();
    assert!((ground.energy() + 3.5).abs() < 1.0e-12);
}
