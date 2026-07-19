use qslib_core::{
    BasisState, Complex64, Hamiltonian, OperatorError, Pauli, PauliString, SiteCount, SiteId,
};

#[test]
fn pauli_string_flips_support_and_tracks_phase() {
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::Y)]).unwrap();
    let (state, coefficient) = operator
        .apply(&BasisState::from_raw_bits(&[0]).unwrap())
        .unwrap();
    assert_eq!(state, BasisState::from_raw_bits(&[1]).unwrap());
    assert_eq!(coefficient, Complex64::new(0.0, 1.0));
}

#[test]
fn canonical_pauli_constructor_rejects_duplicate_support() {
    assert!(
        PauliString::new(vec![(SiteId::new(0), Pauli::X), (SiteId::new(0), Pauli::Z)]).is_err()
    );
}

#[test]
fn ordered_pauli_products_reduce_same_site_factors_with_phase() {
    let (operator, phase) =
        PauliString::product(vec![(SiteId::new(0), Pauli::X), (SiteId::new(0), Pauli::Y)]).unwrap();
    let (state, coefficient) = operator
        .apply(&BasisState::from_raw_bits(&[0]).unwrap())
        .unwrap();
    assert_eq!(state, BasisState::from_raw_bits(&[0]).unwrap());
    assert_eq!(coefficient * phase, Complex64::new(0.0, 1.0));
}

#[test]
fn hermitian_hamiltonian_rejects_nonreal_pauli_coefficients_and_combines_duplicates() {
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap();
    assert!(matches!(
        Hamiltonian::new_hermitian(
            SiteCount::new(1).unwrap(),
            Complex64::new(0.0, 0.0),
            vec![(Complex64::new(0.0, 1.0), operator.clone())]
        ),
        Err(OperatorError::NonHermitianCoefficient { .. })
    ));
    let h = Hamiltonian::new(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![
            (Complex64::new(1.0, 0.0), operator.clone()),
            (Complex64::new(-1.0, 0.0), operator),
        ],
    )
    .unwrap();
    assert!(h.terms().is_empty());
}

#[test]
fn pauli_identity_products_fold_into_the_hamiltonian_constant() {
    let (identity, phase) =
        PauliString::product(vec![(SiteId::new(0), Pauli::X), (SiteId::new(0), Pauli::X)]).unwrap();
    assert!(identity.factors().is_empty());
    assert_eq!(phase, Complex64::new(1.0, 0.0));
    let h = Hamiltonian::new(
        SiteCount::new(1).unwrap(),
        Complex64::new(2.0, 0.0),
        vec![(Complex64::new(3.0, 0.0), identity)],
    )
    .unwrap();
    assert_eq!(h.constant(), Complex64::new(5.0, 0.0));
    assert!(h.terms().is_empty());
}

#[test]
fn local_energy_uses_the_hermitian_matrix_element_for_complex_pauli_action() {
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::Y)]).unwrap();
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(Complex64::new(1.0, 0.0), operator)],
    )
    .unwrap();
    let state_zero = qslib_core::BasisState::from_raw_bits(&[0]).unwrap();
    let state_one = qslib_core::BasisState::from_raw_bits(&[1]).unwrap();
    let amplitudes = vec![
        (state_zero.clone(), Complex64::new(1.0, 0.0)),
        (state_one, Complex64::new(1.0, 0.0)),
    ];
    assert_eq!(
        h.local_energy(&state_zero, &amplitudes).unwrap(),
        Complex64::new(0.0, -1.0)
    );
    let missing = vec![(state_zero.clone(), Complex64::new(1.0, 0.0))];
    assert!(matches!(
        h.local_energy(&state_zero, &missing),
        Err(OperatorError::MissingAmplitude { .. })
    ));
}

#[test]
fn local_energy_preserves_complex_coefficients_and_constants() {
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::X)]).unwrap();
    let h = Hamiltonian::new(
        SiteCount::new(1).unwrap(),
        Complex64::new(2.0, 3.0),
        vec![(Complex64::new(0.0, 1.0), operator)],
    )
    .unwrap();
    let state_zero = qslib_core::BasisState::from_raw_bits(&[0]).unwrap();
    let state_one = qslib_core::BasisState::from_raw_bits(&[1]).unwrap();
    let amplitudes = vec![
        (state_zero.clone(), Complex64::new(1.0, 0.0)),
        (state_one, Complex64::new(1.0, 0.0)),
    ];
    assert_eq!(
        h.local_energy(&state_zero, &amplitudes).unwrap(),
        Complex64::new(2.0, 4.0)
    );
}

#[test]
fn local_energy_cancels_connected_rows_before_requesting_amplitudes() {
    let x = PauliString::new(vec![(SiteId::new(0), Pauli::X)]).unwrap();
    let xz =
        PauliString::new(vec![(SiteId::new(0), Pauli::X), (SiteId::new(1), Pauli::Z)]).unwrap();
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(2).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![
            (Complex64::new(1.0, 0.0), x),
            (Complex64::new(1.0, 0.0), xz),
        ],
    )
    .unwrap();
    let state = qslib_core::BasisState::from_raw_bits(&[0, 1]).unwrap();
    let amplitudes = vec![(state.clone(), Complex64::new(2.0, 0.0))];
    assert_eq!(
        h.local_energy(&state, &amplitudes).unwrap(),
        Complex64::new(0.0, 0.0)
    );
}

#[test]
fn duplicate_term_reduction_is_invariant_under_input_permutation() {
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap();
    let terms = vec![
        (Complex64::new(1.0e16, 0.0), operator.clone()),
        (Complex64::new(-1.0e16, 0.0), operator.clone()),
        (Complex64::new(1.0, 0.0), operator.clone()),
    ];
    let reversed = terms.iter().cloned().rev().collect::<Vec<_>>();
    let first =
        Hamiltonian::new(SiteCount::new(1).unwrap(), Complex64::new(0.0, 0.0), terms).unwrap();
    let second = Hamiltonian::new(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        reversed,
    )
    .unwrap();
    assert_eq!(first, second);
    assert_eq!(first.terms()[0].0, Complex64::new(1.0, 0.0));
}

#[test]
fn tiny_nonzero_coefficients_survive_and_overflow_is_rejected() {
    let operator = PauliString::new(vec![(SiteId::new(0), Pauli::Z)]).unwrap();
    let tiny = Hamiltonian::new(
        SiteCount::new(1).unwrap(),
        Complex64::new(0.0, 0.0),
        vec![(Complex64::new(1.0e-200, 0.0), operator.clone())],
    )
    .unwrap();
    assert_eq!(tiny.terms().len(), 1);
    let state = qslib_core::BasisState::from_raw_bits(&[0]).unwrap();
    assert_eq!(tiny.apply(&state).unwrap().len(), 1);
    assert!(matches!(
        Hamiltonian::new(
            SiteCount::new(1).unwrap(),
            Complex64::new(0.0, 0.0),
            vec![
                (Complex64::new(1.0e308, 0.0), operator.clone()),
                (Complex64::new(1.0e308, 0.0), operator),
            ],
        ),
        Err(OperatorError::NonFiniteCoefficient { .. })
    ));
}
