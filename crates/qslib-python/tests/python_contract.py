"""Specification-first tests for the physicist-facing qslib binding."""

import numpy as np
import pytest
import gc
import threading

import qslib_quantum as qslib


def test_import_contract_and_exception_hierarchy():
    assert qslib.__version__ == "0.1.0"
    assert issubclass(qslib.InputError, qslib.QslibError)
    assert issubclass(qslib.NumericalError, qslib.QslibError)


def test_row_major_geometry_is_owned_and_canonical():
    ids = qslib.row_major_site_ids(2, 3)
    assert ids.dtype == np.uint32
    np.testing.assert_array_equal(ids, [[0, 1], [2, 3], [4, 5]])


def test_pair_resolution_accepts_c_and_fortran_order_and_rejects_bad_values():
    couplings = np.zeros((3, 3), dtype=np.float64, order="F")
    couplings[0, 1] = couplings[1, 0] = 2.0
    pairs, coefficients = qslib.resolve_ising_interactions(couplings)
    np.testing.assert_array_equal(pairs, [[0, 1], [0, 2], [1, 2]])
    np.testing.assert_allclose(coefficients, [2.0, 0.0, 0.0])

    with pytest.raises(qslib.InputError):
        qslib.resolve_ising_interactions(np.zeros((3, 2)))
    with pytest.raises(qslib.InputError):
        qslib.resolve_ising_interactions(np.array([[0.0, np.nan], [np.nan, 0.0]]))


def test_exact_tfim_matrix_and_ground_state_match_one_site_fixture():
    matrix = qslib.tfim_matrix(np.zeros((1, 1)), np.array([2.0]), basis="z")
    np.testing.assert_allclose(matrix, [[0.0, -2.0], [-2.0, 0.0]])
    energy, vector, residual = qslib.tfim_ground_state(
        np.zeros((1, 1)), np.array([2.0]), basis="z"
    )
    assert energy == pytest.approx(-2.0)
    assert vector.dtype == np.complex128
    assert residual < 1.0e-12


def test_four_site_ground_state_fixture_matches_independent_classical_bound():
    couplings = np.zeros((4, 4), dtype=np.float64)
    for first, second, value in ((0, 1, 1.0), (1, 2, 2.0), (2, 3, -0.5)):
        couplings[first, second] = value
        couplings[second, first] = value
    energy, _, residual = qslib.tfim_ground_state(couplings, np.zeros(4))
    assert energy == pytest.approx(-3.5, abs=1.0e-12)
    assert residual < 1.0e-11


def test_observable_contract_reports_totals_densities_and_axis_label():
    result = qslib.tfim_observables(
        np.zeros((1, 1)),
        np.array([2.0]),
        np.array([1.0 + 0.0j, 0.0 + 0.0j]),
        axis="z",
        correlation=(0, 0),
    )
    assert result["energy"] == pytest.approx((0.0, 0.0))
    assert result["energy_density"] == pytest.approx(0.0)
    assert result["energy_variance"] == pytest.approx(4.0)
    assert result["axis"] == "z"
    assert result["magnetization_total"] == pytest.approx(1.0)
    assert result["magnetization_density"] == pytest.approx(1.0)
    assert result["correlation"] == pytest.approx((1.0, 0.0))

    bell = np.array([1.0, 0.0, 0.0, 1.0], dtype=np.complex128) / np.sqrt(2.0)
    distinct = qslib.tfim_observables(
        np.zeros((2, 2)), np.zeros(2), bell, axis="z", correlation=(0, 1)
    )
    assert distinct["correlation"] == pytest.approx((1.0, 0.0))


def test_heisenberg_and_rydberg_matrix_bindings_preserve_hermiticity():
    couplings = np.array([[0.0, 2.0], [2.0, 0.0]])
    heisenberg = qslib.heisenberg_matrix(couplings)
    np.testing.assert_allclose(heisenberg, heisenberg.conj().T)
    rydberg = qslib.rydberg_matrix(couplings, np.array([1.0, 1.0]), np.array([0.0, 0.0]))
    np.testing.assert_allclose(rydberg, rydberg.conj().T)
    np.testing.assert_allclose(
        rydberg,
        [[0.0, -0.5, -0.5, 0.0],
         [-0.5, 0.0, 0.0, -0.5],
         [-0.5, 0.0, 0.0, -0.5],
         [0.0, -0.5, -0.5, 2.0]],
    )

    unit_coupling = np.array([[0.0, 1.0], [1.0, 0.0]])
    np.testing.assert_allclose(
        qslib.tfim_matrix(unit_coupling, np.zeros(2)),
        np.diag([-1.0, 1.0, 1.0, -1.0]),
    )
    np.testing.assert_allclose(
        qslib.heisenberg_matrix(unit_coupling),
        [[0.25, 0.0, 0.0, 0.0],
         [0.0, -0.25, 0.5, 0.0],
         [0.0, 0.5, -0.25, 0.0],
         [0.0, 0.0, 0.0, 0.25]],
    )


def test_tdvp_estimate_returns_mode_specific_owned_statistics():
    weights = np.array([1.0, 1.0])
    local_energies = np.array([1.0 + 0.0j, 3.0 + 0.0j])
    derivatives = np.array([[1.0 + 0.0j], [2.0 + 0.0j]])
    result = qslib.tdvp_estimate(weights, local_energies, derivatives)
    assert result["rhs"].shape == (1,)
    assert result["qgt"].shape == (1, 1)
    assert result["samples"] == 2
    assert result["energy_mean"] == pytest.approx((2.0, 0.0))
    solved = qslib.tdvp_solve(np.array([[2.0]]), np.array([4.0]), shift=0.0)
    assert solved["direction"][0] == pytest.approx(2.0)
    assert solved["converged"]

    complex_energies = np.array([1.0 + 1.0j, 3.0 + 3.0j])
    complex_derivatives = np.array(
        [[1.0 + 0.0j, 0.0 + 0.0j], [3.0 + 0.0j, 2.0 + 0.0j]]
    )
    real_time = qslib.tdvp_estimate(
        weights, complex_energies, complex_derivatives, mode="real_time"
    )
    imaginary_time = qslib.tdvp_estimate(
        weights, complex_energies, complex_derivatives, mode="imaginary_time"
    )
    np.testing.assert_allclose(real_time["qgt"], [[1.0, 1.0], [1.0, 1.0]])
    np.testing.assert_allclose(real_time["rhs"], [1.0, 1.0])
    np.testing.assert_allclose(imaginary_time["rhs"], [-1.0, -1.0])
    assert real_time["energy_variance"] == pytest.approx(2.0)


def test_input_views_are_not_retained_and_nonfinite_fields_fail():
    fields = np.array([2.0, 3.0], dtype=np.float64)
    matrix = np.zeros((2, 2), dtype=np.float64)
    matrix[0, 1] = matrix[1, 0] = 1.0
    qslib.tfim_ground_state(matrix[::-1, ::-1], fields[::-1])
    readonly = np.array(matrix, copy=True)
    readonly.setflags(write=False)
    qslib.tfim_ground_state(readonly, fields)
    with pytest.raises((TypeError, ValueError)):
        qslib.tfim_ground_state(matrix.astype(np.float32), fields)
    with pytest.raises((TypeError, ValueError)):
        qslib.tfim_ground_state(matrix.reshape(4), fields)
    with pytest.raises(qslib.InputError):
        qslib.tfim_ground_state(matrix, np.array([np.inf, 1.0]))


def test_repeated_gc_and_concurrent_calls_do_not_retain_python_buffers():
    matrix = np.zeros((1, 1), dtype=np.float64)
    fields = np.array([2.0], dtype=np.float64)
    failures = []

    def worker():
        try:
            for _ in range(20):
                qslib.tfim_matrix(matrix, fields)
                gc.collect()
        except Exception as error:  # pragma: no cover - failure is asserted below
            failures.append(error)

    threads = [threading.Thread(target=worker) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert failures == []
