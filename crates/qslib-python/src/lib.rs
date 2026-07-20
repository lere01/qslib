//! Physicist-facing Python bindings for validated qslib kernels.
//!
//! The binding owns all Python conversion and returns new NumPy arrays. Rust
//! values are never retained after a call, and the scientific crates remain
//! independent of Python. The initial stable surface covers row-major
//! geometry, dense Ising coupling resolution, and exact TFIM matrix and
//! ground-state inspection.

#![deny(missing_docs)]

use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use qslib::exact::{DenseMatrix, ExactBasis};
use qslib::variational::{DenseQgt, TDVPMode, estimate_tdvp};
use qslib::{
    BasisBit, Boundary, DenseCouplings, InteractionChannel, InteractionTable, PhysicalAxis,
    SiteCount,
};
use std::fmt::Display;

mod exceptions {
    #![allow(missing_docs)]

    use super::PyException;
    use pyo3::create_exception;

    create_exception!(qslib_quantum, QslibError, PyException);
    create_exception!(qslib_quantum, InputError, QslibError);
    create_exception!(qslib_quantum, NumericalError, QslibError);
}

use exceptions::{InputError, NumericalError, QslibError};

type PairResolution = (Py<PyArray2<u32>>, Py<PyArray1<f64>>);

fn input_error(error: impl Display) -> PyErr {
    InputError::new_err(error.to_string())
}

fn numerical_error(error: impl Display) -> PyErr {
    NumericalError::new_err(error.to_string())
}

fn dense_couplings(couplings: &PyReadonlyArray2<'_, f64>) -> PyResult<DenseCouplings> {
    let matrix = couplings.as_array();
    let shape = matrix.shape();
    if shape[0] != shape[1] {
        return Err(input_error("couplings must be square"));
    }
    let n = shape[0];
    let capacity = n
        .checked_mul(n)
        .ok_or_else(|| input_error("coupling matrix dimension overflowed"))?;
    let mut values = Vec::with_capacity(capacity);
    for row in 0..n {
        for column in 0..n {
            values.push(matrix[[row, column]]);
        }
    }
    DenseCouplings::new(SiteCount::new(n).map_err(input_error)?, values).map_err(input_error)
}

fn parse_basis(basis: &str) -> PyResult<qslib::SimulationBasis> {
    basis
        .parse()
        .map_err(|error| input_error(format!("invalid simulation basis: {error}")))
}

fn resolve_tfim(
    couplings: &PyReadonlyArray2<'_, f64>,
    fields: &PyReadonlyArray1<'_, f64>,
    basis: &str,
) -> PyResult<qslib::ResolvedModel> {
    let matrix = couplings.as_array();
    let field_values = fields.as_array();
    let shape = matrix.shape();
    if shape[0] != shape[1] || shape[0] != field_values.len() {
        return Err(input_error(
            "couplings must be square and match fields length",
        ));
    }
    let dense = dense_couplings(couplings)?;
    let site_count = dense.site_count();
    let interactions = dense
        .to_interactions(InteractionChannel::IsingZZ)
        .map_err(input_error)?;
    let table = InteractionTable::new(
        site_count,
        interactions
            .into_iter()
            .map(|term| (term.bond(), term.channel().clone(), term.coefficient()))
            .collect(),
    )
    .map_err(input_error)?;
    let parsed_basis = parse_basis(basis)?;
    qslib::tfim(
        &table,
        field_values.iter().copied().collect::<Vec<_>>().as_slice(),
        parsed_basis,
    )
    .map_err(input_error)
}

fn resolve_heisenberg(
    couplings: &PyReadonlyArray2<'_, f64>,
    basis: &str,
) -> PyResult<qslib::ResolvedModel> {
    let dense = dense_couplings(couplings)?;
    let site_count = dense.site_count();
    let interactions = dense
        .to_interactions(InteractionChannel::HeisenbergExchange)
        .map_err(input_error)?;
    let table = InteractionTable::new(
        site_count,
        interactions
            .into_iter()
            .map(|term| (term.bond(), term.channel().clone(), term.coefficient()))
            .collect(),
    )
    .map_err(input_error)?;
    qslib::heisenberg(&table, parse_basis(basis)?).map_err(input_error)
}

fn resolve_rydberg(
    couplings: &PyReadonlyArray2<'_, f64>,
    omega: &PyReadonlyArray1<'_, f64>,
    detuning: &PyReadonlyArray1<'_, f64>,
) -> PyResult<qslib::ResolvedModel> {
    let dense = dense_couplings(couplings)?;
    if dense.site_count().get() != omega.as_array().len()
        || dense.site_count().get() != detuning.as_array().len()
    {
        return Err(input_error("omega and detuning must match coupling size"));
    }
    let omega_values = omega.as_array().iter().copied().collect::<Vec<_>>();
    let detuning_values = detuning.as_array().iter().copied().collect::<Vec<_>>();
    qslib::rydberg(
        &dense,
        &omega_values,
        &detuning_values,
        qslib::SimulationBasis::Z,
    )
    .map_err(input_error)
}

/// Return canonical row-major site identifiers as a `(Ly, Lx)` owned array.
#[pyfunction]
fn row_major_site_ids(py: Python<'_>, lx: usize, ly: usize) -> PyResult<Py<PyArray2<u32>>> {
    let geometry = qslib::RectangularGeometry::new(lx, ly, Boundary::Open, Boundary::Open)
        .map_err(input_error)?;
    let capacity = lx
        .checked_mul(ly)
        .ok_or_else(|| input_error("geometry dimension overflowed"))?;
    let mut values = Vec::with_capacity(capacity);
    for y in 0..ly {
        for x in 0..lx {
            values.push(geometry.site_id(x, y).map_err(input_error)?.get());
        }
    }
    let array = PyArray2::from_vec2(
        py,
        &(0..ly)
            .map(|row| values[row * lx..(row + 1) * lx].to_vec())
            .collect::<Vec<_>>(),
    )
    .map_err(input_error)?;
    Ok(array.unbind())
}

/// Resolve a symmetric dense Ising coupling matrix to `(pairs, coefficients)`.
///
/// `pairs` is an owned `uint32` array with shape `(terms, 2)` and
/// `coefficients` is an owned `float64` array with one resolved value per row.
#[pyfunction]
fn resolve_ising_interactions(
    py: Python<'_>,
    couplings: PyReadonlyArray2<'_, f64>,
) -> PyResult<PairResolution> {
    let dense = dense_couplings(&couplings)?;
    let interactions = dense
        .to_interactions(InteractionChannel::IsingZZ)
        .map_err(input_error)?;
    let pairs = interactions
        .iter()
        .map(|term| vec![term.bond().first().get(), term.bond().second().get()])
        .collect::<Vec<_>>();
    let coefficients = interactions
        .iter()
        .map(|term| term.coefficient())
        .collect::<Vec<_>>();
    let pair_array = PyArray2::from_vec2(py, &pairs).map_err(input_error)?;
    let coefficient_array = PyArray1::from_vec(py, coefficients).unbind();
    Ok((pair_array.unbind(), coefficient_array))
}

/// Return the canonical full computational basis as a `(dimension, sites)` array.
#[pyfunction]
fn basis_states(
    py: Python<'_>,
    site_count: usize,
    weight: Option<usize>,
) -> PyResult<Py<PyArray2<u8>>> {
    let count = SiteCount::new(site_count).map_err(input_error)?;
    let basis = match weight {
        Some(weight) => ExactBasis::fixed_weight(count, weight),
        None => ExactBasis::full(count),
    }
    .map_err(input_error)?;
    let rows = basis
        .states()
        .iter()
        .map(|state| {
            state
                .bits()
                .iter()
                .map(|bit| if *bit == BasisBit::One { 1 } else { 0 })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    Ok(PyArray2::from_vec2(py, &rows)
        .map_err(input_error)?
        .unbind())
}

fn tfim_matrix_owned(
    couplings: &PyReadonlyArray2<'_, f64>,
    fields: &PyReadonlyArray1<'_, f64>,
    basis: &str,
) -> PyResult<DenseMatrix> {
    let model = resolve_tfim(couplings, fields, basis)?;
    let exact_basis =
        ExactBasis::full(SiteCount::new(fields.as_array().len()).map_err(input_error)?)
            .map_err(input_error)?;
    DenseMatrix::from_hamiltonian(model.hamiltonian(), &exact_basis).map_err(numerical_error)
}

fn matrix_from_model(model: &qslib::ResolvedModel) -> PyResult<DenseMatrix> {
    let exact_basis = ExactBasis::full(model.hamiltonian().site_count()).map_err(input_error)?;
    DenseMatrix::from_hamiltonian(model.hamiltonian(), &exact_basis).map_err(numerical_error)
}

fn matrix_to_rows(matrix: &DenseMatrix) -> Vec<Vec<num_complex::Complex64>> {
    let dimension = matrix.dimension();
    (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix.get(row, column).unwrap_or_default())
                .collect::<Vec<_>>()
        })
        .collect()
}

fn normalized_inner(
    left: &[num_complex::Complex64],
    right: &[num_complex::Complex64],
) -> PyResult<num_complex::Complex64> {
    if left.len() != right.len() {
        return Err(input_error("observable vectors have different dimensions"));
    }
    let norm = left.iter().map(|value| value.norm_sqr()).sum::<f64>();
    if !norm.is_finite() || norm == 0.0 {
        return Err(input_error("state norm must be finite and nonzero"));
    }
    let value = left
        .iter()
        .zip(right)
        .map(|(left, right)| left.conj() * right)
        .sum::<num_complex::Complex64>()
        / norm;
    if !value.re.is_finite() || !value.im.is_finite() {
        return Err(numerical_error("observable expectation is non-finite"));
    }
    Ok(value)
}

/// Return convention-labelled exact TFIM energy moments and one magnetization.
///
/// The state is interpreted in qslib's canonical full-basis order. `axis` is
/// a physical Pauli axis and is independent of the simulation basis. Energies
/// are totals; the returned `energy_density` divides by the number of sites.
/// The magnetization fields contain both total and per-site density values.
#[pyfunction]
#[pyo3(signature=(couplings, fields, state, basis=None, axis="z", correlation=None))]
fn tfim_observables(
    py: Python<'_>,
    couplings: PyReadonlyArray2<'_, f64>,
    fields: PyReadonlyArray1<'_, f64>,
    state: PyReadonlyArray1<'_, num_complex::Complex64>,
    basis: Option<&str>,
    axis: &str,
    correlation: Option<(usize, usize)>,
) -> PyResult<Py<PyDict>> {
    let matrix = tfim_matrix_owned(&couplings, &fields, basis.unwrap_or("z"))?;
    let vector = state.as_array().iter().copied().collect::<Vec<_>>();
    if vector.len() != matrix.dimension() {
        return Err(input_error("state length does not match the exact basis"));
    }
    let applied = matrix.apply(&vector).map_err(numerical_error)?;
    let squared = matrix.apply(&applied).map_err(numerical_error)?;
    let energy = normalized_inner(&vector, &applied)?;
    let second_moment = normalized_inner(&vector, &squared)?;
    let variance = (second_moment - energy * energy).re;
    if !variance.is_finite() || variance < -1.0e-9 {
        return Err(numerical_error("energy variance is invalid"));
    }
    let physical_axis = axis
        .parse::<PhysicalAxis>()
        .map_err(|error| input_error(format!("invalid physical axis: {error}")))?;
    let exact_basis =
        ExactBasis::full(SiteCount::new(fields.as_array().len()).map_err(input_error)?)
            .map_err(input_error)?;
    let magnetization =
        qslib::exact::magnetization_expectation(physical_axis, &exact_basis, &vector)
            .map_err(numerical_error)?;
    let result = PyDict::new(py);
    result.set_item("energy", (energy.re, energy.im))?;
    result.set_item("energy_density", energy.re / fields.as_array().len() as f64)?;
    result.set_item("energy_variance", variance.max(0.0))?;
    result.set_item("axis", physical_axis.as_str())?;
    result.set_item("magnetization_total", magnetization.total())?;
    result.set_item("magnetization_density", magnetization.density())?;
    if let Some((first, second)) = correlation {
        let value = qslib::exact::correlation_expectation(
            physical_axis,
            first,
            second,
            &exact_basis,
            &vector,
        )
        .map_err(input_error)?;
        result.set_item("correlation", (value.re, value.im))?;
    }
    Ok(result.unbind())
}

/// Build an exact complex TFIM matrix from owned or strided NumPy inputs.
#[pyfunction]
#[pyo3(signature=(couplings, fields, basis=None))]
fn tfim_matrix(
    py: Python<'_>,
    couplings: PyReadonlyArray2<'_, f64>,
    fields: PyReadonlyArray1<'_, f64>,
    basis: Option<&str>,
) -> PyResult<Py<PyArray2<num_complex::Complex64>>> {
    let matrix = tfim_matrix_owned(&couplings, &fields, basis.unwrap_or("z"))?;
    let dimension = matrix.dimension();
    let rows = (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| matrix.get(row, column).unwrap_or_default())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    Ok(PyArray2::from_vec2(py, &rows)
        .map_err(numerical_error)?
        .unbind())
}

/// Build an exact complex isotropic Heisenberg matrix from pair couplings.
#[pyfunction]
#[pyo3(signature=(couplings, basis=None))]
fn heisenberg_matrix(
    py: Python<'_>,
    couplings: PyReadonlyArray2<'_, f64>,
    basis: Option<&str>,
) -> PyResult<Py<PyArray2<num_complex::Complex64>>> {
    let matrix = matrix_from_model(&resolve_heisenberg(&couplings, basis.unwrap_or("z"))?)?;
    Ok(PyArray2::from_vec2(py, &matrix_to_rows(&matrix))
        .map_err(numerical_error)?
        .unbind())
}

/// Build an exact complex Rydberg matrix in the canonical z basis.
#[pyfunction]
fn rydberg_matrix(
    py: Python<'_>,
    couplings: PyReadonlyArray2<'_, f64>,
    omega: PyReadonlyArray1<'_, f64>,
    detuning: PyReadonlyArray1<'_, f64>,
) -> PyResult<Py<PyArray2<num_complex::Complex64>>> {
    let matrix = matrix_from_model(&resolve_rydberg(&couplings, &omega, &detuning)?)?;
    Ok(PyArray2::from_vec2(py, &matrix_to_rows(&matrix))
        .map_err(numerical_error)?
        .unbind())
}

/// Estimate real-parameter TDVP statistics from owned sample arrays.
#[pyfunction]
#[pyo3(signature=(weights, local_energies, derivatives, mode="real_time"))]
fn tdvp_estimate(
    py: Python<'_>,
    weights: PyReadonlyArray1<'_, f64>,
    local_energies: PyReadonlyArray1<'_, num_complex::Complex64>,
    derivatives: PyReadonlyArray2<'_, num_complex::Complex64>,
    mode: &str,
) -> PyResult<Py<PyDict>> {
    let weights = weights.as_array();
    let energies = local_energies.as_array();
    let derivative_view = derivatives.as_array();
    if derivative_view.shape()[0] != weights.len() || derivative_view.shape()[0] != energies.len() {
        return Err(input_error("TDVP sample dimensions do not match"));
    }
    let parameter_count = derivative_view.shape()[1];
    let derivative_values = (0..derivative_view.shape()[0])
        .flat_map(|sample| {
            (0..parameter_count).map(move |parameter| derivative_view[[sample, parameter]])
        })
        .collect::<Vec<_>>();
    let tdvp_mode = match mode {
        "real_time" => TDVPMode::RealTime,
        "imaginary_time" => TDVPMode::ImaginaryTime,
        value => return Err(input_error(format!("invalid TDVP mode: {value}"))),
    };
    let statistics = estimate_tdvp(
        &weights.iter().copied().collect::<Vec<_>>(),
        &energies.iter().copied().collect::<Vec<_>>(),
        &derivative_values,
        parameter_count,
        tdvp_mode,
    )
    .map_err(numerical_error)?;
    let qgt_rows = (0..parameter_count)
        .map(|row| {
            statistics.qgt().as_slice()[row * parameter_count..(row + 1) * parameter_count].to_vec()
        })
        .collect::<Vec<_>>();
    let result = PyDict::new(py);
    result.set_item("rhs", PyArray1::from_vec(py, statistics.rhs().to_vec()))?;
    result.set_item(
        "qgt",
        PyArray2::from_vec2(py, &qgt_rows).map_err(numerical_error)?,
    )?;
    result.set_item(
        "energy_mean",
        (statistics.energy_mean().re, statistics.energy_mean().im),
    )?;
    result.set_item("energy_variance", statistics.energy_variance())?;
    result.set_item("samples", statistics.samples())?;
    Ok(result.unbind())
}

/// Solve a checked dense QGT system `(S + shift I) direction = rhs`.
#[pyfunction]
#[pyo3(signature=(qgt, rhs, shift=1.0e-8, tolerance=1.0e-10, max_iterations=256))]
fn tdvp_solve(
    py: Python<'_>,
    qgt: PyReadonlyArray2<'_, f64>,
    rhs: PyReadonlyArray1<'_, f64>,
    shift: f64,
    tolerance: f64,
    max_iterations: usize,
) -> PyResult<Py<PyDict>> {
    let qgt_view = qgt.as_array();
    let shape = qgt_view.shape();
    if shape[0] != shape[1] || shape[0] != rhs.as_array().len() {
        return Err(input_error("QGT and right-hand side shapes do not match"));
    }
    let values = (0..shape[0])
        .flat_map(|row| (0..shape[1]).map(move |column| qgt_view[[row, column]]))
        .collect::<Vec<_>>();
    let matrix = DenseQgt::new(shape[0], values).map_err(numerical_error)?;
    let result = matrix
        .solve(
            &rhs.as_array().iter().copied().collect::<Vec<_>>(),
            shift,
            tolerance,
            max_iterations,
        )
        .map_err(numerical_error)?;
    let output = PyDict::new(py);
    output.set_item(
        "direction",
        PyArray1::from_vec(py, result.direction().to_vec()),
    )?;
    output.set_item("residual_norm", result.residual_norm())?;
    output.set_item("iterations", result.iterations())?;
    output.set_item("converged", result.converged())?;
    Ok(output.unbind())
}

/// Return `(energy, complex_vector, residual)` for the exact TFIM ground state.
#[pyfunction]
#[pyo3(signature=(couplings, fields, basis=None))]
fn tfim_ground_state(
    py: Python<'_>,
    couplings: PyReadonlyArray2<'_, f64>,
    fields: PyReadonlyArray1<'_, f64>,
    basis: Option<&str>,
) -> PyResult<(f64, Py<PyArray1<num_complex::Complex64>>, f64)> {
    let matrix = tfim_matrix_owned(&couplings, &fields, basis.unwrap_or("z"))?;
    let spectrum = qslib::exact::diagonalize_hermitian(&matrix).map_err(numerical_error)?;
    let ground = qslib::exact::GroundState::from_spectrum(&spectrum).map_err(numerical_error)?;
    let vector = PyArray1::from_vec(py, ground.vector().to_vec()).unbind();
    Ok((ground.energy(), vector, ground.residual()))
}

/// Initialize the `qslib_quantum` Python module.
#[pymodule]
fn qslib_quantum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("QslibError", m.py().get_type::<QslibError>())?;
    m.add("InputError", m.py().get_type::<InputError>())?;
    m.add("NumericalError", m.py().get_type::<NumericalError>())?;
    m.add_function(wrap_pyfunction!(row_major_site_ids, m)?)?;
    m.add_function(wrap_pyfunction!(resolve_ising_interactions, m)?)?;
    m.add_function(wrap_pyfunction!(basis_states, m)?)?;
    m.add_function(wrap_pyfunction!(tfim_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tfim_ground_state, m)?)?;
    m.add_function(wrap_pyfunction!(tfim_observables, m)?)?;
    m.add_function(wrap_pyfunction!(heisenberg_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(rydberg_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(tdvp_estimate, m)?)?;
    m.add_function(wrap_pyfunction!(tdvp_solve, m)?)?;
    Ok(())
}
