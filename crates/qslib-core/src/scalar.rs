use crate::BasisError;

/// The reference real scalar for qslib scientific calculations.
pub type Real = f64;

/// The reference complex scalar for qslib scientific calculations.
pub type Complex64 = num_complex::Complex<f64>;

/// Reject a non-finite real value before it enters a scientific calculation.
pub fn ensure_finite(value: Real) -> Result<Real, BasisError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BasisError::NonFiniteScalar { value })
    }
}
