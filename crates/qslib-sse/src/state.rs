//! Canonical bit states and trace-closed padded operator strings.

use crate::{BasisBit, OperatorKind, SseModel, SseModelError};

/// One padded position in an SSE operator string.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Operator {
    /// Unused identity position.
    Identity,
    /// Diagonal term index.
    Diagonal(u32),
    /// Off-diagonal term index.
    OffDiagonal(u32),
}
impl Operator {
    /// Construct an identity position.
    pub fn identity() -> Self {
        Self::Identity
    }
    /// Construct a diagonal position.
    pub fn diagonal(index: u32) -> Self {
        Self::Diagonal(index)
    }
    /// Construct an off-diagonal position.
    pub fn off_diagonal(index: u32) -> Self {
        Self::OffDiagonal(index)
    }
    /// Return the optional term index.
    pub fn term_index(self) -> Option<usize> {
        match self {
            Self::Identity => None,
            Self::Diagonal(index) | Self::OffDiagonal(index) => Some(index as usize),
        }
    }
    /// Return the stored kind.
    pub fn kind(self) -> OperatorKind {
        match self {
            Self::Identity => OperatorKind::Identity,
            Self::Diagonal(_) => OperatorKind::Diagonal,
            Self::OffDiagonal(_) => OperatorKind::OffDiagonal,
        }
    }
}

/// A canonical bit basis state and fixed-length operator string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BasisSseState {
    bits: Vec<BasisBit>,
    operators: Vec<Operator>,
}
impl BasisSseState {
    /// Construct a state after checking non-empty basis and string lengths.
    pub fn new(bits: Vec<BasisBit>, operators: Vec<Operator>) -> Result<Self, SseModelError> {
        if bits.is_empty() || operators.is_empty() {
            return Err(SseModelError::InvalidLength {
                expected: 1,
                actual: 0,
            });
        }
        Ok(Self { bits, operators })
    }
    /// Borrow canonical basis bits.
    pub fn bits(&self) -> &[BasisBit] {
        &self.bits
    }
    /// Mutably borrow canonical basis bits for a validated boundary-state update.
    pub(crate) fn bits_mut(&mut self) -> &mut [BasisBit] {
        &mut self.bits
    }
    /// Borrow the padded operator string.
    pub fn operator_string(&self) -> &[Operator] {
        &self.operators
    }
    /// Mutably borrow the padded string for a validated sampler update.
    pub(crate) fn operator_string_mut(&mut self) -> &mut [Operator] {
        &mut self.operators
    }
    /// Return the current expansion order.
    pub fn expansion_order(&self) -> usize {
        self.operators
            .iter()
            .filter(|operator| !matches!(operator, Operator::Identity))
            .count()
    }
    /// Grow the padded cutoff while preserving all existing operators.
    pub fn grow_operator_string(&mut self, new_length: usize) {
        if new_length > self.operators.len() {
            self.operators.resize(new_length, Operator::identity());
        }
    }
    /// Propagation result including the logarithmic path weight.
    pub fn propagate<M: SseModel>(&self, model: &M) -> Result<PropagationResult, SseModelError> {
        if self.bits.len() != model.num_sites() {
            return Err(SseModelError::InvalidLength {
                expected: model.num_sites(),
                actual: self.bits.len(),
            });
        }
        let mut final_state = self.bits.clone();
        let mut log_weight = 0.0;
        let mut operator_count = 0;
        for operator in &self.operators {
            let Some(index) = operator.term_index() else {
                continue;
            };
            operator_count += 1;
            if operator.kind() != model.operator_kind(index)? {
                return Err(SseModelError::InvalidOperatorKind);
            }
            let value = model.matrix_element(index, &final_state)?;
            if value <= 0.0 {
                return Err(SseModelError::NonPositiveMatrixElement {
                    term_index: index,
                    value,
                });
            }
            log_weight += value.ln();
            if operator.kind() == OperatorKind::OffDiagonal {
                model.apply_off_diagonal(index, &mut final_state)?;
            }
        }
        Ok(PropagationResult {
            final_state: final_state.clone(),
            log_weight,
            trace_closed: final_state == self.bits,
            operator_count,
        })
    }
    /// Validate trace closure after propagation.
    pub fn validate_trace<M: SseModel>(&self, model: &M) -> Result<(), SseModelError> {
        if self.propagate(model)?.trace_closed {
            Ok(())
        } else {
            Err(SseModelError::TraceNotClosed)
        }
    }
}

/// Result of propagating a padded string around imaginary time.
#[derive(Clone, Debug, PartialEq)]
pub struct PropagationResult {
    /// Final canonical bit state.
    pub final_state: Vec<BasisBit>,
    /// Sum of logarithmic local matrix elements.
    pub log_weight: f64,
    /// Whether the trace closed.
    pub trace_closed: bool,
    /// Number of non-identity vertices encountered.
    pub operator_count: usize,
}
