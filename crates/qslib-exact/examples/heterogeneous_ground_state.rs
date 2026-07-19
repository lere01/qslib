//! Exact ground-state gate for a four-site heterogeneous Ising Hamiltonian.

use qslib_core::{Complex64, Hamiltonian, Pauli, PauliString, SiteCount, SiteId};
use qslib_exact::{ExactBasis, GroundState, diagonalize_hermitian};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let h = Hamiltonian::new_hermitian(
        SiteCount::new(4)?,
        Complex64::new(0.0, 0.0),
        vec![
            (-1.0, [(0, Pauli::Z), (1, Pauli::Z)]),
            (-2.0, [(1, Pauli::Z), (2, Pauli::Z)]),
            (0.5, [(2, Pauli::Z), (3, Pauli::Z)]),
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
                .expect("static example support is valid"),
            )
        })
        .collect(),
    )?;
    let basis = ExactBasis::full(SiteCount::new(4)?);
    let basis = basis?;
    let matrix = qslib_exact::DenseMatrix::from_hamiltonian(&h, &basis)?;
    let ground = GroundState::from_spectrum(&diagonalize_hermitian(&matrix)?)?;
    println!(
        "ground energy = {:.12}, residual = {:.3e}",
        ground.energy(),
        ground.residual()
    );
    Ok(())
}
