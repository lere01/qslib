# Geometry and weighted interactions

This guide describes the physical objects represented by qslib's foundational
API. It is written for users specifying a lattice model; the Rust type names
are included only after the scientific convention is clear.

## Sites and coordinates

Rectangular systems use canonical row-major site order

\[
i(x,y)=x+L_x y.
\]

The first coordinate is horizontal and the second vertical. `Boundary::Open`
does not wrap a neighbour, while `Boundary::Periodic` identifies opposite
faces independently in each direction. `LatticeKind::Triangular` uses the
embedding \(r(x,y)=x(1,0)+y(1/2,\sqrt{3}/2)\); this changes distances, not site
order.

```rust
use qslib_core::{Boundary, RectangularGeometry, SiteId};

let geometry = RectangularGeometry::new(3, 2, Boundary::Open, Boundary::Periodic)?;
assert_eq!(geometry.site_id(2, 1)?, SiteId::new(5));
assert_eq!(geometry.coordinate(SiteId::new(4))?.as_tuple(), (1.0, 1.0));
# Ok::<(), Box<dyn std::error::Error>>(())
```

Periodic minimum-image distances are explicit. A shell selector receives a
squared target and either `ShellTolerance::Absolute(epsilon)` or
`ShellTolerance::Relative(epsilon)`. The relative policy is
\(|d^2-d_0^2|\leq\epsilon |d_0^2|\). Invalid targets and tolerances return a
structured error.

`BondMultiplicity::Simple` returns one endpoint-only bond for each unordered
pair. `BondMultiplicity::PeriodicImages` retains the source site, lattice step,
and cell translation, so two physically distinct images never collapse into an
ambiguous duplicate. Extent-one periodic directions simply produce no self-bond.

## Couplings belong to terms

For a pair interaction, qslib resolves the physical coefficient per term:

\[
(i,j,\kappa_{ij},\mathrm{channel},\mathrm{name}).
\]

The coefficient may be positive, negative, or exactly zero. Zero terms remain
in the declared table for provenance but are excluded by
`InteractionTable::active_interactions()` during numerical evaluation. Dense
coupling matrices are row-major, finite, symmetric within an explicitly chosen
tolerance, and have a zero diagonal. Sparse pairs are canonicalized and reject
self-pairs and duplicate identities.

This makes a frustrated or disordered J1-J2 model representable without a
hidden global `J`: each realized \(J_{ij}\) is stored beside its shell or named
term. `InteractionIdentity::named` keeps overlapping physical contributions
distinct when they should not be algebraically merged.

## Disorder provenance

`InteractionTable::realize_uniform_disorder` samples replacement coefficients
from the requested finite interval. The durable result includes the complete
realized interaction table, site count, master seed, distribution bounds,
coefficient semantics, the `chacha20` algorithm, and the keyed `qslib-seed-v1`
child-seed scheme in the `disorder` domain. Derivation is keyed by canonical
interaction identity and logical realization index, so worker scheduling and
insertion order cannot change a realization.

## Legacy order

Legacy x-major arrays (`x * Ly + y`) are not silently reinterpreted. Use
`XMajorAdapter` to convert both state vectors and interaction endpoints before
comparing them with canonical qslib data.
