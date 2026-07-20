# qslib scientific conventions

Specification identifier: `qslib-conventions-v1`

Status: normative draft for the first qslib implementation

## 1. Purpose and authority

This document defines the scientific and data conventions shared by qslib
components. It is the authority for geometry, basis states, Hamiltonians,
symmetries, observables, exact calculations, stochastic estimators, TDVP, SSE,
serialization, and interoperability.

An implementation is conforming only when it follows these conventions or
requires the caller to select a separately named compatibility convention.
Legacy behavior must never be selected silently.

The words **must**, **must not**, **should**, **should not**, and **may** are
normative. A public API may use different Rust type names, but it must preserve
the meanings defined here.

## 2. General principles

1. Physical definitions are independent of simulation algorithms.
2. A basis representation is not itself a physical spin or occupation value.
3. Site ordering, bit ordering, boundary behavior, and normalization are part
   of the scientific input and reproducibility record.
4. Total quantities and densities must have distinct names.
5. Algorithmic shifts, regularizers, truncations, and gauges must not silently
   alter the physical Hamiltonian.
6. Invalid dimensions, indices, coefficients, or states must produce errors.
   They must not be repaired by wrapping, clipping, or dropping data unless the
   requested operation explicitly defines that behavior.

## 3. Units and scalar conventions

### 3.1 Natural units

qslib uses

\[
\hbar = 1, \qquad k_B = 1.
\]

Hamiltonian coefficients have units of energy. Time has units of inverse
energy, inverse temperature is

\[
\beta = \frac{1}{T},
\]

and frequencies specified as Hamiltonian coefficients are angular frequencies.
An API accepting dimensional input must require one internally consistent unit
system. qslib does not infer unit conversions from numerical magnitudes.

### 3.2 Real and complex values

Built-in TFIM, J1-J2, and Rydberg Hamiltonian coefficients are real. General
Hamiltonian terms may be complex only when Hermiticity requirements are
explicitly satisfied or the operator is explicitly declared non-Hermitian.

Reference calculations should use `f64` and `Complex<f64>`. Lower precision may
be used for accelerated execution, but the dtype must be recorded and
accumulations that affect scientific results should use at least `f64` when the
backend supports it.

NaN and infinity are invalid configuration values. Runtime non-finite values
must be reported as numerical failures with context.

### 3.3 Angles and phases

Angles are measured in radians. Phases are equivalent modulo \(2\pi\), but
stored log-wavefunction phases are not required to be reduced to a principal
interval.

## 4. Identifiers and collection sizes

Sites are numbered by a zero-based `SiteId`. The first valid site is `0`, and a
system with \(N\) sites contains exactly

\[
0, 1, \ldots, N-1.
\]

The intended Rust representation is a checked newtype over `u32`. Collection
lengths and memory indices use `usize`. Conversion from `usize` to `SiteId`
must be checked.

Bond and term identifiers are also zero-based. An identifier is stable only
within the owning validated object unless a serialized schema explicitly
guarantees broader stability.

Simulation geometries must contain at least one site. Arithmetic used to
derive site counts, Hilbert-space dimensions, or allocation sizes must be
checked for overflow.

## 5. Lattice geometry and site order

### 5.1 Coordinate axes

Two-dimensional coordinates are written `(x, y)`. `x` increases horizontally
and `y` increases vertically in mathematical descriptions. Display software
may invert the visual y direction, but that must not change site identifiers or
physical coordinates.

### 5.2 Canonical rectangular indexing

Rectangular lattices use row-major indexing:

\[
i(x,y) = x + L_x y,
\]

with

\[
0 \le x < L_x, \qquad 0 \le y < L_y.
\]

The inverse mapping is

\[
x = i \bmod L_x, \qquad
y = \left\lfloor\frac{i}{L_x}\right\rfloor.
\]

Thus a dense array with logical shape `[Ly, Lx]` has site `i` at `[y, x]` and
its C-order flattening is the canonical site vector.

The expression `x * Ly + y` is called **legacy x-major indexing**. It is not a
qslib site order and must be handled through an explicit permutation adapter.

### 5.3 Chains

A chain of length \(L\) is embedded along the positive x axis:

\[
\mathbf r_i = (i,0).
\]

It uses site order `0, 1, ..., L-1`.

### 5.4 Square and rectangular embeddings

The default rectangular embedding has unit lattice vectors

\[
\mathbf a_1=(1,0), \qquad \mathbf a_2=(0,1),
\]

and

\[
\mathbf r(x,y)=x\mathbf a_1+y\mathbf a_2.
\]

Non-unit spacing must be explicit in the geometry.

### 5.5 Triangular embedding

The canonical triangular embedding uses integer cell coordinates with

\[
\mathbf a_1=(1,0), \qquad
\mathbf a_2=\left(\frac12,\frac{\sqrt3}{2}\right).
\]

Site indexing remains `x + Lx*y`. Geometry and site order are separate from the
real-space embedding.

### 5.6 Custom geometry

A custom geometry is an ordered collection of finite Cartesian coordinates.
Coordinate entry `i` belongs to `SiteId(i)`. Duplicate coordinates are invalid
for models containing singular distance-dependent interactions, but may be
accepted by purely graph-defined models when explicitly allowed.

### 5.7 Physical order versus scan order

Autoregressive scan order, patch order, memory tiling, and distributed
partitioning do not change physical site identifiers. Any non-canonical scan
order must be represented by an explicit bijection between sequence positions
and canonical `SiteId` values.

## 6. Boundaries, displacement, and distance

### 6.1 Boundary names

Canonical boundary values are `open` and `periodic`. Abbreviations such as
`obc` and `pbc` are accepted only by named compatibility parsers.

Each independent lattice direction may have its own boundary condition.

### 6.2 Open boundaries

Open boundaries do not wrap coordinates. A neighbour step leaving the finite
domain produces no bond.

### 6.3 Periodic boundaries

Periodic boundaries identify coordinates separated by integer simulation-cell
translations. Wrapping an integer coordinate uses Euclidean modulo and always
returns a coordinate in the half-open interval beginning at zero.

### 6.4 Minimum-image displacement

The displacement from site `i` to site `j` is

\[
\Delta\mathbf r_{ij}=\mathbf r_j-\mathbf r_i
\]

reduced to a minimum-length periodic image when periodic directions exist.

For an orthogonal direction of length \(L\), the canonical scalar reduction is

\[
d' = d - L\left\lfloor\frac{d+L/2}{L}\right\rfloor,
\]

which lies in \([-L/2,L/2)\).

For a skew periodic cell, the selected image must minimize Euclidean length
over valid cell translations. If multiple images have exactly equal length,
the integer cell-translation vector is chosen in lexicographic order. This tie
rule makes oriented displacements deterministic.

Distance is the Euclidean norm of the canonical displacement. Squared distance
should be used when a square root is unnecessary.

### 6.5 Long-range periodic interactions

The default periodic long-range interaction uses one minimum-image separation
for each unordered site pair. It is not an Ewald sum and does not include
multiple periodic images. Ewald, cutoff, image-sum, or experimentally supplied
interaction matrices must be separate named policies.

## 7. Bonds, pair sets, and interaction graphs

An undirected bond has two distinct endpoints and is canonically represented as

\[
(\min(i,j),\max(i,j)).
\]

A bond defines topology. A **weighted interaction** additionally associates a
coefficient and operator channel with that bond. Couplings belong to weighted
interactions, not to the bond set as a whole. The canonical resolved form is

\[
(i,j,\kappa_{ij},\text{channel}),
\]

where `i < j`, `kappa_ij` is finite, and `channel` identifies the operator
acted on by the coefficient. For example, an isotropic Heisenberg interaction
has channel `heisenberg_exchange` and coefficient \(J_{ij}\).

Distinct channels may act on the same unordered pair. Multiple contributions
to the same pair and channel are also allowed when they represent distinct
periodic images or named physical terms. A canonicalization operation may sum
algebraically identical contributions, but it must preserve enough provenance
to reconstruct how the resolved coefficient was obtained.

An interaction identity consists of its canonical endpoints, operator channel,
and, when relevant, periodic-image or named-term identity. Shell membership is
metadata rather than identity unless two shell contributions to the same pair
are intended to remain distinct. Pair-dependent data must be keyed by this
identity, never by incidental insertion or hash-map iteration order.

A canonical simple bond set:

- rejects self-bonds;
- contains each unordered pair at most once;
- is sorted lexicographically when serialized or returned as canonical output.

Tiny periodic lattices may connect the same pair through more than one lattice
direction. qslib distinguishes:

- `simple` multiplicity, which stores the unordered pair once;
- `periodic_images` multiplicity, which stores each physically distinct image
  or directed step with explicit multiplicity.

Built-in model constructors use `simple` multiplicity unless their API or
serialized specification says otherwise. Multiplicity is part of the physical
model and must be included in fingerprints.

A uniform coupling, a coupling per geometric shell, and a coupling per pair
are three input forms for the same resolved weighted-interaction collection.
Uniform and shell-based inputs are convenience constructors. They must be
expanded and validated before Hamiltonian evaluation. No evaluator may assume
that all entries in an interaction collection share one coefficient.

Pair-dependent couplings may be positive, negative, or zero. Canonical output
should omit exactly zero weighted terms from numerical evaluation, while
provenance may retain them when their presence is scientifically meaningful.
Signed zero is not scientifically distinct.

A sparse pair table and a dense coupling matrix are equivalent input forms. A
dense matrix for an undirected scalar channel must have shape `[N,N]`, be
symmetric within a declared validation tolerance, and have a zero diagonal
unless onsite terms are part of that channel. It is resolved from entries with
`i < j`, so symmetric entries are never counted twice. A sparse table must
reject duplicate interaction identities unless its schema explicitly declares
that duplicates are additive.

For quenched disorder, the durable run description must store the realized
coefficient associated with every interaction identity. A distribution name
and random seed are useful provenance but are not substitutes for the realized
coupling table. Generation order must be defined from canonical interaction
order rather than map iteration or worker scheduling.

Distance shells are selected by a target squared distance and a non-negative
tolerance. A shell selector must state whether tolerance is absolute or
relative. The default is absolute tolerance.

## 8. Binary basis states

### 8.1 Canonical bit meaning

Each two-level site is stored as a bit \(b_i\in\{0,1\}\). In a simulation basis
associated with Pauli axis \(a\),

\[
\sigma_i^a |b_i\rangle = (1-2b_i)|b_i\rangle.
\]

Therefore:

| Bit | Axis eigenvalue | Canonical label |
| ---: | ---: | --- |
| `0` | `+1` | `Plus` |
| `1` | `-1` | `Minus` |

APIs should prefer `Plus` and `Minus` over `Up` and `Down`, because the latter
are used inconsistently across physics communities.

The Pauli eigenvalue conversion is

\[
z(b)=1-2b=(-1)^b.
\]

### 8.2 Rydberg occupation

The canonical Rydberg occupation is

\[
n_i=b_i=\frac{1-\sigma_i^z}{2}.
\]

Thus bit `0` is the unoccupied ground state and bit `1` is the occupied Rydberg
state. This convention is intentionally separate from the names `Plus` and
`Minus`.

An alternate convention \(n=(1+Z)/2\) must be represented by an explicit basis
or compatibility transform. It must not reuse the same serialized state label.

### 8.3 Dense state layout

A batch of basis states has logical shape `[batch, site]`. The site dimension is
in canonical site order. Accepted scalar storage types may include booleans and
unsigned integers, but every value must be exactly zero or one.

### 8.4 Bit-packed state and endianness

Bit-packed states are little-endian by site identifier:

\[
m(b)=\sum_{i=0}^{N-1} b_i 2^i.
\]

Site `0` is the least significant bit. Full computational-basis state vectors
use mask `m` as the vector index. The basis order is therefore

\[
|00\ldots0\rangle, |10\ldots0\rangle,
|01\ldots0\rangle, \ldots
\]

when site `0` is written first inside the ket.

For more sites than fit in one machine word, words are ordered from least to
most significant and each word is itself little-endian by site identifier.
Serialized packed states must record the word width.

### 8.5 Fixed-occupation sectors

The Hamming weight is

\[
K(b)=\sum_i b_i.
\]

A fixed sector contains states with one declared value of `K`. Sector basis
states are ordered by increasing packed integer value unless another order is
explicitly named. A Hamiltonian may use this sector only if it conserves
Hamming weight in the selected simulation basis.

## 9. Physical axes and simulation bases

Physical operators are expressed in a fixed physical frame with axes `x`, `y`,
and `z`. A simulation basis states which physical Pauli operator is diagonal in
the stored bits.

In simulation basis `z`, bit values encode physical \(Z\) eigenstates. In
simulation basis `x`, they encode physical \(X\) eigenstates. A basis change
does not relabel the physical Hamiltonian or observables.

The Hadamard rotation satisfies

\[
H Z H = X, \qquad H X H = Z.
\]

Consequently, the physical TFIM

\[
H_{\mathrm{TFIM}}=-J\sum_{\langle i,j\rangle}Z_iZ_j-h\sum_iX_i
\]

is represented in simulation basis `x` as

\[
H_{\mathrm{TFIM},x}=-J\sum_{\langle i,j\rangle}X_iX_j-h\sum_iZ_i.
\]

Initial-state specifications must distinguish physical polarization from the
simulation basis. For example, physical `PlusZ` is the all-zero bitstring only
in simulation basis `z`; in simulation basis `x` it is a superposition.

## 10. Wavefunctions and state vectors

A wavefunction assigns a complex amplitude \(\psi(b)\) to each basis state.
Unless explicitly normalized,

\[
Z_\psi=\sum_b |\psi(b)|^2
\]

need not equal one. Sampling probability is

\[
p(b)=\frac{|\psi(b)|^2}{Z_\psi}.
\]

The canonical complex log wavefunction is

\[
\log\psi(b)=\log|\psi(b)|+i\arg\psi(b).
\]

Its imaginary part may be unwrapped. Ratios are evaluated as

\[
\frac{\psi(b')}{\psi(b)}=
\exp\left(\log\psi(b')-\log\psi(b)\right).
\]

A zero amplitude has real log amplitude `-infinity`. Algorithms must handle
such states deliberately and must not create NaN through an unguarded
`-infinity - -infinity` operation.

Global complex phase does not affect physical observables. State-vector
comparisons that are intended to be phase insensitive must align or eliminate
global phase explicitly.

## 11. Pauli and spin operators

The Pauli matrices are

\[
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
Y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

Spin-half operators use

\[
S^a=\frac12\sigma^a.
\]

In simulation basis `z`, `X_i` flips bit `i`, while `Z_i` multiplies a state by
`1 - 2*b_i`. A product of `X` operators flips every listed site exactly once.
Repeated Pauli factors on the same site must be algebraically reduced rather
than represented as duplicate flips.

Operator support is an ordered or canonical set of distinct `SiteId` values as
required by the operator. Support arity is validated.

## 12. Hamiltonian conventions

### 12.1 General term representation

A physical Hamiltonian is

\[
H=cI+\sum_a c_a O_a.
\]

The scalar constant `c` is part of the physical energy. Algorithmic shifts used
for sampling or solving are stored separately.

Term order has no physical meaning. Duplicate algebraically identical terms
may be combined, but the combination procedure must be deterministic.

Every resolved term carries its own coefficient \(c_a\). Writing one coupling
outside a sum is only shorthand for assigning that value to every term in the
sum. Public model representations and Hamiltonian evaluators must support
term-dependent coefficients. A representation may factor a common numerical
scale for storage or performance, but its resolved physical coefficient is the
product of that scale and the term weight and must be used in serialization,
fingerprinting, validation, and cross-backend comparisons.

### 12.2 Transverse-field Ising model

The general pair-dependent TFIM is

\[
H=-\sum_{(i,j)\in E}J_{ij}Z_iZ_j-\sum_i h_iX_i.
\]

The homogeneous notation \(-J\sum Z_iZ_j-h\sum X_i\) is a constructor shorthand
for `J_ij=J` and `h_i=h`. Each weighted bond contributes once, including any
explicit multiplicity. Positive `J_ij` is ferromagnetic in this convention.
qslib model construction permits any finite real `J_ij` and `h_i`; a particular
algorithm may impose additional sign restrictions.

In simulation basis `z`:

- diagonal contribution from pair `(i,j)`: \(-J_{ij} z_i z_j\);
- connected state for each field term: flip site `i` with matrix element
  \(-h_i\).

In simulation basis `x`:

- diagonal contribution: \(-h_i z_i\), where `z_i=1-2*b_i` now labels the
  physical X eigenvalue;
- connected state for each bond: flip sites `i` and `j` with matrix element
  \(-J_{ij}\).

### 12.3 Heisenberg exchange and the J1-J2 model

The general isotropic pair-dependent Heisenberg exchange model is

\[
H=\sum_{(i,j)\in E}J_{ij}\,\mathbf S_i\cdot\mathbf S_j.
\]

Positive `J_ij` is antiferromagnetic and negative `J_ij` is ferromagnetic.
With \(\hbar=1\) and simulation basis `z`, one bond has

\[
\langle b|H_{ij}|b\rangle=\frac{J_{ij}}{4}z_i z_j.
\]

If the two bits are different, the pair-flipped state has matrix element

\[
\langle b^{ij}|H_{ij}|b\rangle=\frac{J_{ij}}{2}.
\]

Parallel bits have no off-diagonal connection from that bond. This Hamiltonian
conserves Hamming weight in simulation basis `z`.

The homogeneous J1-J2 model is the specialization

\[
H=J_1\sum_{(i,j)\in E_1}\mathbf S_i\cdot\mathbf S_j
 +J_2\sum_{(i,j)\in E_2}\mathbf S_i\cdot\mathbf S_j,
\]

expanded to `J_ij=J1` on `E1` and `J_ij=J2` on `E2`. The sets `E1` and `E2`
must state their geometric meaning. For the unit square lattice, `E1` is the
axial nearest-neighbour shell and `E2` is the diagonal next-nearest-neighbour
shell.

A disordered J1-J2 model may instead assign a distinct realized coupling to
each member of either shell, for example \(J_{ij}^{(1)}\) on `E1` and
\(J_{ij}^{(2)}\) on `E2`. This remains an isotropic Heisenberg exchange model
with J1-J2 geometry, but it is not the homogeneous two-parameter J1-J2 model.
Its input must say whether values are supplied directly, drawn independently,
or generated with spatial correlations.

Pair-dependent values constitute bond disorder. They produce frustration only
when the preferred local bond constraints cannot all be satisfied, as can
happen through mixed signs, competing interaction shells, or lattice geometry.
Documentation and result metadata should not use `disordered` and `frustrated`
as synonyms.

If the same pair occurs in more than one shell or periodic image, the physical
contributions are additive unless the model explicitly treats the images as
separate channels. The resolved representation must make this choice visible.

### 12.4 Rydberg model

The canonical driven Rydberg Hamiltonian is

\[
H=-\frac{\Omega}{2}\sum_i X_i
  -\Delta\sum_i n_i
  +\sum_{i<j}V_{ij}n_i n_j,
\qquad n_i=b_i=\frac{1-Z_i}{2}.
\]

For van der Waals interactions,

\[
V_{ij}=\frac{C_6}{r_{ij}^6}.
\]

An interaction matrix must be symmetric, finite, and have a zero diagonal
unless a separately defined onsite term is intended. When the complete matrix
is used, the diagonal energy is evaluated either as

\[
\sum_{i<j}V_{ij}n_i n_j
\]

or equivalently

\[
\frac12\sum_{i,j}V_{ij}n_i n_j.
\]

Each site flip has matrix element \(-\Omega/2\). `Omega`, `Delta`, and `C6` are
finite real values. Algorithms may restrict their signs, but the physical
model does not silently take absolute values.

### 12.5 Basis-aware construction

A model specification describes a physical Hamiltonian. Conversion to a
simulation basis produces an equivalent operator representation. The original
physical specification and selected simulation basis must both remain
available for provenance and observable evaluation.

## 13. Hamiltonian action and local energy

For a basis state `b`, a Hamiltonian action consists of:

- its diagonal matrix element \(H_{bb}\);
- connected states \(b'\ne b\) with nonzero matrix elements \(H_{bb'}\).

Connected states must be unique after applying all flips. Contributions from
different terms leading to the same state may be returned separately or
combined, but the selected behavior must be explicit.

The variational local energy is

\[
E_{\mathrm{loc}}(b)=
\frac{\langle b|H|\psi\rangle}{\langle b|\psi\rangle}
=H_{bb}+\sum_{b'\ne b}H_{bb'}\frac{\psi(b')}{\psi(b)}.
\]

The variational energy is

\[
E=\mathbb E_{p(b)}[E_{\mathrm{loc}}(b)].
\]

For a Hermitian Hamiltonian the exact expectation is real. A measured imaginary
component must be retained as a diagnostic rather than silently discarded.

Energy variance is

\[
\operatorname{Var}(H)=
\mathbb E_p\left[|E_{\mathrm{loc}}-E|^2\right].
\]

`energy_variance_density` means `Var(H)/N`. A quantity divided by `N^2` must use
a different explicit name.

## 14. Symmetry conventions

### 14.1 Site permutations

A site permutation stores `source_for_destination`:

\[
(g\cdot b)_j=b_{p_g(j)}.
\]

This definition matches gather-style array indexing. APIs must not describe the
same array ambiguously as both `old_to_new` and `new_to_old`. Conversion to a
destination-for-source representation requires an explicit inverse.

A permutation must be a bijection over `0..N`. Composition and inversion must
be tested against the action above.

### 14.2 Lattice transformations

A translation by `(dx, dy)` maps physical coordinates to

\[
(x,y)\mapsto(x+dx,y+dy)
\]

with periodic wrapping where defined. Point-group operations act on physical
coordinates and are then converted to canonical site permutations.

### 14.3 Group projection

For a finite group `G` and one-dimensional unitary character `chi`, the
projected wavefunction convention is

\[
(P_\chi\psi)(b)=\frac{1}{|G|}
\sum_{g\in G}\chi(g)^*\psi(g^{-1}\cdot b).
\]

The trivial representation has `chi(g)=1`. Any implementation using an
equivalent change of summation variable must demonstrate parity with this
definition.

### 14.4 Spin inversion

Global spin inversion maps every bit as

\[
b_i\mapsto1-b_i.
\]

It is a symmetry only when it commutes with the physical Hamiltonian in the
selected parameters and basis.

### 14.5 Diagonal sign and phase gauges

A diagonal gauge acts as

\[
\psi'(b)=e^{i\phi(b)}\psi(b).
\]

The standard linear gauge form is

\[
\phi(b)=\pi\sum_i a_i b_i.
\]

Integer `a_i` values produce signs. Non-integer values produce general phases
and must not be described as a sign-only Marshall transformation. A gauge
changes the representation of amplitudes and operator matrix elements but not
physical predictions when applied consistently.

## 15. Observable normalization

### 15.1 Energy

`energy` is total energy \(\langle H\rangle\). `energy_density` and
`energy_per_site` both mean

\[
e=\frac{\langle H\rangle}{N}.
\]

Public schemas should choose one canonical field name and treat the other as a
documented alias only at compatibility boundaries.

### 15.2 Magnetization

Total and per-site magnetization are

\[
M_a=\sum_i\sigma_i^a, \qquad m_a=\frac{M_a}{N}.
\]

A field named `sigma_a` in aggregate output means \(m_a\), not \(M_a\).

### 15.3 Correlations

The two-site Pauli correlation is

\[
C_{aa}(i,j)=\langle\sigma_i^a\sigma_j^a\rangle.
\]

The connected correlation is

\[
C^{\mathrm c}_{aa}(i,j)=C_{aa}(i,j)
-\langle\sigma_i^a\rangle\langle\sigma_j^a\rangle.
\]

A shell or distance average is the arithmetic mean over the explicitly
selected pair set. Ordered and unordered pair averages must not share a field
name unless their normalization makes them identical.

### 15.4 Structure factor

The raw static structure factor is

\[
S_{aa}(\mathbf q)=\frac1N\sum_{i,j}
e^{i\mathbf q\cdot(\mathbf r_i-\mathbf r_j)}
\langle\sigma_i^a\sigma_j^a\rangle.
\]

The connected structure factor replaces the correlation with its connected
form and must be named accordingly. Wavevectors are Cartesian and satisfy that
`q dot r` is dimensionless.

### 15.5 Sublattice order

For real weights `w_i`,

\[
m_w=\frac1N\sum_i w_i\sigma_i^z.
\]

Outputs must distinguish \(\langle m_w\rangle\),
\(\langle|m_w|\rangle\), and \(\langle m_w^2\rangle\). They are not
interchangeable order parameters.

### 15.6 Quantum Fisher information

For the pure-state generator

\[
G=\frac12\sum_i w_i\sigma_i^z,
\]

qslib defines

\[
F_Q=4\operatorname{Var}(G)
=\operatorname{Var}\left(\sum_i w_i\sigma_i^z\right),
\qquad f_Q=\frac{F_Q}{N}.
\]

The factor of four is therefore already absorbed by using Pauli values in the
last expression. APIs using a different generator must report it explicitly.

### 15.7 Thermal observables

With `kB=1`, heat capacity is

\[
C=\beta^2\left(\langle H^2\rangle-\langle H\rangle^2\right),
\qquad c=C/N.
\]

## 16. Exact bases and exact solvers

A full two-level Hilbert space has dimension \(2^N\). A fixed-Hamming-weight
sector has dimension

\[
\binom NK.
\]

Exact basis enumeration follows increasing little-endian packed-state value.
Sparse matrices use basis indices, not raw packed masks, when a restricted
sector is active.

Eigenvalues are returned in ascending algebraic order unless the request names
another ordering. Eigenvectors are columns and correspond one-to-one with the
returned eigenvalues.

Degenerate eigenspaces do not have unique eigenvectors. Tests must compare the
invariant subspace or projector rather than individual vector phases and bases.

Real-time exact evolution uses

\[
|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle.
\]

Imaginary-time evolution uses \(e^{-\tau H}\) followed by explicit
normalization when a normalized state is required.

## 17. Stochastic samples and statistics

A sample batch has shape `[sample, site]`. Sample count refers to the leading
dimension after all distributed partitions are combined.

Means use normalized sample weights. Unweighted samples have equal weight.
Exact enumeration weights must sum to one after normalization.

An IID standard error for real values is

\[
\mathrm{SE}=\frac{s}{\sqrt B},
\]

where `s` is the sample standard deviation with denominator `B-1`. This formula
must not be reported for correlated Markov-chain samples without an effective
sample-size correction.

Complex estimators retain real and imaginary sample statistics. A real-only
summary may be produced only when the discarded imaginary component remains
available as a diagnostic.

Independent-chain aggregation must distinguish within-chain Monte Carlo error
from between-chain variation. R-hat, effective sample size, and integrated
autocorrelation time must name the algorithm used because multiple estimators
exist.

For quenched disorder, expectation at fixed realization and the disorder
ensemble average are distinct operations. Using `r` for a realization,

\[
\overline{\langle A\rangle}
=\sum_r w_r\langle A\rangle_r,
\qquad \sum_r w_r=1.
\]

Outputs must retain the realization identifier and fixed-realization estimate
before ensemble reduction. Uncertainty from stochastic sampling within each
realization and uncertainty from the finite realization ensemble must be
reported separately. A run with one realization has no empirical estimate of
disorder-ensemble uncertainty.

## 18. Variational geometry and TDVP

For real parameters \(\theta_k\), define complex log derivatives

\[
O_k(b)=\frac{\partial\log\psi(b)}{\partial\theta_k},
\qquad \Delta O_k=O_k-\langle O_k\rangle.
\]

The real-parameter quantum geometric tensor is

\[
S_{kl}=\operatorname{Re}\left\langle
\Delta O_k^*\Delta O_l
\right\rangle.
\]

Define the complex energy covariance

\[
F_k=\left\langle\Delta O_k^*
\left(E_{\mathrm{loc}}-\langle E_{\mathrm{loc}}\rangle\right)
\right\rangle.
\]

The canonical equations are

\[
S\dot\theta=\operatorname{Im}F
\quad\text{for real time},
\]

and

\[
S\dot\theta=-\operatorname{Re}F
\quad\text{for imaginary time}.
\]

Empirical covariances use denominator `B`, matching expectation-value
estimators. Statistical uncertainty estimates may use `B-1` where appropriate.

Regularization is part of the numerical solve, not the physical QGT. A result
must record the requested and effective backend, regularization method,
regularization magnitude, spectral cutoff, convergence tolerance, iteration
count, and residual.

For a solved direction `v` and mode-specific right-hand side `f`, the projected
residual diagnostic is

\[
r^2=\operatorname{Var}(H)+v^T S v-2v^T f.
\]

The normalized residual is `r2 / Var(H)` when variance is nonzero. Numerical
flooring used at zero variance must be recorded.

Parameter flattening order must be deterministic and fingerprinted. Restoring a
checkpoint with a different parameter-layout fingerprint is an error.

## 19. Time integration

An integrator state records at least physical time, parameters, proposed next
step, accepted-step count, and rejected-step count.

Adaptive attempts are transactional. A rejected attempt must not change
physical time, accepted parameters, accepted trajectory, or the random streams
assigned to accepted time nodes.

Step-error norms must state their metric and normalization. A QGT norm uses

\[
\|\delta\|_S=\sqrt{\delta^T S\delta}.
\]

If this value is divided by parameter count or another scale, the output and
tolerance specification must say so explicitly.

Trajectory observations are associated with accepted states. Attempt logs may
include rejected states but must not present them as physical trajectory
points.

## 20. SSE decomposition conventions

An SSE decomposition may write the physical Hamiltonian as

\[
H=E_{\mathrm{shift}}-\sum_a B_a,
\]

where sampled matrix elements of `B_a` are non-negative. `E_shift` is an
algorithmic decomposition value derived from the physical Hamiltonian. It does
not modify that Hamiltonian.

For expansion order `n`, inverse temperature `beta`, and this sign convention,

\[
E=E_{\mathrm{shift}}-\frac{\langle n\rangle}{\beta},
\]

and

\[
C=\langle n^2\rangle-\langle n\rangle^2-\langle n\rangle.
\]

Reported physical results must restore the shift. SSE-specific operator kinds,
padding identities, and update vertices must not be exposed as physical
Hamiltonian operators.

## 21. Determinism and numerical validation

Canonical collections such as bond sets, site permutations, sectors, and
connected-state lists must have deterministic ordering.

A seed alone is not a complete reproducibility contract. Reproducible output
must identify the random-number algorithm, qslib version, convention version,
threading or distributed policy where it affects reduction order, and relevant
floating-point dtype.

Parallel scheduling must not change per-chain random streams. Chain or job
seeds should be derived from a master seed and stable logical index rather than
worker identity.

Floating-point equality tolerances belong to the tested quantity. qslib must
not define one global epsilon for geometry, Hermiticity, solver convergence,
and physical parity tests.

Reference parity tests should use `f64`, deterministic inputs, and small systems
whose full matrices or state spaces can be independently enumerated.

## 22. Serialization and schema evolution

Every durable scientific document must include a schema version. Documents
governed by this specification include

```text
convention_schema: qslib-conventions-v1
```

Canonical serialized enum values use lowercase snake case. Examples include
`open`, `periodic`, `row_major`, `z`, `real_time`, and `imaginary_time`.

Versioned input schemas should reject unknown fields by default. Compatibility
loaders may accept legacy fields, but they must resolve them into a canonical
document before simulation.

Serialized multidimensional arrays must record shape, dtype, and flattening
order. Native Rust, NumPy, Torch, or BLAS memory layout must not be inferred
from an unlabelled byte sequence.

Serialized complex JSON scalars use an object with named real and imaginary
parts:

```json
{"re": 1.0, "im": -0.25}
```

Binary formats may use interleaved or split complex storage only when the
format metadata identifies it.

Scientific fingerprints must include all convention-sensitive fields,
including site order, bit convention, simulation basis, boundary policy, bond
multiplicity, every resolved Hamiltonian coefficient, and convention schema.
For disordered models they must also include the canonical interaction
identities and realized coupling table. Distribution parameters, correlation
rules, and random-generator details must be stored as provenance when qslib
generated that table.

A breaking change to any normative mapping or normalization requires a new
convention schema identifier.

## 23. Interoperability with existing projects

### 23.1 ncli dynamics and modern row-major paths

The modern ncli dynamics path uses `site = x + Lx*y`, little-endian exact basis
masks, `bit 0 -> +1`, `bit 1 -> -1`, and Rydberg `occupation = bit`. These map
directly to qslib after validation of Hamiltonian signs and boundary policy.

### 23.2 ncli legacy x-major paths

Legacy J1-J2 and symmetry paths may use

\[
i_{\mathrm{legacy}}=xL_y+y.
\]

They require an explicit permutation to qslib row-major order. Stripe labels
and directional observables must be transformed with the state. Merely
reinterpreting the same flat array is invalid.

### 23.3 Existing Rust SSE states

The existing SSE crate uses row-major geometry, but its `Spin` type combines
two conventions:

- TFIM interprets `Up` as Pauli-Z `+1`;
- Rydberg interprets `Up` as occupied.

Under qslib, occupied Rydberg state has Pauli-Z value `-1`. Therefore there is
no universal raw enum conversion from legacy SSE `Spin` to qslib bits.

Adapters must be model aware:

| Meaning | qslib bit | legacy SSE state |
| --- | ---: | --- |
| Pauli-Z `+1` | `0` | `Up` |
| Pauli-Z `-1` | `1` | `Down` |
| Rydberg unoccupied | `0` | `Down` |
| Rydberg occupied | `1` | `Up` |

New shared code must use qslib basis and occupation types rather than extending
the ambiguous legacy enum.

## 24. Required conformance vectors

Implementations must include equivalent automated tests for the following
vectors.

### 24.1 Rectangular indexing

For `Lx=3`, `Ly=2`:

| `(x,y)` | Site |
| --- | ---: |
| `(0,0)` | `0` |
| `(1,0)` | `1` |
| `(2,0)` | `2` |
| `(0,1)` | `3` |
| `(1,1)` | `4` |
| `(2,1)` | `5` |

The open square-lattice nearest-neighbour simple bond set is

```text
(0,1), (0,3), (1,2), (1,4), (2,5), (3,4), (4,5)
```

The periodic simple bond set is

```text
(0,1), (0,2), (0,3), (1,2), (1,4),
(2,5), (3,4), (3,5), (4,5)
```

### 24.2 Bit packing and local values

For site-ordered bits `[1, 0, 1, 1]`:

```text
packed mask      = 13
Pauli-Z values   = [-1, +1, -1, -1]
Rydberg occupancy = [1, 0, 1, 1]
Hamming weight   = 3
```

### 24.3 Two-site TFIM

For one bond `(0,1)` in simulation basis `z`:

- state `00` has diagonal energy `-J`;
- state `01` has diagonal energy `+J`;
- every state connects to its site-0 flip with matrix element `-h`;
- every state connects to its site-1 flip with matrix element `-h`.

### 24.4 One Heisenberg bond

For one bond of strength `J`:

- `00` and `11` have diagonal energy `+J/4` and no pair-flip connection;
- `01` and `10` have diagonal energy `-J/4`;
- `01` connects to `10`, and `10` connects to `01`, with matrix element `J/2`.

### 24.5 Pair-dependent Heisenberg couplings

For three sites with weighted bonds `J_01=2` and `J_12=-3`, state `010` has
Pauli values `[+1, -1, +1]`. Its diagonal exchange energy is

\[
\frac{2}{4}(+1)(-1)+\frac{-3}{4}(-1)(+1)=\frac14.
\]

It connects to `100` by flipping pair `(0,1)` with matrix element `+1`, and it
connects to `001` by flipping pair `(1,2)` with matrix element `-3/2`. A
backend that factors out one global exchange coupling cannot pass this vector.

### 24.6 Two-site Rydberg model

For interaction `V`:

| State | Diagonal energy |
| --- | ---: |
| `00` | `0` |
| `01` | `-Delta` |
| `10` | `-Delta` |
| `11` | `-2*Delta + V` |

Every single-site flip has matrix element `-Omega/2`.

### 24.7 Observable normalization

For the all-zero product state in simulation basis `z`:

```text
sigma_z per site = +1
raw zz correlation = +1 for every pair
connected zz correlation = 0
uniform Pauli-generator Fisher density = 0
```

### 24.8 Basis rotation parity

For every small TFIM system used as a reference, the spectra of the canonical
`z`-basis and Hadamard-rotated `x`-basis matrices must agree within the declared
`f64` eigensolver tolerance.

## 25. Open decisions for later specifications

The following are deliberately not fixed by this first convention document:

- the public crate and module hierarchy inside the qslib workspace;
- the sparse matrix storage format;
- the default dense and sparse eigensolver implementations;
- a particular random-number generator;
- binary artifact and checkpoint container formats;
- GPU tensor interoperability;
- model-specific Ewald or cutoff conventions;
- the public stability level of experimental algorithms.

Each may be specified separately without changing this document, provided it
does not alter the scientific meanings defined above.
