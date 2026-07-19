use crate::{BasisError, Real, SiteCount, SiteId};
use std::cmp::Ordering;
use std::fmt::{self, Display, Formatter};

/// Errors raised while constructing or querying a simulation geometry.
#[derive(Clone, Debug, PartialEq)]
pub enum GeometryError {
    /// A lattice extent was zero.
    InvalidExtent {
        /// Extent name.
        axis: &'static str,
        /// Supplied extent.
        value: usize,
    },
    /// A checked site-count or index calculation overflowed.
    DimensionOverflow {
        /// Checked operation that overflowed.
        operation: &'static str,
    },
    /// A coordinate is outside the finite lattice domain.
    CoordinateOutOfRange {
        /// Horizontal coordinate.
        x: isize,
        /// Vertical coordinate.
        y: isize,
    },
    /// A bond was requested from a site to itself.
    SelfBond {
        /// Site used at both endpoints.
        site: SiteId,
    },
    /// A serialized bond is not in canonical endpoint or source form.
    NonCanonicalBond {
        /// First endpoint supplied by the caller.
        first: SiteId,
        /// Second endpoint supplied by the caller.
        second: SiteId,
    },
    /// A site identifier is outside this geometry.
    SiteOutOfRange {
        /// Invalid site identifier.
        site: u32,
        /// Number of sites in the geometry.
        site_count: usize,
    },
    /// A custom geometry contained no sites.
    EmptyCustomGeometry,
    /// A custom coordinate was not finite.
    NonFiniteCoordinate {
        /// Coordinate index when known.
        index: usize,
    },
    /// A shell target or tolerance was invalid.
    InvalidShell {
        /// Invalid target or tolerance value.
        value: Real,
    },
}

impl Display for GeometryError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidExtent { axis, value } => {
                write!(formatter, "{axis} extent {value} is invalid")
            }
            Self::DimensionOverflow { operation } => {
                write!(formatter, "geometry {operation} overflowed")
            }
            Self::CoordinateOutOfRange { x, y } => {
                write!(formatter, "coordinate ({x}, {y}) is outside the geometry")
            }
            Self::SelfBond { site } => {
                write!(formatter, "site {} cannot form a self-bond", site.get())
            }
            Self::NonCanonicalBond { first, second } => write!(
                formatter,
                "bond endpoints {} and {} are not in canonical form",
                first.get(),
                second.get()
            ),
            Self::SiteOutOfRange { site, site_count } => write!(
                formatter,
                "site {site} is outside a {site_count}-site geometry"
            ),
            Self::EmptyCustomGeometry => {
                formatter.write_str("a custom geometry must contain at least one site")
            }
            Self::NonFiniteCoordinate { index } => {
                write!(formatter, "custom coordinate {index} is non-finite")
            }
            Self::InvalidShell { value } => {
                write!(formatter, "shell parameter {value:?} is invalid")
            }
        }
    }
}

impl std::error::Error for GeometryError {}

/// Boundary condition applied independently to a lattice direction.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum Boundary {
    /// Do not wrap a neighbour leaving the finite domain.
    Open,
    /// Identify opposite faces of the finite domain.
    Periodic,
}

/// Multiplicity policy for generated lattice bonds.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum BondMultiplicity {
    /// Keep each unordered endpoint pair once.
    Simple,
    /// Keep every generated periodic image or directed lattice step.
    PeriodicImages,
}

/// Tolerance policy for selecting a squared-distance shell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ShellTolerance {
    /// Accept `|d² - target²| <= epsilon`.
    Absolute(Real),
    /// Accept `|d² - target²| <= epsilon * |target²|`.
    Relative(Real),
}

/// Real-space embedding used by a rectangular coordinate grid.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum LatticeKind {
    /// Orthogonal square or rectangular embedding.
    Square,
    /// The canonical triangular embedding.
    Triangular,
}

/// Explicit adapter for legacy x-major flat indexing.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct XMajorAdapter {
    lx: usize,
    ly: usize,
}

impl XMajorAdapter {
    /// Construct an adapter for an `Lx` by `Ly` rectangular array.
    pub fn new(lx: usize, ly: usize) -> Result<Self, GeometryError> {
        if lx == 0 {
            return Err(GeometryError::InvalidExtent {
                axis: "x",
                value: lx,
            });
        }
        if ly == 0 {
            return Err(GeometryError::InvalidExtent {
                axis: "y",
                value: ly,
            });
        }
        let total = lx.checked_mul(ly).ok_or(GeometryError::DimensionOverflow {
            operation: "site count",
        })?;
        if total > u32::MAX as usize {
            return Err(GeometryError::DimensionOverflow {
                operation: "site identifier",
            });
        }
        Ok(Self { lx, ly })
    }

    /// Convert a legacy x-major identifier to canonical row-major order.
    pub fn to_canonical(&self, legacy: SiteId) -> Result<SiteId, GeometryError> {
        let index = legacy.get() as usize;
        let total = self.lx * self.ly;
        if index >= total {
            return Err(GeometryError::SiteOutOfRange {
                site: legacy.get(),
                site_count: total,
            });
        }
        let x = index / self.ly;
        let y = index % self.ly;
        SiteId::try_from_usize(x + self.lx * y).map_err(|_| GeometryError::DimensionOverflow {
            operation: "site identifier",
        })
    }

    /// Convert a canonical row-major identifier to legacy x-major order.
    pub fn from_canonical(&self, canonical: SiteId) -> Result<SiteId, GeometryError> {
        let index = canonical.get() as usize;
        let total = self.lx * self.ly;
        if index >= total {
            return Err(GeometryError::SiteOutOfRange {
                site: canonical.get(),
                site_count: total,
            });
        }
        let x = index % self.lx;
        let y = index / self.lx;
        SiteId::try_from_usize(x * self.ly + y).map_err(|_| GeometryError::DimensionOverflow {
            operation: "site identifier",
        })
    }
}

/// A finite Cartesian coordinate associated with one canonical site.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coordinate {
    x: Real,
    y: Real,
}

impl Coordinate {
    /// Construct a finite Cartesian coordinate.
    pub fn new(x: Real, y: Real) -> Result<Self, GeometryError> {
        if x.is_finite() && y.is_finite() {
            Ok(Self { x, y })
        } else {
            Err(GeometryError::NonFiniteCoordinate { index: 0 })
        }
    }

    /// Return the horizontal coordinate.
    pub const fn x(self) -> Real {
        self.x
    }

    /// Return the vertical coordinate.
    pub const fn y(self) -> Real {
        self.y
    }

    /// Return `(x, y)` for compact comparisons and serialization adapters.
    pub const fn as_tuple(self) -> (Real, Real) {
        (self.x, self.y)
    }
}

/// An undirected bond with an optional periodic-image translation.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Bond {
    first: SiteId,
    second: SiteId,
    image_x: i32,
    image_y: i32,
    source: SiteId,
    direction_x: i8,
    direction_y: i8,
}

impl Bond {
    /// Construct a canonical simple bond between two distinct sites.
    pub fn new(a: SiteId, b: SiteId) -> Result<Self, GeometryError> {
        if a == b {
            return Err(GeometryError::SelfBond { site: a });
        }
        let (first, second) = if a < b { (a, b) } else { (b, a) };
        Ok(Self {
            first,
            second,
            image_x: 0,
            image_y: 0,
            source: first,
            direction_x: 0,
            direction_y: 0,
        })
    }

    /// Reconstruct a canonical bond from its complete serialized identity.
    ///
    /// This constructor is intended for durable artifact adapters. The
    /// endpoints must already be in canonical order; no site-order conversion
    /// or periodic-image reinterpretation is performed implicitly.
    pub fn from_parts(
        first: SiteId,
        second: SiteId,
        image_x: i32,
        image_y: i32,
        source: SiteId,
        direction_x: i8,
        direction_y: i8,
    ) -> Result<Self, GeometryError> {
        if first >= second || (source != first && source != second) {
            return Err(GeometryError::NonCanonicalBond { first, second });
        }
        Ok(Self {
            first,
            second,
            image_x,
            image_y,
            source,
            direction_x,
            direction_y,
        })
    }

    /// Return a canonical bond carrying a periodic image translation.
    pub(crate) fn with_image(
        a: SiteId,
        b: SiteId,
        image_x: i32,
        image_y: i32,
        direction_x: isize,
        direction_y: isize,
    ) -> Result<Self, GeometryError> {
        let mut bond = Self::new(a, b)?;
        if a > b {
            bond.image_x = -image_x;
            bond.image_y = -image_y;
        } else {
            bond.image_x = image_x;
            bond.image_y = image_y;
        }
        bond.source = a;
        bond.direction_x =
            i8::try_from(direction_x).map_err(|_| GeometryError::DimensionOverflow {
                operation: "bond direction",
            })?;
        bond.direction_y =
            i8::try_from(direction_y).map_err(|_| GeometryError::DimensionOverflow {
                operation: "bond direction",
            })?;
        Ok(bond)
    }

    /// Return the first canonical endpoint.
    pub const fn first(self) -> SiteId {
        self.first
    }

    /// Return the second canonical endpoint.
    pub const fn second(self) -> SiteId {
        self.second
    }

    /// Return the periodic image translation in lattice-cell units.
    pub const fn image_translation(self) -> (i32, i32) {
        (self.image_x, self.image_y)
    }

    /// Return the directed source site used to generate this image.
    pub const fn source(self) -> SiteId {
        self.source
    }

    /// Return the lattice step that generated this image.
    pub const fn direction(self) -> (i8, i8) {
        (self.direction_x, self.direction_y)
    }
}

impl Ord for Bond {
    fn cmp(&self, other: &Self) -> Ordering {
        (
            self.first,
            self.second,
            self.image_x,
            self.image_y,
            self.source,
            self.direction_x,
            self.direction_y,
        )
            .cmp(&(
                other.first,
                other.second,
                other.image_x,
                other.image_y,
                other.source,
                other.direction_x,
                other.direction_y,
            ))
    }
}

impl PartialOrd for Bond {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A rectangular lattice with explicit row-major site order and boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct RectangularGeometry {
    lx: usize,
    ly: usize,
    boundary_x: Boundary,
    boundary_y: Boundary,
    kind: LatticeKind,
    site_count: SiteCount,
}

impl RectangularGeometry {
    /// Construct an orthogonal rectangular geometry.
    pub fn new(
        lx: usize,
        ly: usize,
        boundary_x: Boundary,
        boundary_y: Boundary,
    ) -> Result<Self, GeometryError> {
        Self::with_kind(lx, ly, boundary_x, boundary_y, LatticeKind::Square)
    }

    /// Construct a rectangular coordinate grid with an explicit embedding.
    pub fn with_kind(
        lx: usize,
        ly: usize,
        boundary_x: Boundary,
        boundary_y: Boundary,
        kind: LatticeKind,
    ) -> Result<Self, GeometryError> {
        if lx == 0 {
            return Err(GeometryError::InvalidExtent {
                axis: "x",
                value: lx,
            });
        }
        if ly == 0 {
            return Err(GeometryError::InvalidExtent {
                axis: "y",
                value: ly,
            });
        }
        let site_total = lx.checked_mul(ly).ok_or(GeometryError::DimensionOverflow {
            operation: "site count",
        })?;
        if site_total > u32::MAX as usize {
            return Err(GeometryError::DimensionOverflow {
                operation: "site identifier",
            });
        }
        let site_count = SiteCount::new(site_total).map_err(|_| GeometryError::InvalidExtent {
            axis: "sites",
            value: site_total,
        })?;
        Ok(Self {
            lx,
            ly,
            boundary_x,
            boundary_y,
            kind,
            site_count,
        })
    }

    /// Number of sites.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Return the horizontal extent.
    pub const fn lx(&self) -> usize {
        self.lx
    }

    /// Return the vertical extent.
    pub const fn ly(&self) -> usize {
        self.ly
    }

    /// Return the horizontal boundary condition.
    pub const fn boundary_x(&self) -> Boundary {
        self.boundary_x
    }

    /// Return the vertical boundary condition.
    pub const fn boundary_y(&self) -> Boundary {
        self.boundary_y
    }

    /// Return the selected embedding.
    pub const fn kind(&self) -> LatticeKind {
        self.kind
    }

    /// Map row-major `(x, y)` coordinates to a site identifier.
    pub fn site_id(&self, x: usize, y: usize) -> Result<SiteId, GeometryError> {
        if x >= self.lx || y >= self.ly {
            return Err(GeometryError::CoordinateOutOfRange {
                x: x as isize,
                y: y as isize,
            });
        }
        let index = x
            .checked_add(
                y.checked_mul(self.lx)
                    .ok_or(GeometryError::DimensionOverflow {
                        operation: "row-major index",
                    })?,
            )
            .ok_or(GeometryError::DimensionOverflow {
                operation: "row-major index",
            })?;
        SiteId::try_from_usize(index).map_err(|error| match error {
            BasisError::IdentifierOverflow { value } => GeometryError::DimensionOverflow {
                operation: if value > u32::MAX as usize {
                    "site identifier"
                } else {
                    "row-major index"
                },
            },
            _ => GeometryError::DimensionOverflow {
                operation: "site identifier",
            },
        })
    }

    /// Map a canonical site identifier to its row-major coordinate.
    pub fn coordinate(&self, site: SiteId) -> Result<Coordinate, GeometryError> {
        if site.get() as usize >= self.lx * self.ly {
            return Err(GeometryError::SiteOutOfRange {
                site: site.get(),
                site_count: self.lx * self.ly,
            });
        }
        let x = (site.get() as usize) % self.lx;
        let y = (site.get() as usize) / self.lx;
        let (x, y) = (x as Real, y as Real);
        match self.kind {
            LatticeKind::Square => Coordinate::new(x, y),
            LatticeKind::Triangular => Coordinate::new(x + 0.5 * y, 3.0_f64.sqrt() * 0.5 * y),
        }
    }

    /// Return the canonical minimum-image Cartesian displacement `r_b - r_a`.
    pub fn minimum_image_displacement(
        &self,
        a: SiteId,
        b: SiteId,
    ) -> Result<Coordinate, GeometryError> {
        let ca = self.coordinate(a)?;
        let cb = self.coordinate(b)?;
        let raw = (cb.x - ca.x, cb.y - ca.y);
        if self.kind == LatticeKind::Square {
            let dx = if self.boundary_x == Boundary::Periodic {
                raw.0
                    - self.lx as Real * ((raw.0 + 0.5 * self.lx as Real) / self.lx as Real).floor()
            } else {
                raw.0
            };
            let dy = if self.boundary_y == Boundary::Periodic {
                raw.1
                    - self.ly as Real * ((raw.1 + 0.5 * self.ly as Real) / self.ly as Real).floor()
            } else {
                raw.1
            };
            return Coordinate::new(dx, dy);
        }
        let mut best = raw;
        let mut best_norm = raw.0 * raw.0 + raw.1 * raw.1;
        let mut best_translation = (0_i64, 0_i64);
        let cell_y = 3.0_f64.sqrt() * 0.5 * self.ly as Real;
        let radius = best_norm.sqrt();
        let (y_start, y_end) = if self.boundary_y == Boundary::Periodic {
            (
                ((raw.1 - radius) / cell_y).floor() as i64,
                ((raw.1 + radius) / cell_y).ceil() as i64,
            )
        } else {
            (0, 0)
        };
        for ty in y_start..=y_end {
            let residual_x = raw.0 - 0.5 * self.ly as Real * ty as Real;
            let x_radius = (best_norm - (raw.1 - cell_y * ty as Real).powi(2))
                .max(0.0)
                .sqrt();
            let (x_start, x_end) = if self.boundary_x == Boundary::Periodic {
                (
                    ((residual_x - x_radius) / self.lx as Real).floor() as i64,
                    ((residual_x + x_radius) / self.lx as Real).ceil() as i64,
                )
            } else {
                (0, 0)
            };
            for tx in x_start..=x_end {
                if tx == 0 && ty == 0 {
                    continue;
                }
                let sx = tx as Real * self.lx as Real + 0.5 * ty as Real * self.ly as Real;
                let sy = cell_y * ty as Real;
                let candidate = (raw.0 - sx, raw.1 - sy);
                let norm = candidate.0 * candidate.0 + candidate.1 * candidate.1;
                if norm < best_norm || (norm == best_norm && (tx, ty) < best_translation) {
                    best = candidate;
                    best_norm = norm;
                    best_translation = (tx, ty);
                }
            }
        }
        Coordinate::new(best.0, best.1)
    }

    /// Return unordered pairs whose minimum-image squared distance lies in a shell.
    pub fn pairs_at_squared_distance(
        &self,
        target_squared: Real,
        tolerance: ShellTolerance,
    ) -> Result<Vec<Bond>, GeometryError> {
        if !target_squared.is_finite() || target_squared < 0.0 {
            return Err(GeometryError::InvalidShell {
                value: target_squared,
            });
        }
        let epsilon = match tolerance {
            ShellTolerance::Absolute(value) | ShellTolerance::Relative(value)
                if value.is_finite() && value >= 0.0 =>
            {
                value
            }
            ShellTolerance::Absolute(value) | ShellTolerance::Relative(value) => {
                return Err(GeometryError::InvalidShell { value });
            }
        };
        let mut result = Vec::new();
        for i in 0..self.site_count.get() {
            for j in (i + 1)..self.site_count.get() {
                let displacement =
                    self.minimum_image_displacement(SiteId::new(i as u32), SiteId::new(j as u32))?;
                let squared = displacement.x * displacement.x + displacement.y * displacement.y;
                let allowed = match tolerance {
                    ShellTolerance::Absolute(_) => epsilon,
                    ShellTolerance::Relative(_) => epsilon * target_squared.abs(),
                };
                if (squared - target_squared).abs() <= allowed {
                    result.push(Bond::new(SiteId::new(i as u32), SiteId::new(j as u32))?);
                }
            }
        }
        Ok(result)
    }

    /// Generate canonical bonds under the requested multiplicity policy.
    pub fn bonds(&self, multiplicity: BondMultiplicity) -> Result<Vec<Bond>, GeometryError> {
        let mut bonds = Vec::new();
        let directions: &[(isize, isize)] = match self.kind {
            LatticeKind::Square => &[(1, 0), (0, 1)],
            LatticeKind::Triangular => &[(1, 0), (0, 1), (1, -1)],
        };
        for y in 0..self.ly {
            for x in 0..self.lx {
                let source = self.site_id(x, y)?;
                for &(dx, dy) in directions {
                    if let Some((target, ix, iy)) = self.step(x, y, dx, dy) {
                        let candidate = if multiplicity == BondMultiplicity::Simple {
                            Bond::new(source, target)?
                        } else {
                            Bond::with_image(source, target, ix, iy, dx, dy)?
                        };
                        if multiplicity == BondMultiplicity::PeriodicImages
                            || !bonds.iter().any(|bond: &Bond| {
                                bond.first() == candidate.first()
                                    && bond.second() == candidate.second()
                            })
                        {
                            bonds.push(candidate);
                        }
                    }
                }
            }
        }
        bonds.sort_unstable();
        Ok(bonds)
    }

    fn step(&self, x: usize, y: usize, dx: isize, dy: isize) -> Option<(SiteId, i32, i32)> {
        let nx = x as isize + dx;
        let ny = y as isize + dy;
        let (nx, ix) = wrap(nx, self.lx, self.boundary_x)?;
        let (ny, iy) = wrap(ny, self.ly, self.boundary_y)?;
        let target = self.site_id(nx, ny).ok()?;
        if target == self.site_id(x, y).ok()? {
            return None;
        }
        Some((target, ix, iy))
    }
}

fn wrap(value: isize, extent: usize, boundary: Boundary) -> Option<(usize, i32)> {
    if (0..extent as isize).contains(&value) {
        return Some((value as usize, 0));
    }
    match boundary {
        Boundary::Open => None,
        Boundary::Periodic => {
            let extent = extent as isize;
            let quotient = value.div_euclid(extent);
            Some((value.rem_euclid(extent) as usize, quotient as i32))
        }
    }
}

/// A custom ordered collection of finite Cartesian sites.
#[derive(Clone, Debug, PartialEq)]
pub struct CustomGeometry {
    coordinates: Vec<Coordinate>,
    site_count: SiteCount,
}

impl CustomGeometry {
    /// Construct a custom geometry; vector position is the canonical site id.
    pub fn new(coordinates: Vec<Coordinate>) -> Result<Self, GeometryError> {
        if coordinates.is_empty() {
            return Err(GeometryError::EmptyCustomGeometry);
        }
        if coordinates.len() > u32::MAX as usize {
            return Err(GeometryError::DimensionOverflow {
                operation: "site identifier",
            });
        }
        let site_count =
            SiteCount::new(coordinates.len()).map_err(|_| GeometryError::EmptyCustomGeometry)?;
        Ok(Self {
            coordinates,
            site_count,
        })
    }

    /// Number of custom sites.
    pub const fn site_count(&self) -> SiteCount {
        self.site_count
    }

    /// Return a custom site's coordinate.
    pub fn coordinate(&self, site: SiteId) -> Result<Coordinate, GeometryError> {
        self.coordinates
            .get(site.get() as usize)
            .copied()
            .ok_or(GeometryError::SiteOutOfRange {
                site: site.get(),
                site_count: self.coordinates.len(),
            })
    }

    /// Generate the complete simple graph of unordered custom-site pairs.
    pub fn complete_bonds(&self) -> Result<Vec<Bond>, GeometryError> {
        let mut result = Vec::new();
        for i in 0..self.coordinates.len() {
            for j in (i + 1)..self.coordinates.len() {
                result.push(Bond::new(
                    SiteId::try_from_usize(i).map_err(|_| GeometryError::DimensionOverflow {
                        operation: "site identifier",
                    })?,
                    SiteId::try_from_usize(j).map_err(|_| GeometryError::DimensionOverflow {
                        operation: "site identifier",
                    })?,
                )?);
            }
        }
        Ok(result)
    }
}
