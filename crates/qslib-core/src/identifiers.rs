use crate::BasisError;
use std::fmt::{self, Display, Formatter};
use std::str::FromStr;

/// A checked number of sites in a finite simulation system.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Hash)]
pub struct SiteCount(usize);

impl SiteCount {
    /// Construct a non-empty site count.
    pub fn new(value: usize) -> Result<Self, BasisError> {
        if value == 0 {
            Err(BasisError::EmptySystem)
        } else {
            Ok(Self(value))
        }
    }

    /// Return the number of sites.
    pub const fn get(self) -> usize {
        self.0
    }

    /// Validate a site identifier against this count.
    pub fn validate(self, site: SiteId) -> Result<(), BasisError> {
        if (site.0 as usize) < self.0 {
            Ok(())
        } else {
            Err(BasisError::SiteOutOfRange {
                site: site.0,
                site_count: self.0,
            })
        }
    }
}

/// A zero-based site identifier with a stable serialized representation.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Hash)]
pub struct SiteId(u32);

impl SiteId {
    /// Construct an identifier from its already checked `u32` value.
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// Convert a platform-sized index without truncation.
    pub fn try_from_usize(value: usize) -> Result<Self, BasisError> {
        u32::try_from(value)
            .map(Self)
            .map_err(|_| BasisError::IdentifierOverflow { value })
    }

    /// Return the numeric identifier.
    pub const fn get(self) -> u32 {
        self.0
    }
}

/// A physical Pauli axis named independently from a stored simulation basis.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PhysicalAxis {
    /// Physical x axis.
    X,
    /// Physical y axis.
    Y,
    /// Physical z axis.
    Z,
}

impl PhysicalAxis {
    /// Return the canonical lowercase serialized spelling.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
        }
    }
}

impl Display for PhysicalAxis {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl FromStr for PhysicalAxis {
    type Err = BasisError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "x" => Ok(Self::X),
            "y" => Ok(Self::Y),
            "z" => Ok(Self::Z),
            value => Err(BasisError::InvalidAxis {
                value: value.to_owned(),
            }),
        }
    }
}

/// The physical Pauli axis diagonal in the stored binary basis.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum SimulationBasis {
    /// The physical x operator is diagonal.
    X,
    /// The physical y operator is diagonal.
    Y,
    /// The physical z operator is diagonal.
    Z,
}

impl SimulationBasis {
    /// Return the canonical lowercase serialized spelling.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
        }
    }
}

impl Display for SimulationBasis {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl FromStr for SimulationBasis {
    type Err = BasisError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "x" => Ok(Self::X),
            "y" => Ok(Self::Y),
            "z" => Ok(Self::Z),
            value => Err(BasisError::InvalidAxis {
                value: value.to_owned(),
            }),
        }
    }
}
