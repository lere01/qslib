use std::fmt::{self, Display, Formatter};

/// Errors raised while validating identifiers, basis states, and packed data.
#[derive(Clone, Debug, PartialEq)]
pub enum BasisError {
    /// A finite lattice or state was requested with zero sites.
    EmptySystem,
    /// A NaN or infinity was supplied where a physical scalar is required.
    NonFiniteScalar {
        /// The supplied non-finite value.
        value: f64,
    },
    /// A checked size calculation would overflow the platform representation.
    DimensionOverflow {
        /// Name of the calculation that overflowed.
        operation: &'static str,
    },
    /// An axis or simulation basis used a non-canonical spelling.
    InvalidAxis {
        /// The supplied non-canonical spelling.
        value: String,
    },
    /// A platform-sized integer could not be represented by a `SiteId`.
    IdentifierOverflow {
        /// The platform-sized value that could not be represented.
        value: usize,
    },
    /// A site identifier does not belong to the declared system.
    SiteOutOfRange {
        /// The invalid site identifier.
        site: u32,
        /// The exclusive upper bound for valid sites.
        site_count: usize,
    },
    /// A raw dense state contained a value other than zero or one.
    InvalidBit {
        /// Position of the invalid value in the dense input.
        index: usize,
        /// The supplied raw value.
        value: u8,
    },
    /// The supplied number of machine words differs from the checked size.
    InvalidWordCount {
        /// Number of words required by the site count.
        expected: usize,
        /// Number of words supplied by the caller.
        actual: usize,
    },
    /// Bits above the declared site count are set in the final word.
    NonCanonicalHighBits {
        /// Index of the word containing padding bits.
        word_index: usize,
        /// Supplied word value.
        value: u64,
        /// Mask of bits that are valid for the final word.
        valid_mask: u64,
    },
    /// Serialized bytes do not contain exactly the declared number of words.
    SerializedLength {
        /// Number of bytes required by the selected word width.
        expected: usize,
        /// Number of bytes supplied by the caller.
        actual: usize,
    },
    /// A fixed-weight sector requests more occupied bits than sites.
    WeightOutOfRange {
        /// Requested number of one bits.
        weight: usize,
        /// Number of sites available.
        site_count: usize,
    },
}

impl Display for BasisError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySystem => formatter.write_str("a system must contain at least one site"),
            Self::DimensionOverflow { operation } => {
                write!(formatter, "checked {operation} dimension overflowed")
            }
            Self::InvalidAxis { value } => {
                write!(
                    formatter,
                    "axis or simulation basis {value:?} is not canonical"
                )
            }
            Self::NonFiniteScalar { value } => write!(formatter, "non-finite scalar {value:?}"),
            Self::IdentifierOverflow { value } => {
                write!(formatter, "site identifier {value} does not fit in u32")
            }
            Self::SiteOutOfRange { site, site_count } => {
                write!(
                    formatter,
                    "site {site} is outside a {site_count}-site system"
                )
            }
            Self::InvalidBit { index, value } => {
                write!(
                    formatter,
                    "bit {value} at position {index} is not zero or one"
                )
            }
            Self::InvalidWordCount { expected, actual } => write!(
                formatter,
                "packed state requires {expected} words, received {actual}"
            ),
            Self::NonCanonicalHighBits {
                word_index,
                value,
                valid_mask,
            } => write!(
                formatter,
                "packed word {word_index} has non-canonical high bits: {value:#x} (valid mask {valid_mask:#x})"
            ),
            Self::SerializedLength { expected, actual } => write!(
                formatter,
                "serialized packed state requires {expected} bytes, received {actual}"
            ),
            Self::WeightOutOfRange { weight, site_count } => write!(
                formatter,
                "Hamming weight {weight} exceeds the {site_count}-site system"
            ),
        }
    }
}

impl std::error::Error for BasisError {}
