use crate::{BasisError, SiteCount};

/// A binary value in canonical site order. It is not a model-specific spin
/// label or a Rydberg occupation until an explicit model interprets it.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum BasisBit {
    /// Bit zero, with eigenvalue `+1` for the diagonal simulation axis.
    Zero,
    /// Bit one, with eigenvalue `-1` for the diagonal simulation axis.
    One,
}

impl BasisBit {
    /// Return the packed numeric value.
    pub const fn as_u8(self) -> u8 {
        match self {
            Self::Zero => 0,
            Self::One => 1,
        }
    }

    /// Return the Pauli eigenvalue associated with this bit in the simulation
    /// basis.
    pub const fn pauli_eigenvalue(self) -> i8 {
        match self {
            Self::Zero => 1,
            Self::One => -1,
        }
    }
}

impl TryFrom<u8> for BasisBit {
    type Error = BasisError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Zero),
            1 => Ok(Self::One),
            value => Err(BasisError::InvalidBit { index: 0, value }),
        }
    }
}

/// An owned dense state in canonical site order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BasisState {
    bits: Vec<BasisBit>,
}

impl BasisState {
    /// Construct a state from raw zero/one values with checked validation.
    pub fn from_raw_bits(values: &[u8]) -> Result<Self, BasisError> {
        if values.is_empty() {
            return Err(BasisError::EmptySystem);
        }
        let bits = values
            .iter()
            .copied()
            .enumerate()
            .map(|(index, value)| match value {
                0 => Ok(BasisBit::Zero),
                1 => Ok(BasisBit::One),
                value => Err(BasisError::InvalidBit { index, value }),
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { bits })
    }

    /// Construct a state from already validated binary values.
    pub fn from_bits(bits: &[BasisBit]) -> Result<Self, BasisError> {
        if bits.is_empty() {
            Err(BasisError::EmptySystem)
        } else {
            Ok(Self {
                bits: bits.to_vec(),
            })
        }
    }

    /// Return the number of sites.
    pub fn len(&self) -> usize {
        self.bits.len()
    }

    /// Return whether this state has no sites. Valid states are never empty;
    /// this method makes generic collection code explicit.
    pub fn is_empty(&self) -> bool {
        self.bits.is_empty()
    }

    /// Borrow bits in canonical site order.
    pub fn bits(&self) -> &[BasisBit] {
        &self.bits
    }

    /// Borrow this state without copying its dense bit storage.
    pub fn as_view(&self) -> BasisStateView<'_> {
        BasisStateView { bits: &self.bits }
    }

    /// Return the Hamming weight.
    pub fn hamming_weight(&self) -> usize {
        self.bits
            .iter()
            .filter(|bit| **bit == BasisBit::One)
            .count()
    }

    /// Return Pauli eigenvalues for the diagonal simulation axis.
    pub fn pauli_eigenvalues(&self) -> Vec<i8> {
        self.bits.iter().map(|bit| bit.pauli_eigenvalue()).collect()
    }

    /// Pack this state little-endian by site identifier.
    pub fn pack(&self) -> Result<PackedState, BasisError> {
        PackedState::from_bits(&self.bits)
    }
}

/// An immutable dense view over a canonical binary state.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BasisStateView<'a> {
    bits: &'a [BasisBit],
}

impl<'a> BasisStateView<'a> {
    /// Return the number of sites in the view.
    pub const fn len(self) -> usize {
        self.bits.len()
    }

    /// Return whether the view is empty. Valid qslib states are non-empty.
    pub const fn is_empty(self) -> bool {
        self.bits.is_empty()
    }

    /// Borrow bits in canonical site order.
    pub const fn bits(self) -> &'a [BasisBit] {
        self.bits
    }

    /// Return the Hamming weight without allocating.
    pub fn hamming_weight(self) -> usize {
        self.bits
            .iter()
            .filter(|bit| **bit == BasisBit::One)
            .count()
    }
}

/// Width of each word in a serialized packed-state byte stream.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum WordWidth {
    /// Eight bits per serialized word.
    U8,
    /// Sixteen bits per serialized word.
    U16,
    /// Thirty-two bits per serialized word.
    U32,
    /// Sixty-four bits per serialized word.
    U64,
}

impl WordWidth {
    /// Return the number of bits in one serialized word.
    pub const fn bits(self) -> usize {
        match self {
            Self::U8 => 8,
            Self::U16 => 16,
            Self::U32 => 32,
            Self::U64 => 64,
        }
    }

    /// Return the number of bytes in one serialized word.
    pub const fn bytes(self) -> usize {
        self.bits() / 8
    }
}

/// A little-endian-by-site packed binary state.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PackedState {
    site_count: SiteCount,
    words: Vec<u64>,
}

impl PackedState {
    /// Construct a packed state from dense bits.
    pub fn from_bits(bits: &[BasisBit]) -> Result<Self, BasisError> {
        let site_count = SiteCount::new(bits.len())?;
        let mut words = vec![0; word_count(site_count.get())];
        for (index, bit) in bits.iter().copied().enumerate() {
            if bit == BasisBit::One {
                words[index / 64] |= 1_u64 << (index % 64);
            }
        }
        Ok(Self { site_count, words })
    }

    /// Construct from little-endian machine words and validate unused bits.
    pub fn from_words(site_count: usize, words: &[u64]) -> Result<Self, BasisError> {
        let site_count = SiteCount::new(site_count)?;
        let expected = word_count(site_count.get());
        if words.len() != expected {
            return Err(BasisError::InvalidWordCount {
                expected,
                actual: words.len(),
            });
        }
        let mut owned = words.to_vec();
        validate_high_bits(site_count.get(), &owned)?;
        if let Some(last) = owned.last_mut() {
            *last &= valid_mask(site_count.get());
        }
        Ok(Self {
            site_count,
            words: owned,
        })
    }

    /// Construct from width-labelled little-endian serialized words.
    pub fn from_bytes(
        site_count: usize,
        width: WordWidth,
        bytes: &[u8],
    ) -> Result<Self, BasisError> {
        let sites = SiteCount::new(site_count)?;
        let serialized_words = word_count_for_width(sites.get(), width);
        let expected_bytes =
            serialized_words
                .checked_mul(width.bytes())
                .ok_or(BasisError::DimensionOverflow {
                    operation: "serialized byte length",
                })?;
        if bytes.len() != expected_bytes {
            return Err(BasisError::SerializedLength {
                expected: expected_bytes,
                actual: bytes.len(),
            });
        }
        let mut words = vec![0; word_count(sites.get())];
        for serialized_index in 0..serialized_words {
            let start = serialized_index * width.bytes();
            let mut value = 0_u64;
            for (offset, byte) in bytes[start..start + width.bytes()]
                .iter()
                .copied()
                .enumerate()
            {
                value |= u64::from(byte) << (offset * 8);
            }
            let base = serialized_index.checked_mul(width.bits()).ok_or(
                BasisError::DimensionOverflow {
                    operation: "serialized bit offset",
                },
            )?;
            for bit in 0..width.bits() {
                let position = base.checked_add(bit).ok_or(BasisError::DimensionOverflow {
                    operation: "serialized bit position",
                })?;
                if position >= sites.get() {
                    if value & (1_u64 << bit) != 0 {
                        return Err(BasisError::NonCanonicalHighBits {
                            word_index: serialized_index,
                            value,
                            valid_mask: valid_mask_for_width(sites.get() - base, width),
                        });
                    }
                } else if value & (1_u64 << bit) != 0 {
                    words[position / 64] |= 1_u64 << (position % 64);
                }
            }
        }
        Self::from_words(sites.get(), &words)
    }

    /// Return the number of sites represented by this state.
    pub const fn site_count(&self) -> usize {
        self.site_count.get()
    }

    /// Borrow internal words ordered from least to most significant.
    pub fn words_le(&self) -> &[u64] {
        &self.words
    }

    /// Read one site by its zero-based index.
    pub fn bit(&self, site: usize) -> Result<BasisBit, BasisError> {
        if site >= self.site_count() {
            return if let Ok(site) = u32::try_from(site) {
                Err(BasisError::SiteOutOfRange {
                    site,
                    site_count: self.site_count(),
                })
            } else {
                Err(BasisError::IdentifierOverflow { value: site })
            };
        }
        let value = (self.words[site / 64] >> (site % 64)) & 1;
        Ok(if value == 0 {
            BasisBit::Zero
        } else {
            BasisBit::One
        })
    }

    /// Return the Hamming weight.
    pub fn hamming_weight(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }

    /// Serialize into fixed-width little-endian words.
    pub fn to_bytes(&self, width: WordWidth) -> Result<Vec<u8>, BasisError> {
        let words = word_count_for_width(self.site_count(), width);
        let byte_count = words
            .checked_mul(width.bytes())
            .ok_or(BasisError::DimensionOverflow {
                operation: "serialized byte length",
            })?;
        let mut bytes = vec![0; byte_count];
        for word_index in 0..words {
            let base =
                word_index
                    .checked_mul(width.bits())
                    .ok_or(BasisError::DimensionOverflow {
                        operation: "serialized bit offset",
                    })?;
            let mut value = 0_u64;
            for bit in 0..width.bits() {
                let position = base.checked_add(bit).ok_or(BasisError::DimensionOverflow {
                    operation: "serialized bit position",
                })?;
                if position < self.site_count() && self.bit(position)? == BasisBit::One {
                    value |= 1_u64 << bit;
                }
            }
            let start = word_index * width.bytes();
            let encoded = value.to_le_bytes();
            bytes[start..start + width.bytes()].copy_from_slice(&encoded[..width.bytes()]);
        }
        Ok(bytes)
    }

    fn increment(&mut self) -> bool {
        let mut carry = true;
        for word in &mut self.words {
            if carry {
                let (next, overflow) = word.overflowing_add(1);
                *word = next;
                carry = overflow;
            }
        }
        if carry
            || self
                .words
                .last()
                .is_some_and(|word| *word & !valid_mask(self.site_count()) != 0)
        {
            self.words.fill(0);
            true
        } else {
            false
        }
    }
}

/// Canonical full computational basis in increasing packed-integer order.
pub struct FullBasis {
    current: PackedState,
    done: bool,
}

impl FullBasis {
    /// Construct a full basis iterator for a non-empty system.
    pub fn new(site_count: SiteCount) -> Result<Self, BasisError> {
        if site_count.get() >= usize::BITS as usize {
            return Err(BasisError::DimensionOverflow {
                operation: "full basis dimension",
            });
        }
        Ok(Self {
            current: PackedState {
                site_count,
                words: vec![0; word_count(site_count.get())],
            },
            done: false,
        })
    }
}

impl Iterator for FullBasis {
    type Item = PackedState;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            None
        } else {
            let output = self.current.clone();
            self.done = self.current.increment();
            Some(output)
        }
    }
}

/// Canonical fixed-Hamming-weight basis in increasing packed-integer order.
pub struct SectorBasis {
    site_count: SiteCount,
    positions: Vec<usize>,
    done: bool,
}

impl SectorBasis {
    /// Construct the basis of states with exactly `weight` one bits.
    pub fn new(site_count: SiteCount, weight: usize) -> Result<Self, BasisError> {
        if weight > site_count.get() {
            return Err(BasisError::WeightOutOfRange {
                weight,
                site_count: site_count.get(),
            });
        }
        Ok(Self {
            site_count,
            positions: (0..weight).collect(),
            done: false,
        })
    }
}

impl Iterator for SectorBasis {
    type Item = PackedState;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let mut words = vec![0; word_count(self.site_count.get())];
        for position in &self.positions {
            words[*position / 64] |= 1_u64 << (*position % 64);
        }
        let output = PackedState {
            site_count: self.site_count,
            words,
        };
        let mut advanced = false;
        for index in 0..self.positions.len() {
            let upper = self
                .positions
                .get(index + 1)
                .copied()
                .unwrap_or(self.site_count.get());
            if self.positions[index] + 1 < upper {
                self.positions[index] += 1;
                for lower in 0..index {
                    self.positions[lower] = lower;
                }
                advanced = true;
                break;
            }
        }
        if !advanced {
            self.done = true;
        }
        Some(output)
    }
}

fn word_count(site_count: usize) -> usize {
    site_count.div_ceil(64)
}

fn word_count_for_width(site_count: usize, width: WordWidth) -> usize {
    site_count.div_ceil(width.bits())
}

fn valid_mask(site_count: usize) -> u64 {
    let remainder = site_count % 64;
    if remainder == 0 {
        u64::MAX
    } else {
        (1_u64 << remainder) - 1
    }
}

fn valid_mask_for_width(valid_bits: usize, width: WordWidth) -> u64 {
    if valid_bits >= width.bits() {
        u64::MAX
    } else {
        (1_u64 << valid_bits) - 1
    }
}

fn validate_high_bits(site_count: usize, words: &[u64]) -> Result<(), BasisError> {
    let mask = valid_mask(site_count);
    if let Some((word_index, value)) = words
        .iter()
        .copied()
        .enumerate()
        .next_back()
        .filter(|(_, value)| *value & !mask != 0)
    {
        return Err(BasisError::NonCanonicalHighBits {
            word_index,
            value,
            valid_mask: mask,
        });
    }
    Ok(())
}
