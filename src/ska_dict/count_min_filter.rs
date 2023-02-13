//! A Countmin filter for counting k-mers from reads
//!
//! This is to support a minimum count from FASTQ files, for basic error
//! correction, while keeping memory use low and predictable. Use the
//! [`CountMin`] struct.
//!
//! See <https://en.wikipedia.org/wiki/Count-min_sketch> for more
//! details on this data structure.

use super::bit_encoding::UInt;
use super::split_kmer::SplitKmer;

/// A Countmin table of specified width and height, counts input k-mers, returns
/// whether they have passed a count threshold.
///
/// Table can be reset to avoid reallocation for every sample.
#[derive(Debug, Clone)]
pub struct CountMin {
    /// Table width (estimated number of unique k-mers)
    width: usize,
    /// Number of bits to shift masked value to get table entry
    width_shift: u32,
    /// Table height (number of hashes)
    height: usize,
    /// Mask to convert hash into table column
    mask: u64,
    /// Table of counts
    counts: Vec<u16>,
    /// Minimum count to pass filter
    min_count: u16,
}

impl CountMin {
    /// Creates a new countmin table of specified size, and pass threshold
    ///
    /// Note:
    /// - Must call [`CountMin::init()`] before using.
    /// - Width will be rounded down to the closest power of two
    /// - Maximum count is [`u16::MAX`] i.e. 65535
    pub fn empty(width: usize, height: usize, min_count: u16) -> Self {
        // Consistent with consts used, but ensures a power of two
        let width_bits: usize = f64::floor(f64::log2(width as f64)) as usize;
        let width = 1 << (width_bits + 1);
        // Use MSB rather than LSB
        let width_shift = u64::BITS - width_bits as u32;
        let mask = (width as u64 - 1) << width_shift;

        Self {
            width,
            width_shift,
            height,
            mask,
            counts: Vec::new(),
            min_count,
        }
    }

    /// Initialises table so it is ready for use.
    ///
    /// Allocates memory and sets hash functions.
    pub fn init(&mut self) {
        if self.counts.is_empty() {
            self.counts = vec![0; self.width * self.height];
        }
    }

    /// Check if [`CountMin::init()`] has been called and table is ready for use
    pub fn is_init(&self) -> bool {
        !self.counts.is_empty()
    }

    /// Reset all counts to zero
    pub fn reset(&mut self) {
        self.counts = vec![0; self.width * self.height];
    }

    /// Add an observation of a k-mer and middle base to the filter, and return if it passed
    /// minimum count filtering criterion.
    pub fn filter<IntT: for<'a> UInt<'a>>(&mut self, kmer: &SplitKmer<IntT>) -> bool {
        // This is possible because of the k-mer size restriction, the top two
        // bit are always zero
        let mut count = 0;
        for row_idx in 0..self.height {
            let hash_val = kmer.get_hash(row_idx);
            let table_idx =
                row_idx * self.width + (((hash_val & self.mask) >> self.width_shift) as usize);
            self.counts[table_idx] = self.counts[table_idx].saturating_add(1);
            if row_idx == 0 {
                count = self.counts[table_idx];
            } else {
                count = u16::min(count, self.counts[table_idx]);
            }
        }
        count >= self.min_count
    }
}
