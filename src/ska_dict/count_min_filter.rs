//! A Countmin filter for counting k-mers from reads
//!
//! This is to support a minimum count from FASTQ files, for basic error
//! correction, while keeping memory use low and predictable. Use the
//! [`CountMin`] struct.
//!
//! See <https://en.wikipedia.org/wiki/Count-min_sketch> for more
//! details on this data structure.

use ahash::RandomState;

/// A Countmin table of specified width and height, counts input k-mers, returns
/// whether they have passed a count threshold.
///
/// Table can be reset to avoid reallocation for every sample.
///
/// # Examples
///
/// ```
/// use ska::ska_dict::count_min_filter::CountMin;
///
/// // Create the filter
/// let width: usize = 1 << 26;
/// let height = 6;
/// let min_count = 2;
/// let mut cm_filter = CountMin::empty(width, height, min_count);
/// cm_filter.init();
///
/// // Use it
/// let mut passed = cm_filter.filter(0 as u64); // false
/// passed = cm_filter.filter(0 as u64);         // true
///
/// // Reset for use with a new file
/// cm_filter.reset();
/// ```
#[derive(Debug, Clone)]
pub struct CountMin {
    /// Table width (estimated number of unique k-mers)
    width: usize,
    /// Table height (number of hashes)
    height: usize,
    /// Hash generators
    hash_factory: Vec<RandomState>,
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
        let mask = width as u64 - 1;

        // Reserve for these gets call by the vec! macro used in init
        let hash_factory = Vec::new();
        let counts = Vec::new();

        Self {
            width,
            height,
            hash_factory,
            mask,
            counts,
            min_count,
        }
    }

    /// Initialises table so it is ready for use.
    ///
    /// Allocates memory and sets hash functions.
    pub fn init(&mut self) {
        if self.counts.is_empty() {
            self.hash_factory = vec![RandomState::new(); self.height];
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

    /// Add an observation of a k-mer to the filter, and return if it passed
    /// minimum count filtering criterion.
    pub fn filter(&mut self, kmer: u64) -> bool {
        let mut count = 0;
        for hash_it in self.hash_factory.iter().enumerate() {
            let (row_idx, hash) = hash_it;
            let table_idx = row_idx * self.width + ((hash.hash_one(kmer) & self.mask) as usize);
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
