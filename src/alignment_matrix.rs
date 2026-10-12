//! Sparse, sample-ID/base packed alignment storage for the WASM workflow.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::hash::Hash;
use std::io::{self, Write};

use hashbrown::hash_map::Entry;
use hashbrown::HashMap;
use num_traits::{AsPrimitive, PrimInt};

const GAP: u8 = 15; // The symbol for gaps in alignments, "-". Currently unused
const SAMPLE_ID_BITS: usize = 11; // Sample ID bits, these determine the maximum amount of samples that can be added
pub(crate) const MAX_SAMPLES: usize = 1 << SAMPLE_ID_BITS;

// These are masks and bit positions for the data structures while processing. The basic structure is the following,
// with squared brackets ([]) delimiting u16s.
//
// First u16:  [4 bits for sample, 1 bit for extending or not, 11 bits for the total number of samples inside]
// Second u16: [10 bits for the number of ranges, 4 bits for the padding information, 2 bits that are used to store the actual data]
// Note the actual data is, in general, first the ranges (pairs of sample ids as 11bits), then the individual IDs (just 11 bits in sequence). However, if no ranges are present, we directly start with individual bits.

const SAMPLE_ID_MASK: u16 = (MAX_SAMPLES as u16) - 1; // This is a mask of the lower-value bits of the sample id ones.
const BASE_MASK: u16 = 0x000f; // this is a mask of four bits (1111).
const EXTEND_MASK: u16 = 0x0010; // This is a mask of just a single bit located in the fifth position, 10000. This determines whether we should try or not to extend the last range there is when adding new samples.
const RANGE_COUNT_MASK: u16 = 0x03ff; // This is a mask of the 10 lower bits (from 0 to 9): 1111111111.
const PADDING_MASK: u16 = 0x3c00; // This extracts four bits from the middle of a u16: 11110000000000.
const PAYLOAD_BIT: usize = 30; // this is the bit where data starts: the 30th.

// Compactification info
const MIN_RANGE_SAMPLES: usize = 3; // Min. number of consecutive samples we should have to create a range from them
const MAX_GROUP_WORDS: usize = (PAYLOAD_BIT + SAMPLE_ID_BITS * MAX_SAMPLES).div_ceil(16); // Total number of "words" (u16s in a vector) we can have ever.

/// Concrete calls share the same codes in extracted and imported alignments.
#[inline]
pub(crate) fn concrete_base_code(symbol: u8) -> Option<u8> {
    match symbol.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

#[inline]
fn is_concrete_code(code: u8) -> bool {
    code < 4
}

/// Encode an ASCII IUPAC symbol as a four-bit value.
#[inline]
pub(crate) fn encode_symbol(symbol: u8) -> u8 {
    match symbol.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        b'R' => 4,
        b'Y' => 5,
        b'S' => 6,
        b'W' => 7,
        b'K' => 8,
        b'M' => 9,
        b'B' => 10,
        b'D' => 11,
        b'H' => 12,
        b'V' => 13,
        b'N' => 14,
        b'-' => GAP,
        _ => panic!("Unsupported alignment symbol: {}", symbol as char),
    }
}

/// Decode a four-bit value to its ASCII IUPAC symbol.
#[inline]
pub(crate) fn decode_symbol(code: u8) -> u8 {
    match code {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        4 => b'R',
        5 => b'Y',
        6 => b'S',
        7 => b'W',
        8 => b'K',
        9 => b'M',
        10 => b'B',
        11 => b'D',
        12 => b'H',
        13 => b'V',
        14 => b'N',
        GAP => b'-',
        _ => unreachable!("alignment code is four bits"),
    }
}

/// Pack one split-kmer key and its IUPAC middle base into two little-endian
/// words. Split-kmer keys use at most 60 bits, leaving the high nibble free.
pub(crate) fn pack_entry_u64(kmer: u64, base: u8) -> [u32; 2] {
    debug_assert_eq!(kmer >> 60, 0);
    let packed = kmer | (u64::from(encode_symbol(base)) << 60);
    [packed as u32, (packed >> 32) as u32]
}

/// Pack one split-kmer key and its IUPAC middle base into four little-endian
/// words. Split-kmer keys use at most 124 bits, leaving the high nibble free.
pub(crate) fn pack_entry_u128(kmer: u128, base: u8) -> [u32; 4] {
    debug_assert_eq!(kmer >> 124, 0);
    let packed = kmer | (u128::from(encode_symbol(base)) << 124);
    [
        packed as u32,
        (packed >> 32) as u32,
        (packed >> 64) as u32,
        (packed >> 96) as u32,
    ]
}

#[inline]
fn should_compact(samples: usize) -> bool {
    matches!(samples, 32 | 64 | 128)
        || ((256..=MAX_SAMPLES).contains(&samples) && samples % 256 == 0)
}

#[inline]
fn eligible_for_compaction(words: &[u16], samples: usize) -> bool {
    if range_count(words) != 0 {
        return individual_count(words) != 0;
    }
    let multiplier = if samples < 256 { 2 } else { 4 };
    multiplier * group_len(words) >= samples
}

/// A final individual before the last range cannot change during recompression.
#[inline]
fn has_individual_tail(words: &[u16]) -> bool {
    let singles = individual_count(words);
    if singles == 0 {
        return false;
    }
    let ranges = range_count(words);
    ranges == 0
        || read_id(
            words,
            individual_start(words) + (singles - 1) * SAMPLE_ID_BITS,
        ) > range(words, ranges - 1).1
}

/// Number of entries in a strict lower-triangle matrix of `samples` rows.
pub(crate) fn strict_triangle_len(samples: usize) -> Option<usize> {
    if samples < 2 {
        return Some(0);
    }

    let preceding = samples - 1;
    let (left, right) = if samples.is_multiple_of(2) {
        (samples / 2, preceding)
    } else {
        (samples, preceding / 2)
    };
    left.checked_mul(right)
}

/// Index `(row, column)` in a strict lower triangle, where `row > column`.
pub(crate) fn strict_triangle_index(row: usize, column: usize) -> Option<usize> {
    if row <= column {
        return None;
    }
    strict_triangle_len(row)?.checked_add(column)
}

#[inline]
fn unpack_entry<IntT>(words: &[u32]) -> (IntT, u8)
where
    IntT: Copy + 'static,
    u128: AsPrimitive<IntT>,
{
    let mut packed = 0u128;
    for (index, word) in words.iter().enumerate() {
        packed |= u128::from(*word) << (index * 32);
    }
    let shift = words.len() * 32 - 4;
    let code = ((packed >> shift) & 0xF) as u8;
    // Packet width is derived from IntT; the masked key always fits.
    let key = (packed & ((1u128 << shift) - 1)).as_();
    (key, code)
}

/// Uncompressed upper bound, including the range/padding header.
#[inline]
fn group_words(samples: usize) -> usize {
    (PAYLOAD_BIT + SAMPLE_ID_BITS * samples).div_ceil(16)
}

/// Get the nucleotide base of a vector
#[inline]
fn group_base(words: &[u16]) -> u8 {
    (words[0] & BASE_MASK) as u8
}

/// Get the number of samples in a vector
#[inline]
fn group_len(words: &[u16]) -> usize {
    usize::from(words[0] >> 5) + 1
}

/// Get the number of ranges in a vector
#[inline]
fn range_count(words: &[u16]) -> usize {
    usize::from(words[1] & RANGE_COUNT_MASK)
}

/// Bit where all the set content of a vector ends (including the contents before the actual data, and the data itself)
#[inline]
fn payload_end(words: &[u16]) -> usize {
    words.len() * 16 - usize::from((words[1] & PADDING_MASK) >> 10)
}

/// Bit where the individual sample_ids start
#[inline]
fn individual_start(words: &[u16]) -> usize {
    PAYLOAD_BIT + 2 * SAMPLE_ID_BITS * range_count(words)
}

/// How many individual sample_ids we've got
#[inline]
fn individual_count(words: &[u16]) -> usize {
    (payload_end(words) - individual_start(words)) / SAMPLE_ID_BITS
}

/// Returns a sample_id as u16 that begins at the "bit" bit
#[inline]
fn read_id(words: &[u16], bit: usize) -> u16 {
    let word = bit / 16;
    let shift = bit % 16;
    let mut value = u32::from(words[word]);
    if shift + SAMPLE_ID_BITS > 16 {
        value |= u32::from(words[word + 1]) << 16;
    }
    ((value >> shift) as u16) & SAMPLE_ID_MASK
}

/// Writes the "id" sample_id starting at the "bit" bit
#[inline]
fn write_id(words: &mut [u16], bit: usize, id: u16) {
    debug_assert!(id <= SAMPLE_ID_MASK);
    let word = bit / 16;
    let shift = bit % 16;
    let mask = u32::from(SAMPLE_ID_MASK) << shift;
    let value = u32::from(id) << shift;
    words[word] = (words[word] & !(mask as u16)) | value as u16;
    if shift + SAMPLE_ID_BITS > 16 {
        words[word + 1] = (words[word + 1] & !((mask >> 16) as u16)) | (value >> 16) as u16;
    }
}

/// Sets the padding info
#[inline]
fn set_padding(words: &mut [u16], end: usize) {
    let padding = words.len() * 16 - end;
    debug_assert!(padding < 16);
    words[1] = (words[1] & !PADDING_MASK) | ((padding as u16) << 10);
}

/// Read the "index"th range in this vector
#[inline]
fn range(words: &[u16], index: usize) -> (u16, u16) {
    let bit = PAYLOAD_BIT + 2 * SAMPLE_ID_BITS * index;
    (read_id(words, bit), read_id(words, bit + SAMPLE_ID_BITS))
}

/// Add a sample_id to an existing vector
#[inline]
fn append_group(
    words: &mut Vec<u16>,
    sample_id: u16,
    threshold_words: usize,
    final_words: usize,
) -> Result<(), String> {
    let count = group_len(words);
    if count >= MAX_SAMPLES || usize::from(sample_id) >= MAX_SAMPLES {
        return Err("Site group exceeds the 2,048-sample limit.".into());
    }
    let ranges = range_count(words);
    if words[0] & EXTEND_MASK != 0 {
        debug_assert!(ranges > 0);
        let end_bit = PAYLOAD_BIT + (2 * ranges - 1) * SAMPLE_ID_BITS;
        if read_id(words, end_bit) + 1 == sample_id {
            write_id(words, end_bit, sample_id);
            words[0] = (words[0] & (BASE_MASK | EXTEND_MASK)) | ((count as u16) << 5);
            return Ok(());
        }
    }

    let end = payload_end(words);
    let required = (end + SAMPLE_ID_BITS).div_ceil(16);
    let capacity = words.capacity();
    let target = if ranges == 0 && words.len() >= threshold_words && capacity < final_words {
        final_words.max(required)
    } else if required > capacity {
        let geometric = capacity.saturating_mul(2).max(required);
        if ranges == 0 {
            geometric.min(final_words).max(required)
        } else {
            geometric.min(required + 32).min(final_words).max(required)
        }
    } else {
        capacity
    };
    if target > capacity {
        words
            .try_reserve_exact(target - words.len())
            .map_err(|error| format!("Could not extend site group: {error}"))?;
    }
    words.resize(required, 0);
    write_id(words, end, sample_id);
    set_padding(words, end + SAMPLE_ID_BITS);
    words[0] = (words[0] & BASE_MASK) | ((count as u16) << 5);
    Ok(())
}

/// Create new vector from a single sample_id
#[inline]
fn new_group(base: u8, sample_id: u16) -> Result<Vec<u16>, String> {
    if base > GAP || usize::from(sample_id) >= MAX_SAMPLES {
        return Err("Invalid site base or sample ID.".into());
    }
    let mut words = Vec::new();
    words
        .try_reserve_exact(group_words(1))
        .map_err(|error| format!("Could not allocate site group: {error}"))?;
    words.resize(group_words(1), 0);
    words[0] = u16::from(base);
    write_id(&mut words, PAYLOAD_BIT, sample_id);
    set_padding(&mut words, PAYLOAD_BIT + SAMPLE_ID_BITS);
    Ok(words)
}

/// Test-only decoder for checking all sample IDs independently of compression.
#[cfg(test)]
struct GroupIds<'a> {
    words: &'a [u16],
    range_index: usize,
    range_id: u16,
    individual_index: usize,
}

#[cfg(test)]
impl<'a> GroupIds<'a> {
    fn new(words: &'a [u16]) -> Self {
        Self {
            words,
            range_index: 0,
            range_id: if range_count(words) == 0 {
                0
            } else {
                range(words, 0).0
            },
            individual_index: 0,
        }
    }
}

#[cfg(test)]
impl Iterator for GroupIds<'_> {
    type Item = u16;

    fn next(&mut self) -> Option<u16> {
        let ranged = (self.range_index < range_count(self.words)).then_some(self.range_id);
        let single = (self.individual_index < individual_count(self.words)).then(|| {
            read_id(
                self.words,
                individual_start(self.words) + self.individual_index * SAMPLE_ID_BITS,
            )
        });
        if let Some(id) = ranged {
            if single.is_none_or(|single| id < single) {
                let end = range(self.words, self.range_index).1;
                if id == end {
                    self.range_index += 1;
                    if self.range_index < range_count(self.words) {
                        self.range_id = range(self.words, self.range_index).0;
                    }
                } else {
                    self.range_id += 1;
                }
                return Some(id);
            }
        }
        if single.is_some() {
            self.individual_index += 1;
        }
        single
    }
}

// This struct is used for storing a fixed size to be used for ALL split k-mers, and allocating it only once
struct CompactionScratch {
    ids: [u16; MAX_SAMPLES],
    encoded: [u16; MAX_GROUP_WORDS],
}

impl CompactionScratch {
    fn new() -> Self {
        Self {
            ids: [0; MAX_SAMPLES],
            encoded: [0; MAX_GROUP_WORDS],
        }
    }

    fn encode(&mut self, words: &[u16]) -> usize {
        let count = group_len(words);
        let old_ranges = range_count(words);
        let individuals = individual_count(words);
        for index in 0..individuals {
            self.ids[index] = read_id(words, individual_start(words) + index * SAMPLE_ID_BITS);
            debug_assert!(index == 0 || self.ids[index - 1] < self.ids[index]);
        }

        // Earlier singletons cannot gain missing IDs: new samples arrive in order.
        let mut last_end = if old_ranges == 0 {
            None
        } else {
            Some(range(words, old_ranges - 1).1)
        };
        let prefix = last_end.map_or(0, |end| {
            self.ids[..individuals].partition_point(|id| *id < end)
        });
        let mut tail = prefix;
        if let Some(end) = last_end.as_mut() {
            while tail < individuals && self.ids[tail] == *end + 1 {
                *end = self.ids[tail];
                tail += 1;
            }
        }

        let mut ranges = old_ranges;
        let mut singles = prefix;
        let mut index = tail;
        let mut ends_in_range = old_ranges != 0 && tail == individuals;
        while index < individuals {
            let start = index;
            index += 1;
            while index < individuals && self.ids[index] == self.ids[index - 1] + 1 {
                index += 1;
            }
            ends_in_range = index - start >= MIN_RANGE_SAMPLES;
            if ends_in_range {
                ranges += 1;
            } else {
                singles += index - start;
            }
        }
        debug_assert!(ranges <= usize::from(RANGE_COUNT_MASK));
        let end = PAYLOAD_BIT + SAMPLE_ID_BITS * (2 * ranges + singles);
        let length = end.div_ceil(16);
        self.encoded[..length].fill(0);

        // Copy the packed range prefix without decoding its IDs or endpoints.
        let prefix_end = individual_start(words);
        let whole_words = prefix_end / 16;
        self.encoded[..whole_words].copy_from_slice(&words[..whole_words]);
        if prefix_end % 16 != 0 {
            self.encoded[whole_words] = words[whole_words] & ((1u16 << (prefix_end % 16)) - 1);
        }
        self.encoded[0] = u16::from(group_base(words)) | (((count - 1) as u16) << 5);
        self.encoded[1] = (self.encoded[1] & !(RANGE_COUNT_MASK | PADDING_MASK)) | ranges as u16;
        if let Some(last_end) = last_end {
            write_id(
                &mut self.encoded[..length],
                PAYLOAD_BIT + (2 * old_ranges - 1) * SAMPLE_ID_BITS,
                last_end,
            );
        }
        let mut range_bit = prefix_end;
        let mut single_bit = PAYLOAD_BIT + 2 * SAMPLE_ID_BITS * ranges;
        for id in &self.ids[..prefix] {
            write_id(&mut self.encoded[..length], single_bit, *id);
            single_bit += SAMPLE_ID_BITS;
        }
        index = tail;
        while index < individuals {
            let start = index;
            index += 1;
            while index < individuals && self.ids[index] == self.ids[index - 1] + 1 {
                index += 1;
            }
            if index - start >= MIN_RANGE_SAMPLES {
                write_id(&mut self.encoded[..length], range_bit, self.ids[start]);
                write_id(
                    &mut self.encoded[..length],
                    range_bit + SAMPLE_ID_BITS,
                    self.ids[index - 1],
                );
                range_bit += 2 * SAMPLE_ID_BITS;
            } else {
                for id in &self.ids[start..index] {
                    write_id(&mut self.encoded[..length], single_bit, *id);
                    single_bit += SAMPLE_ID_BITS;
                }
            }
        }
        if ends_in_range {
            self.encoded[0] |= EXTEND_MASK;
        }
        set_padding(&mut self.encoded[..length], end);
        length
    }
}

/// Reserve before swapping, so an allocation failure leaves the live group intact.
fn compact_group_with<F>(
    words: &mut Vec<u16>,
    scratch: &mut CompactionScratch,
    reserve: F,
) -> Result<bool, ()>
where
    F: FnOnce(usize) -> Result<Vec<u16>, ()>,
{
    let length = scratch.encode(words);
    if length * 4 > words.len() * 3 {
        return Ok(false);
    }
    let mut replacement = reserve(length)?;
    if replacement.capacity() < length {
        return Err(());
    }
    replacement.extend_from_slice(&scratch.encoded[..length]);
    *words = replacement;
    Ok(true)
}

/// Progress required by the website during a compression scan.
#[derive(Debug, Clone, Default)]
pub(crate) struct CompactionProgress {
    pub processed_samples: usize,
    pub sites_total: usize,
    pub sites_processed: usize,
    pub complete: bool,
}

// Calculation of distances
#[inline]
fn accumulate_group_distances(group: &[u16], distance_row: &mut [u32]) {
    // Groups partition prior samples; every differing call contributes exactly once.
    for index in 0..range_count(group) {
        let (start, end) = range(group, index);
        for distance in &mut distance_row[usize::from(start)..=usize::from(end)] {
            *distance += 1;
        }
    }
    let start = individual_start(group);
    let count = individual_count(group);
    let full = count / 16 * 16;
    for block in (0..full).step_by(16) {
        // This might help vectorisation (?)
        for index in 0..16 {
            let id = read_id(group, start + (block + index) * SAMPLE_ID_BITS);
            distance_row[usize::from(id)] += 1;
        }
    }
    for index in full..count {
        let id = read_id(group, start + index * SAMPLE_ID_BITS);
        distance_row[usize::from(id)] += 1;
    }
}

/// A site while samples are still being added.
#[derive(Debug, Clone)]
enum GrowingSite {
    Monomorphic(Vec<u16>),
    Polymorphic(Vec<Vec<u16>>),
}

/// When extracting alignment
#[derive(Debug, Clone)]
struct FrozenGroup {
    words: Vec<u16>,
    range_cursor: u16,
    individual_cursor: u16,
}

// The second u16 cursor fits the old padding on native, WASM32, and WASM64.
const _: () = assert!(std::mem::size_of::<FrozenGroup>() == std::mem::size_of::<(Vec<u16>, u16)>());

impl FrozenGroup {
    fn new(words: Vec<u16>) -> Self {
        Self {
            words,
            range_cursor: 0,
            individual_cursor: 0,
        }
    }

    fn reset(&mut self) {
        self.range_cursor = 0;
        self.individual_cursor = 0;
    }

    fn take_sample(&mut self, sample_id: u16) -> Option<u8> {
        while usize::from(self.range_cursor) < range_count(&self.words) {
            let (start, end) = range(&self.words, usize::from(self.range_cursor));
            if end < sample_id {
                self.range_cursor += 1;
            } else {
                if start <= sample_id {
                    return Some(group_base(&self.words));
                }
                break;
            }
        }
        let start = individual_start(&self.words);
        while usize::from(self.individual_cursor) < individual_count(&self.words) {
            let id = read_id(
                &self.words,
                start + usize::from(self.individual_cursor) * SAMPLE_ID_BITS,
            );
            if id > sample_id {
                break;
            }
            self.individual_cursor += 1;
            if id == sample_id {
                return Some(group_base(&self.words));
            }
        }
        None
    }
}

#[derive(Debug, Clone)]
enum FrozenSite {
    Monomorphic(FrozenGroup),
    Polymorphic { start: usize, len: usize },
}

#[derive(Debug, Clone)]
struct FrozenMatrix {
    sites: Vec<FrozenSite>,
    polymorphic_groups: Vec<FrozenGroup>,
}

impl FrozenMatrix {
    fn try_with_capacity(site_count: usize, group_count: usize) -> Result<Self, String> {
        let mut sites = Vec::new();
        sites
            .try_reserve_exact(site_count)
            .map_err(|error| format!("Could not allocate frozen site array: {error}"))?;

        let mut polymorphic_groups = Vec::new();
        polymorphic_groups
            .try_reserve_exact(group_count)
            .map_err(|error| format!("Could not allocate frozen group array: {error}"))?;

        Ok(Self {
            sites,
            polymorphic_groups,
        })
    }
}

#[derive(Debug, Clone)]
enum SiteStorage<IntT> {
    Growing(HashMap<IntT, GrowingSite>),
    Frozen(FrozenMatrix),
    Released,
}

#[derive(Debug, Clone)]
pub(crate) struct AlignmentMatrix<IntT> {
    storage: SiteStorage<IntT>,
    sample_count: usize,
    sealed: bool,
    /// Packed-word threshold for switching to the planned final capacity.
    growth_threshold_words: usize,
    final_group_words: usize,
}

impl<IntT> Default for AlignmentMatrix<IntT> {
    fn default() -> Self {
        Self {
            storage: SiteStorage::Growing(HashMap::new()),
            sample_count: 0,
            sealed: false,
            growth_threshold_words: 1024,
            final_group_words: group_words(MAX_SAMPLES),
        }
    }
}

impl<IntT> AlignmentMatrix<IntT>
where
    IntT: Copy + Eq + Hash + PrimInt + 'static,
    u128: AsPrimitive<IntT>,
{
    pub(crate) fn is_sealed(&self) -> bool {
        self.sealed
    }

    /// Prepare group capacities once for a batch whose total was validated by plan_samples.
    /// This does NOT remove the single check whenever adding samples, that prevents allocating
    /// more memory than needed.
    pub(crate) fn set_expected_samples(&mut self, total: usize) -> Result<(), String> {
        if total == 0 || total > MAX_SAMPLES || total < self.sample_count {
            return Err("Invalid planned sample count.".into());
        }
        let final_words = group_words(total);
        let threshold_words = 1usize << (final_words - 1).ilog2();
        if self.sample_count != 0 {
            let SiteStorage::Growing(postings) = &mut self.storage else {
                return Err("The alignment is closed to new samples.".into());
            };
            let prepare = |words: &mut Vec<u16>| -> Result<(), String> {
                if range_count(words) != 0 {
                    return Ok(());
                }
                let capacity = words.capacity();
                let target = if words.len() >= threshold_words {
                    final_words
                } else if capacity < threshold_words {
                    capacity.next_power_of_two().min(threshold_words)
                } else {
                    capacity
                };
                if target > capacity {
                    words
                        .try_reserve_exact(target - words.len())
                        .map_err(|error| format!("Could not prepare site group: {error}"))?;
                }
                Ok(())
            };
            for site in postings.values_mut() {
                match site {
                    GrowingSite::Monomorphic(words) => prepare(words)?,
                    GrowingSite::Polymorphic(groups) => {
                        for words in groups {
                            prepare(words)?;
                        }
                    }
                }
            }
        }
        self.growth_threshold_words = threshold_words;
        self.final_group_words = final_words;
        Ok(())
    }

    /// Store one sample and write its distances to earlier samples into `distance_row`,
    /// visiting each split k-mer once. A failure after storing has begun seals the matrix.
    pub(crate) fn add_sample_words(
        &mut self,
        words: Vec<u32>,
        distance_row: &mut [u32],
        on_compaction: &dyn Fn(&CompactionProgress),
    ) -> Result<(), String> {
        if self.sealed {
            return Err("The alignment is closed to new samples.".into());
        }
        if distance_row.len() != self.sample_count {
            return Err("Distance row length does not match the existing sample count.".into());
        }
        if self.sample_count >= MAX_SAMPLES {
            return Err("Sparse alignment supports at most 2,048 samples.".into());
        }
        let width = std::mem::size_of::<IntT>() / std::mem::size_of::<u32>();
        if words.len() % width != 0 {
            return Err("Incomplete packed split-kmer entry.".into());
        }
        let result = self.store_and_compare(&words, distance_row);
        if result.is_err() {
            self.sealed = true;
        } else if should_compact(self.sample_count) {
            self.compact_groups(on_compaction);
        }
        result
    }

    fn compact_groups(&mut self, on_progress: &dyn Fn(&CompactionProgress)) {
        self.compact_groups_with(on_progress, &mut |length| {
            let mut replacement = Vec::new();
            replacement.try_reserve_exact(length).map_err(|_| ())?;
            Ok(replacement)
        });
    }

    /// Final pass uses a strict quarter-support threshold at every sample count.
    pub(crate) fn compact_final(&mut self, on_progress: &dyn Fn(&CompactionProgress)) {
        self.compact_groups_with_policy(
            on_progress,
            &mut |length| {
                let mut replacement = Vec::new();
                replacement.try_reserve_exact(length).map_err(|_| ())?;
                Ok(replacement)
            },
            &|words, samples| 4 * group_len(words) > samples && individual_count(words) != 0,
        );
    }

    fn compact_groups_with<R>(&mut self, on_progress: &dyn Fn(&CompactionProgress), reserve: &mut R)
    where
        R: FnMut(usize) -> Result<Vec<u16>, ()>,
    {
        self.compact_groups_with_policy(on_progress, reserve, &eligible_for_compaction);
    }

    fn compact_groups_with_policy<R>(
        &mut self,
        on_progress: &dyn Fn(&CompactionProgress),
        reserve: &mut R,
        eligible: &dyn Fn(&[u16], usize) -> bool,
    ) where
        R: FnMut(usize) -> Result<Vec<u16>, ()>,
    {
        let SiteStorage::Growing(postings) = &mut self.storage else {
            return;
        };
        let mut progress = CompactionProgress {
            processed_samples: self.sample_count,
            sites_total: postings.len(),
            ..CompactionProgress::default()
        };
        on_progress(&progress);
        let mut scratch = CompactionScratch::new();
        for site in postings.values_mut() {
            let mut visit = |words: &mut Vec<u16>| {
                if has_individual_tail(words) && eligible(words, self.sample_count) {
                    // A failed replacement allocation leaves the original vector intact.
                    let _ = compact_group_with(words, &mut scratch, &mut *reserve);
                }
            };
            match site {
                GrowingSite::Monomorphic(words) => visit(words),
                GrowingSite::Polymorphic(groups) => {
                    for words in groups {
                        visit(words);
                    }
                }
            }
            progress.sites_processed += 1;
            if progress.sites_processed % 65_536 == 0 {
                on_progress(&progress);
            }
        }
        progress.complete = true;
        on_progress(&progress);
    }

    fn store_and_compare(&mut self, words: &[u32], distance_row: &mut [u32]) -> Result<(), String> {
        let SiteStorage::Growing(postings) = &mut self.storage else {
            return Err("The alignment is closed to new samples.".into());
        };
        let sample_id = self.sample_count as u16;
        let words_per_entry = std::mem::size_of::<IntT>() / std::mem::size_of::<u32>();
        for entry in words.chunks_exact(words_per_entry) {
            let (kmer, code) = unpack_entry::<IntT>(entry);
            match postings.entry(kmer) {
                Entry::Occupied(mut entry) => {
                    let site = entry.get_mut();
                    match site {
                        GrowingSite::Monomorphic(existing) if group_base(existing) == code => {
                            append_group(
                                existing,
                                sample_id,
                                self.growth_threshold_words,
                                self.final_group_words,
                            )?;
                        }
                        GrowingSite::Monomorphic(existing) => {
                            let new_words = new_group(code, sample_id)?;

                            let mut groups = Vec::new(); // This should be equivalent to with_capacity
                            groups.try_reserve_exact(2).map_err(|error| {
                                format!("Could not allocate site groups: {error}")
                            })?;

                            // Native defaults ignore ambiguous and missing calls in either sample.
                            if is_concrete_code(code) && is_concrete_code(group_base(existing)) {
                                accumulate_group_distances(existing, distance_row);
                            }

                            // Move data out, without recopying
                            let old_site =
                                std::mem::replace(site, GrowingSite::Monomorphic(Vec::new()));
                            let GrowingSite::Monomorphic(old_words) = old_site else {
                                unreachable!()
                            };
                            groups.push(old_words);
                            groups.push(new_words);
                            *site = GrowingSite::Polymorphic(groups);
                        }
                        GrowingSite::Polymorphic(groups) => {
                            let mut found = false;
                            for group in groups.iter_mut() {
                                let base = group_base(group);
                                if base == code {
                                    append_group(
                                        group,
                                        sample_id,
                                        self.growth_threshold_words,
                                        self.final_group_words,
                                    )?;
                                    found = true;
                                } else if is_concrete_code(code) && is_concrete_code(base) {
                                    accumulate_group_distances(group, distance_row);
                                }
                            }
                            if !found {
                                let new_words = new_group(code, sample_id)?;
                                groups.try_reserve(1).map_err(|error| {
                                    format!("Could not extend site groups: {error}")
                                })?;
                                groups.push(new_words);
                            }
                        }
                    }
                }
                Entry::Vacant(site) => {
                    site.insert(GrowingSite::Monomorphic(new_group(code, sample_id)?));
                }
            }
        }
        self.sample_count += 1;
        Ok(())
    }

    /// Move packed groups without copying their payloads, then discard k-mer keys.
    fn freeze(&mut self) -> Result<(), String> {
        self.freeze_with(FrozenMatrix::try_with_capacity)
    }

    fn freeze_with<F>(&mut self, reserve: F) -> Result<(), String>
    where
        F: FnOnce(usize, usize) -> Result<FrozenMatrix, String>,
    {
        self.sealed = true;
        let (site_count, group_count) = match &self.storage {
            SiteStorage::Growing(postings) => {
                let mut group_count = 0usize;
                for site in postings.values() {
                    if let GrowingSite::Polymorphic(groups) = site {
                        group_count = group_count
                            .checked_add(groups.len())
                            .ok_or_else(|| "Frozen group count overflowed.".to_string())?;
                    }
                }
                (postings.len(), group_count)
            }
            SiteStorage::Frozen(_) => return Ok(()),
            SiteStorage::Released => {
                return Err("Alignment sites have already been exported and released.".into())
            }
        };

        // Complete all fallible reservations while the original map remains available.
        let mut frozen = reserve(site_count, group_count)?;
        let SiteStorage::Growing(postings) =
            std::mem::replace(&mut self.storage, SiteStorage::Released)
        else {
            unreachable!()
        };

        // Both vectors were reserved to their final lengths above; moving the groups cannot allocate.
        for (_key, site) in postings {
            match site {
                GrowingSite::Monomorphic(words) => {
                    frozen
                        .sites
                        .push(FrozenSite::Monomorphic(FrozenGroup::new(words)));
                }
                GrowingSite::Polymorphic(groups) => {
                    let start = frozen.polymorphic_groups.len();
                    let len = groups.len();
                    for words in groups {
                        frozen.polymorphic_groups.push(FrozenGroup::new(words));
                    }
                    frozen.sites.push(FrozenSite::Polymorphic { start, len });
                }
            }
        }

        self.storage = SiteStorage::Frozen(frozen);
        Ok(())
    }

    /// Reset cursors on every attempt, permitting retries after writer errors.
    fn write_fasta<W: Write>(
        &mut self,
        names: &[String],
        writer: &mut W,
        on_sample: &dyn Fn(usize),
    ) -> io::Result<()> {
        let SiteStorage::Frozen(frozen) = &mut self.storage else {
            return Err(io::Error::other("The sparse array is not available."));
        };
        if names.len() != self.sample_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Sample names do not match the matrix.",
            ));
        }
        for site in &mut frozen.sites {
            match site {
                FrozenSite::Monomorphic(group) => group.reset(),
                FrozenSite::Polymorphic { start, len } => {
                    for group in &mut frozen.polymorphic_groups[*start..(*start + *len)] {
                        group.reset();
                    }
                }
            }
        }
        let mut buffer = [0u8; 8192];
        for (sample, name) in names.iter().enumerate() {
            writer.write_all(b"> ")?;
            writer.write_all(name.as_bytes())?;
            writer.write_all(b"\n")?;
            let mut buffered = 0;
            for site in &mut frozen.sites {
                let code = match site {
                    FrozenSite::Monomorphic(group) => {
                        group.take_sample(sample as u16).unwrap_or(GAP)
                    }
                    FrozenSite::Polymorphic { start, len } => frozen.polymorphic_groups
                        [*start..(*start + *len)]
                        .iter_mut()
                        .find_map(|group| group.take_sample(sample as u16))
                        .unwrap_or(GAP),
                };
                buffer[buffered] = decode_symbol(code);
                buffered += 1;
                if buffered == buffer.len() {
                    writer.write_all(&buffer)?;
                    buffered = 0;
                }
            }
            if buffered > 0 {
                writer.write_all(&buffer[..buffered])?;
            }
            writer.write_all(b"\n")?;
            on_sample(sample + 1);
        }
        Ok(())
    }

    pub(crate) fn export_fasta_gzip(
        &mut self,
        names: &[String],
        on_sample: &dyn Fn(usize),
    ) -> Result<Vec<u8>, String> {
        self.freeze()?;
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        self.write_fasta(names, &mut encoder, on_sample)
            .map_err(|error| format!("Could not write alignment: {error}"))?;
        let bytes = encoder
            .finish()
            .map_err(|error| format!("Could not gzip alignment: {error}"))?;
        self.storage = SiteStorage::Released;
        Ok(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    fn visit_test_groups<IntT>(matrix: &AlignmentMatrix<IntT>, mut visit: impl FnMut(&Vec<u16>)) {
        match &matrix.storage {
            SiteStorage::Growing(sites) => {
                for site in sites.values() {
                    match site {
                        GrowingSite::Monomorphic(words) => visit(words),
                        GrowingSite::Polymorphic(groups) => {
                            for words in groups {
                                visit(words);
                            }
                        }
                    }
                }
            }
            SiteStorage::Frozen(frozen) => {
                for site in &frozen.sites {
                    if let FrozenSite::Monomorphic(group) = site {
                        visit(&group.words);
                    }
                }
                for group in &frozen.polymorphic_groups {
                    visit(&group.words);
                }
            }
            SiteStorage::Released => {}
        }
    }

    fn packed_capacity<IntT>(matrix: &AlignmentMatrix<IntT>) -> usize {
        let mut capacity = 0;
        visit_test_groups(matrix, |words| capacity += words.capacity());
        capacity
    }

    fn contains_ranges<IntT>(matrix: &AlignmentMatrix<IntT>) -> bool {
        let mut found = false;
        visit_test_groups(matrix, |words| found |= range_count(words) != 0);
        found
    }

    fn compact(words: &mut Vec<u16>) -> bool {
        compact_group_with(words, &mut CompactionScratch::new(), |length| {
            let mut replacement = Vec::new();
            replacement.try_reserve_exact(length).map_err(|_| ())?;
            Ok(replacement)
        })
        .unwrap()
    }

    fn group(ids: &[u16]) -> Vec<u16> {
        let mut words = new_group(4, ids[0]).unwrap();
        for id in &ids[1..] {
            append_group(&mut words, *id, 1024, group_words(MAX_SAMPLES)).unwrap();
        }
        words
    }

    #[test]
    fn fields_clear_previous_bits_at_every_offset() {
        for offset in 0..32 {
            let mut words = [u16::MAX; 4];
            for id in [2047, 255, 256, 1023, 1024, 0] {
                write_id(&mut words, offset, id);
                assert_eq!(read_id(&words, offset), id);
                for bit in 0..64 {
                    if !(offset..offset + 11).contains(&bit) {
                        assert_ne!(words[bit / 16] & (1 << (bit % 16)), 0);
                    }
                }
            }
        }
    }

    #[test]
    fn uncompressed_groups_cover_padding_and_maximum_count() {
        let mut words = new_group(4, 0).unwrap();
        let mut paddings = [false; 16];
        assert_eq!(words.len(), 3);
        for id in 1..MAX_SAMPLES {
            append_group(&mut words, id as u16, 1024, group_words(MAX_SAMPLES)).unwrap();
            let padding = (words[1] & PADDING_MASK) >> 10;
            paddings[usize::from(padding)] = true;
            assert_eq!(group_len(&words), id + 1);
            assert_eq!(read_id(&words, payload_end(&words) - 11), id as u16);
            if padding != 0 {
                assert_eq!(words.last().unwrap() >> (16 - padding), 0);
            }
        }
        assert!(paddings.iter().all(|seen| *seen));
        assert_eq!(
            GroupIds::new(&words).collect::<Vec<_>>(),
            (0..2048).collect::<Vec<u16>>()
        );
        let snapshot = words.clone();
        assert!(append_group(&mut words, 2047, 1024, group_words(MAX_SAMPLES)).is_err());
        assert_eq!(words, snapshot);
        assert!(new_group(0, 2048).is_err());
    }

    #[test]
    fn ranges_require_three_samples_and_merge_with_individuals() {
        let ids: Vec<u16> = [0, 1]
            .into_iter()
            .chain(4..90)
            .chain([94, 95, 99])
            .chain(110..200)
            .chain([205, 207])
            .collect();
        let mut words = group(&ids);
        assert!(compact(&mut words));
        assert_eq!(range_count(&words), 2);
        assert_eq!(range(&words, 0), (4, 89));
        assert_eq!(range(&words, 1), (110, 199));
        assert_eq!(individual_count(&words), 7);
        assert_eq!(GroupIds::new(&words).collect::<Vec<_>>(), ids);
        assert_eq!(words[0] & EXTEND_MASK, 0);
        let encoded = words.clone();
        assert!(!compact(&mut words));
        assert_eq!(words, encoded);

        let two = group(&[0, 1]);
        let three = group(&[0, 1, 2]);
        let mut scratch = CompactionScratch::new();
        let length = scratch.encode(&two);
        assert_eq!(range_count(&scratch.encoded[..length]), 0);
        let length = scratch.encode(&three);
        assert_eq!(range(&scratch.encoded[..length], 0), (0, 2));
    }

    #[test]
    fn extension_overwrites_boundaries_and_recompaction_reactivates_it() {
        let mut words = group(&(0..255).collect::<Vec<_>>());
        assert!(compact(&mut words));
        let capacity = words.capacity();
        for id in 255..=256 {
            append_group(&mut words, id, 128, group_words(2048)).unwrap();
        }
        assert_eq!(range(&words, 0), (0, 256));
        assert_eq!(words.capacity(), capacity);
        append_group(&mut words, 300, 128, group_words(2048)).unwrap();
        assert_eq!(words[0] & EXTEND_MASK, 0);
        for id in 301..400 {
            append_group(&mut words, id, 128, group_words(2048)).unwrap();
        }
        assert!(compact(&mut words));
        assert_ne!(words[0] & EXTEND_MASK, 0);
        let capacity = words.capacity();
        for id in 400..=1024 {
            append_group(&mut words, id, 128, group_words(2048)).unwrap();
        }
        assert_eq!(range(&words, 1), (300, 1024));
        assert_eq!(words.capacity(), capacity);
        for id in 1025..=2047 {
            append_group(&mut words, id, 128, group_words(2048)).unwrap();
        }
        assert_eq!(range(&words, 1), (300, 2047));
        assert!(GroupIds::new(&words).is_sorted());
    }

    #[test]
    fn failed_compaction_keeps_group_and_allocations_unchanged() {
        let mut words = group(&(0..255).collect::<Vec<_>>());
        let old = words.clone();
        let capacity = words.capacity();
        assert_eq!(
            compact_group_with(&mut words, &mut CompactionScratch::new(), |_| Err(())),
            Err(())
        );
        assert_eq!(words, old);
        assert_eq!(words.capacity(), capacity);
        append_group(&mut words, 255, 128, group_words(2048)).unwrap();
        assert_eq!(group_len(&words), 256);
    }

    #[test]
    fn distances_visit_all_range_members_and_singletons_once() {
        let ids: Vec<u16> = (0..90).chain([94, 95, 99]).chain(110..200).collect();
        let mut words = group(&ids);
        assert!(compact(&mut words));
        let mut row = vec![7; 200];
        accumulate_group_distances(&words, &mut row);
        for (id, value) in row.iter().enumerate() {
            assert_eq!(*value, 7 + u32::from(ids.contains(&(id as u16))));
        }
    }

    #[test]
    fn frozen_cursors_handle_skips_and_fit_existing_padding() {
        struct OldGroup {
            _words: Vec<u16>,
            _cursor: u16,
        }
        assert_eq!(
            std::mem::size_of::<FrozenGroup>(),
            std::mem::size_of::<OldGroup>()
        );
        let ids: Vec<u16> = [0, 1]
            .into_iter()
            .chain(4..90)
            .chain([94, 95, 99])
            .chain(110..200)
            .chain([205, 207])
            .collect();
        let mut words = group(&ids);
        assert!(compact(&mut words));
        let mut frozen = FrozenGroup::new(words);
        for id in [0, 2, 5, 90, 95, 100, 180, 205, 206, 207, 208] {
            assert_eq!(frozen.take_sample(id), ids.contains(&id).then_some(4));
        }
        frozen.reset();
        for id in 0..210 {
            assert_eq!(frozen.take_sample(id), ids.contains(&id).then_some(4));
        }
    }

    #[test]
    fn checkpoint_selection_and_batch_reservation() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.sample_count = 256;
        let SiteStorage::Growing(sites) = &mut matrix.storage else {
            unreachable!()
        };
        sites.insert(
            0,
            GrowingSite::Monomorphic(group(&(0..63).collect::<Vec<_>>())),
        );
        sites.insert(
            1,
            GrowingSite::Monomorphic(group(&(0..64).collect::<Vec<_>>())),
        );
        let events = RefCell::new(Vec::new());
        matrix.compact_groups(&|p| events.borrow_mut().push(p.clone()));
        assert!(events.borrow().last().unwrap().complete);
        let SiteStorage::Growing(sites) = &matrix.storage else {
            unreachable!()
        };
        let GrowingSite::Monomorphic(compacted) = &sites[&1] else {
            unreachable!()
        };
        assert_eq!(range_count(compacted), 1);
        let GrowingSite::Monomorphic(ineligible) = &sites[&0] else {
            unreachable!()
        };
        assert_eq!(range_count(ineligible), 0);
        let capacity = compacted.capacity();
        matrix.set_expected_samples(2048).unwrap();
        let SiteStorage::Growing(sites) = &matrix.storage else {
            unreachable!()
        };
        let GrowingSite::Monomorphic(compacted) = &sites[&1] else {
            unreachable!()
        };
        assert_eq!(compacted.capacity(), capacity);
    }

    #[test]
    fn compaction_allocation_failure_preserves_vectors_without_sealing() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.sample_count = 255;
        let SiteStorage::Growing(sites) = &mut matrix.storage else {
            unreachable!()
        };
        sites.insert(
            0,
            GrowingSite::Monomorphic(group(&(0..129).collect::<Vec<_>>())),
        );
        let old_capacity = packed_capacity(&matrix);
        let events = RefCell::new(Vec::new());
        matrix.compact_groups_with(&|p| events.borrow_mut().push(p.clone()), &mut |_| Err(()));
        assert!(events.borrow().last().unwrap().complete);
        assert_eq!(packed_capacity(&matrix), old_capacity);
        assert!(!contains_ranges(&matrix));
        assert!(!matrix.is_sealed());
        let mut row = vec![0; 255];
        matrix
            .add_sample_words(pack_entry_u64(0, b'R').to_vec(), &mut row, &|_| {})
            .unwrap();
        assert_eq!(matrix.sample_count, 256);
        assert!(row.iter().all(|distance| *distance == 0));
        assert!(!matrix.is_sealed());
    }

    #[test]
    fn scans_report_bounded_progress_and_skip_canonical_range_only_groups() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.sample_count = 511;
        let SiteStorage::Growing(sites) = &mut matrix.storage else {
            unreachable!()
        };
        for key in 0..65_537 {
            sites.insert(key, GrowingSite::Monomorphic(new_group(0, 0).unwrap()));
        }
        let mut dense = group(&(0..511).collect::<Vec<_>>());
        assert!(compact(&mut dense));
        sites.insert(65_537, GrowingSite::Monomorphic(dense));
        let events = RefCell::new(Vec::new());
        matrix.compact_groups(&|p| events.borrow_mut().push(p.clone()));
        let events = events.borrow();
        assert_eq!(
            events.iter().map(|p| p.sites_processed).collect::<Vec<_>>(),
            [0, 65_536, 65_538]
        );
        assert!(events.last().unwrap().complete);
        assert_eq!(events.last().unwrap().sites_total, 65_538);
        assert!(contains_ranges(&matrix));
    }

    #[test]
    fn ranged_groups_grow_with_at_most_32_spare_words() {
        let mut words = group(&(0..100).collect::<Vec<_>>());
        assert!(compact(&mut words));
        for id in (120..2048).step_by(2) {
            append_group(&mut words, id, 128, group_words(2048)).unwrap();
            assert!(words.capacity() <= words.len() + 32);
        }
        assert_eq!(range(&words, 0), (0, 99));
        assert_eq!(individual_count(&words), 964);
    }

    #[test]
    fn checkpoints_include_every_256_samples_after_the_early_passes() {
        let checkpoints: Vec<_> = (0..=MAX_SAMPLES + 256)
            .filter(|n| should_compact(*n))
            .collect();
        assert_eq!(
            checkpoints,
            [32, 64, 128, 256, 512, 768, 1024, 1280, 1536, 1792, 2048]
        );
    }

    #[test]
    fn first_compression_support_thresholds_are_inclusive() {
        for (samples, minimum) in [
            (32, 16),
            (64, 32),
            (128, 64),
            (255, 128),
            (256, 64),
            (512, 128),
        ] {
            assert!(!eligible_for_compaction(
                &group(&(0..minimum - 1).collect::<Vec<_>>()),
                samples
            ));
            assert!(eligible_for_compaction(
                &group(&(0..minimum).collect::<Vec<_>>()),
                samples
            ));
        }
    }

    #[test]
    fn recompression_accepts_small_capacity_savings_below_support_threshold() {
        let mut words = group(&(0..16).collect::<Vec<_>>());
        assert!(compact(&mut words));
        for id in 20..25 {
            append_group(&mut words, id, 128, group_words(MAX_SAMPLES)).unwrap();
        }
        assert!(4 * group_len(&words) < 256);
        assert!(eligible_for_compaction(&words, 256));
        let ids = GroupIds::new(&words).collect::<Vec<_>>();
        let old_length = words.len();
        let old_capacity = words.capacity();
        assert!(compact(&mut words));
        let saved_bytes = (old_capacity - words.capacity()) * 2;
        assert!(saved_bytes > 0 && saved_bytes < 16);
        assert!(words.len() * 4 <= old_length * 3);
        assert_eq!(range(&words, 0), (0, 15));
        assert_eq!(range(&words, 1), (20, 24));
        assert_eq!(GroupIds::new(&words).collect::<Vec<_>>(), ids);
        assert!(!eligible_for_compaction(&words, 256));
        append_group(&mut words, 25, 128, group_words(MAX_SAMPLES)).unwrap();
        assert_eq!(range(&words, 1), (20, 25));
    }

    #[test]
    fn recompression_matches_full_encoding_without_expanding_old_ranges() {
        for offset in 0..16 {
            let mut ids: Vec<u16> = (offset..offset + 80)
                .chain([offset + 85, offset + 86])
                .chain(offset + 90..offset + 140)
                .chain([offset + 145, offset + 147])
                .collect();
            let mut words = group(&ids);
            assert!(compact(&mut words));
            for id in offset + 150..2048 {
                if id % 13 < 9 {
                    append_group(&mut words, id, 128, group_words(MAX_SAMPLES)).unwrap();
                    ids.push(id);
                }
                if id % 31 == 0 || id == 2047 {
                    let mut incremental = CompactionScratch::new();
                    let length = incremental.encode(&words);
                    let mut full = CompactionScratch::new();
                    let expected_length = full.encode(&group(&ids));
                    assert_eq!(length, expected_length);
                    assert_eq!(incremental.encoded[..length], full.encoded[..length]);
                    compact(&mut words);
                    assert_eq!(GroupIds::new(&words).collect::<Vec<_>>(), ids);
                    assert_eq!(range(&words, 0), (offset, offset + 79));
                    assert_eq!(range(&words, 1), (offset + 90, offset + 139));
                }
            }
        }
    }

    #[test]
    fn recompression_extends_only_last_range_and_preserves_short_runs() {
        // A valid packed group with trailing singletons touching its final range.
        let ranges = [(10, 19), (30, 39)];
        let singles: Vec<u16> = [0, 1, 25]
            .into_iter()
            .chain(40..60)
            .chain([65, 66])
            .collect();
        let count = 20 + singles.len();
        let end = PAYLOAD_BIT + SAMPLE_ID_BITS * (4 + singles.len());
        let mut words = vec![0; end.div_ceil(16)];
        words[0] = 4 | (((count - 1) as u16) << 5);
        words[1] = 2;
        for (index, (start, end)) in ranges.into_iter().enumerate() {
            let bit = PAYLOAD_BIT + 2 * SAMPLE_ID_BITS * index;
            write_id(&mut words, bit, start);
            write_id(&mut words, bit + SAMPLE_ID_BITS, end);
        }
        for (index, id) in singles.iter().enumerate() {
            write_id(&mut words, PAYLOAD_BIT + (4 + index) * SAMPLE_ID_BITS, *id);
        }
        set_padding(&mut words, end);
        let ids = GroupIds::new(&words).collect::<Vec<_>>();
        assert!(compact(&mut words));
        assert_eq!(range_count(&words), 2);
        assert_eq!(range(&words, 0), (10, 19));
        assert_eq!(range(&words, 1), (30, 59));
        assert_eq!(individual_count(&words), 5);
        assert_eq!(words[0] & EXTEND_MASK, 0);
        assert_eq!(GroupIds::new(&words).collect::<Vec<_>>(), ids);
        let snapshot = words.clone();
        assert!(!compact(&mut words));
        assert_eq!(words, snapshot);
    }

    #[test]
    fn polymorphic_compression_uses_each_groups_support_and_global_sample_count() {
        for samples in [32, 64, 128, 256, 512, 768] {
            let minimum = samples / if samples < 256 { 2 } else { 4 };
            let mut matrix = AlignmentMatrix::<u64>::default();
            matrix.sample_count = samples;
            let mut enough = group(&(0..minimum as u16).collect::<Vec<_>>());
            enough[0] &= !BASE_MASK;
            let below = group(&(minimum as u16..(2 * minimum - 1) as u16).collect::<Vec<_>>());
            let SiteStorage::Growing(sites) = &mut matrix.storage else {
                unreachable!()
            };
            sites.insert(0, GrowingSite::Polymorphic(vec![enough, below]));
            let events = RefCell::new(Vec::new());
            matrix.compact_groups(&|p| events.borrow_mut().push(p.clone()));
            let SiteStorage::Growing(sites) = &matrix.storage else {
                unreachable!()
            };
            let GrowingSite::Polymorphic(groups) = &sites[&0] else {
                unreachable!()
            };
            assert_eq!(range_count(&groups[0]), 1);
            assert_eq!(range_count(&groups[1]), 0);
            assert_eq!(group_base(&groups[0]), 0);
            assert_eq!(group_base(&groups[1]), 4);
            assert_eq!(events.borrow().last().unwrap().sites_processed, 1);
            assert!(events.borrow().last().unwrap().complete);
        }
    }

    #[test]
    fn polymorphic_recompression_keeps_scheduled_support_exemption() {
        let mut ranged = group(&(0..16).collect::<Vec<_>>());
        assert!(compact(&mut ranged));
        for id in 20..25 {
            append_group(&mut ranged, id, 128, group_words(MAX_SAMPLES)).unwrap();
        }
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.sample_count = 256;
        let SiteStorage::Growing(sites) = &mut matrix.storage else {
            unreachable!()
        };
        sites.insert(
            0,
            GrowingSite::Polymorphic(vec![ranged, new_group(0, 30).unwrap()]),
        );
        matrix.compact_groups(&|_| {});
        let SiteStorage::Growing(sites) = &matrix.storage else {
            unreachable!()
        };
        let GrowingSite::Polymorphic(groups) = &sites[&0] else {
            unreachable!()
        };
        assert_eq!(range_count(&groups[0]), 2);
        assert_eq!(range(&groups[0], 1), (20, 24));
    }

    #[test]
    fn final_compression_is_strictly_above_quarter_support_for_all_group_kinds() {
        for samples in [20, 32, 128, 256, 768] {
            let quarter = samples / 4;
            let below = group(&(0..quarter as u16).collect::<Vec<_>>());
            let above = group(&(quarter as u16..(2 * quarter + 1) as u16).collect::<Vec<_>>());
            let mut low_ranged = group(&(0..16).collect::<Vec<_>>());
            assert!(compact(&mut low_ranged));
            for id in 20..25 {
                append_group(&mut low_ranged, id, 128, group_words(MAX_SAMPLES)).unwrap();
            }
            let mut matrix = AlignmentMatrix::<u64>::default();
            matrix.sample_count = samples;
            let SiteStorage::Growing(sites) = &mut matrix.storage else {
                unreachable!()
            };
            sites.insert(0, GrowingSite::Monomorphic(below.clone()));
            sites.insert(1, GrowingSite::Monomorphic(above.clone()));
            sites.insert(2, GrowingSite::Polymorphic(vec![below, above]));
            if samples >= 128 {
                sites.insert(3, GrowingSite::Monomorphic(low_ranged.clone()));
            }
            matrix.compact_final(&|_| {});
            let SiteStorage::Growing(sites) = &matrix.storage else {
                unreachable!()
            };
            let GrowingSite::Monomorphic(below) = &sites[&0] else {
                unreachable!()
            };
            let GrowingSite::Monomorphic(above) = &sites[&1] else {
                unreachable!()
            };
            assert_eq!(range_count(below), 0);
            assert_eq!(range_count(above), 1);
            let GrowingSite::Polymorphic(groups) = &sites[&2] else {
                unreachable!()
            };
            assert_eq!(range_count(&groups[0]), 0);
            assert_eq!(range_count(&groups[1]), 1);
            if samples >= 128 {
                let GrowingSite::Monomorphic(words) = &sites[&3] else {
                    unreachable!()
                };
                assert_eq!(*words, low_ranged);
            }
        }
    }

    #[test]
    fn extracted_and_imported_concrete_classifiers_agree() {
        for code in 0..=GAP {
            assert_eq!(
                is_concrete_code(code),
                concrete_base_code(decode_symbol(code)).is_some()
            );
        }
        for (symbol, code) in [
            (b'a', 0),
            (b'c', 1),
            (b'g', 2),
            (b't', 3),
            (b'u', 3),
            (b'U', 3),
        ] {
            assert_eq!(concrete_base_code(symbol), Some(code));
            assert_eq!(encode_symbol(symbol), code);
        }
        for symbol in b"RYSWKMBDHVN-ryswkmbdhvn?" {
            assert_eq!(concrete_base_code(*symbol), None);
        }
    }

    fn native_distance_fixture<IntT>(k: usize, width: usize)
    where
        IntT:
            Copy + Eq + Hash + PrimInt + for<'a> crate::ska_dict::bit_encoding::UInt<'a> + 'static,
        u128: AsPrimitive<IntT>,
    {
        use crate::merge_ska_array::MergeSkaArray;
        use crate::merge_ska_dict::MergeSkaDict;
        let samples = 128;
        let mut matrix = AlignmentMatrix::<IntT>::default();
        matrix.set_expected_samples(samples).unwrap();
        let mut names: Vec<String> = (0..samples).map(|n| format!("sample_{n}")).collect();
        let mut native_calls: HashMap<IntT, Vec<u8>> = HashMap::new();
        let mut rows = Vec::new();
        for sample in 0..samples {
            let calls = [
                Some(if sample < 32 {
                    0
                } else if sample < 64 {
                    4
                } else {
                    2
                }),
                Some(if sample < 32 {
                    4
                } else if sample < 64 {
                    0
                } else {
                    3
                }),
                Some((4 + sample % 12) as u8),
                Some((sample % 16) as u8),
                (sample % 5 != 0).then_some((sample % 4) as u8),
                Some(if sample % 2 == 0 { 0 } else { 3 }),
            ];
            let mut packet = Vec::new();
            for (key, code) in calls.into_iter().enumerate() {
                if let Some(code) = code {
                    let key = key as u128;
                    native_calls
                        .entry(key.as_())
                        .or_insert_with(|| vec![b'-'; samples])[sample] = decode_symbol(code);
                    if width == 2 {
                        packet.extend(pack_entry_u64(key as u64, decode_symbol(code)));
                    } else {
                        packet.extend(pack_entry_u128(key, decode_symbol(code)));
                    }
                }
            }
            let mut row = vec![0; sample];
            matrix.add_sample_words(packet, &mut row, &|_| {}).unwrap();
            rows.push(row);
        }
        assert!(contains_ranges(&matrix));
        let mut dynamic = MergeSkaDict::<IntT>::new(k, samples, true);
        dynamic.build_from_array(&mut names, &mut native_calls);
        let array = MergeSkaArray::new(&dynamic);
        let native = array.distance(0.0, true);
        for (row, distances) in rows.iter().enumerate() {
            for (column, distance) in distances.iter().enumerate() {
                assert_eq!(
                    f64::from(*distance),
                    native[column][row - column - 1].distance()
                );
            }
        }
    }

    #[test]
    fn u64_distances_match_native_default_with_all_iupac_calls() {
        native_distance_fixture::<u64>(31, 2);
    }

    #[test]
    fn u128_distances_match_native_default_with_all_iupac_calls() {
        native_distance_fixture::<u128>(33, 4);
    }

    fn oracle<IntT>(samples: usize, width: usize)
    where
        IntT: Copy + Eq + Hash + PrimInt + 'static,
        u128: AsPrimitive<IntT>,
    {
        let mut matrix = AlignmentMatrix::<IntT>::default();
        matrix.set_expected_samples(255).unwrap();
        let mut plain = AlignmentMatrix::<IntT>::default();
        plain.set_expected_samples(samples).unwrap();
        let events = RefCell::new(Vec::new());
        let mut history: Vec<[Option<u8>; 8]> = Vec::new();
        for sample in 0..samples {
            if sample == 255 || sample == 511 || sample == 767 || sample == 1023 {
                matrix
                    .set_expected_samples(samples.min(sample + 256))
                    .unwrap();
            }
            let mut calls = [None; 8];
            calls[0] = Some(if sample < 1024 { 0 } else { 4 });
            calls[1] = (sample % 7 != 3).then_some(1);
            calls[2] = Some((sample % 3) as u8);
            calls[3] = Some((sample % 15) as u8);
            calls[4] = (sample % 12 < 9).then_some(3);
            calls[5] = (sample % 2 == 0).then_some(2);
            calls[6] = (sample.wrapping_mul(1103515245).wrapping_add(12345) % 17 < 10)
                .then_some((sample % 4) as u8);
            calls[7] = (sample > 700 && sample % 5 < 3).then_some(14);
            let mut packet = Vec::new();
            for (key, base) in calls.iter().enumerate() {
                if let Some(base) = base {
                    let key = key as u128 + if width == 4 { 1u128 << 80 } else { 0 };
                    if width == 2 {
                        packet.extend(pack_entry_u64(key as u64, decode_symbol(*base)));
                    } else {
                        packet.extend(pack_entry_u128(key, decode_symbol(*base)));
                    }
                }
            }
            let mut row = vec![7; sample];
            let mut plain_row = row.clone();
            plain.store_and_compare(&packet, &mut plain_row).unwrap();
            matrix
                .add_sample_words(packet, &mut row, &|p| {
                    if p.complete {
                        events.borrow_mut().push(p.clone());
                    }
                })
                .unwrap();
            let expected: Vec<u32> = history
                .iter()
                .map(|earlier| {
                    7 + earlier
                        .iter()
                        .zip(calls)
                        .filter(|(a, b)| matches!((a, b), (Some(a), Some(b)) if *a < 4 && *b < 4 && *a != *b))
                        .count() as u32
                })
                .collect();
            assert_eq!(row, expected, "sample {sample}");
            assert_eq!(plain_row, expected);
            history.push(calls);
        }
        let expected_checkpoints: Vec<usize> =
            (1..=samples).filter(|n| should_compact(*n)).collect();
        assert_eq!(
            events
                .borrow()
                .iter()
                .map(|p| p.processed_samples)
                .collect::<Vec<_>>(),
            expected_checkpoints
        );
        assert!(contains_ranges(&matrix));
        assert!(!matrix.is_sealed());
        let capacity = packed_capacity(&matrix);
        let mut row = vec![0; samples];
        if samples == MAX_SAMPLES {
            assert!(matrix
                .add_sample_words(Vec::new(), &mut row, &|_| {})
                .is_err());
            assert_eq!(matrix.sample_count, samples);
            assert!(!matrix.is_sealed());
        }
        matrix.freeze().unwrap();
        assert_eq!(packed_capacity(&matrix), capacity);
        let names: Vec<String> = (0..samples).map(|n| format!("sample_{n}")).collect();
        struct FailWriter {
            remaining: usize,
        }
        impl Write for FailWriter {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                if self.remaining == 0 {
                    return Err(io::Error::other("injected"));
                }
                let written = bytes.len().min(self.remaining);
                self.remaining -= written;
                Ok(written)
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        assert!(matrix
            .write_fasta(&names, &mut FailWriter { remaining: 32 }, &|_| {})
            .is_err());
        let mut fasta = Vec::new();
        matrix.write_fasta(&names, &mut fasta, &|_| {}).unwrap();
        let mut columns = vec![Vec::new(); 8];
        let text = String::from_utf8(fasta).unwrap();
        for (sample, sequence) in text
            .lines()
            .filter(|line| !line.starts_with('>'))
            .enumerate()
        {
            for (column, base) in sequence.bytes().enumerate() {
                columns[column].push(base);
            }
            assert_eq!(sequence.len(), 8);
            assert!(sample < samples);
        }
        let mut expected_columns: Vec<Vec<u8>> = (0..8)
            .map(|key| {
                history
                    .iter()
                    .map(|calls| calls[key].map_or(b'-', decode_symbol))
                    .collect()
            })
            .collect();
        columns.sort();
        expected_columns.sort();
        assert_eq!(columns, expected_columns);
        matrix.export_fasta_gzip(&names, &|_| {}).unwrap();
        assert!(matches!(matrix.storage, SiteStorage::Released));
        assert_eq!(packed_capacity(&matrix), 0);
    }

    #[test]
    fn every_u64_distance_row_matches_independent_oracle() {
        oracle::<u64>(2048, 2);
    }

    #[test]
    fn every_u128_distance_row_matches_independent_oracle() {
        oracle::<u128>(1025, 4);
    }
}
