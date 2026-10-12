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

// Range creation threshold
const MIN_RANGE_SAMPLES: usize = 3; // Min. number of consecutive samples we should have to create a range from them

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

/// Move an overlapping packed span right, preserving every bit outside it.
/// Read each chunk before overwriting it, and work backwards through the span.
fn shift_packed_bits_right(words: &mut [u16], start: usize, bits: usize, shift: usize) {
    debug_assert!(start + bits + shift <= words.len() * 16);
    let mut remaining = bits;
    while remaining != 0 {
        let width = remaining.min(16);
        remaining -= width;
        let source = start + remaining;
        let source_word = source / 16;
        let source_shift = source % 16;
        let mask = (1u32 << width) - 1;
        let mut value = u32::from(words[source_word]);
        if source_shift + width > 16 {
            value |= u32::from(words[source_word + 1]) << 16;
        }
        value = (value >> source_shift) & mask;
        let destination = source + shift;
        let destination_word = destination / 16;
        let destination_shift = destination % 16;
        let shifted_mask = mask << destination_shift;
        let shifted_value = value << destination_shift;
        words[destination_word] =
            (words[destination_word] & !(shifted_mask as u16)) | shifted_value as u16;
        if destination_shift + width > 16 {
            words[destination_word + 1] = (words[destination_word + 1]
                & !((shifted_mask >> 16) as u16))
                | (shifted_value >> 16) as u16;
        }
    }
}

/// Add a sample_id to an existing vector
#[inline]
fn append_group(
    words: &mut Vec<u16>,
    sample_id: u16,
    threshold_words: usize,
    final_words: usize,
) -> Result<(), String> {
    append_group_with_reserve(
        words,
        sample_id,
        threshold_words,
        final_words,
        |words, target| {
            words
                .try_reserve_exact(target - words.len())
                .map_err(|error| format!("Could not extend site group: {error}"))
        },
    )
}

#[inline]
fn append_group_with_reserve<F>(
    words: &mut Vec<u16>,
    sample_id: u16,
    threshold_words: usize,
    final_words: usize,
    reserve: F,
) -> Result<(), String>
where
    F: FnOnce(&mut Vec<u16>, usize) -> Result<(), String>,
{
    let count = group_len(words);
    if count >= MAX_SAMPLES || usize::from(sample_id) >= MAX_SAMPLES {
        return Err("Site group exceeds the 2,048-sample limit.".into());
    }
    let ranges = range_count(words);
    // Can we extend the last range?
    let end;
    if words[0] & EXTEND_MASK != 0 {
        debug_assert!(ranges > 0);
        let end_bit = PAYLOAD_BIT + (2 * ranges - 1) * SAMPLE_ID_BITS;
        if read_id(words, end_bit) + 1 == sample_id {
            write_id(words, end_bit, sample_id);
            words[0] = (words[0] & (BASE_MASK | EXTEND_MASK)) | ((count as u16) << 5);
            return Ok(());
        }
        // A failed extension falls through to the individual append below.
        end = payload_end(words);
    } else {
        end = payload_end(words);
        let singles = individual_count(words);
        let promoted_ids = MIN_RANGE_SAMPLES - 1;
        if singles >= promoted_ids {
            let previous = read_id(words, end - SAMPLE_ID_BITS);
            let start = read_id(words, end - promoted_ids * SAMPLE_ID_BITS);
            if start + 1 == previous && previous + 1 == sample_id {
                let range_bit = individual_start(words);
                // Two individual fields become two endpoints; no growth is needed.
                shift_packed_bits_right(
                    words,
                    range_bit,
                    (singles - promoted_ids) * SAMPLE_ID_BITS,
                    2 * SAMPLE_ID_BITS,
                );
                write_id(words, range_bit, start);
                write_id(words, range_bit + SAMPLE_ID_BITS, sample_id);
                debug_assert!(ranges + 1 <= usize::from(RANGE_COUNT_MASK));
                words[1] = (words[1] & !RANGE_COUNT_MASK) | (ranges as u16 + 1);
                words[0] = (words[0] & BASE_MASK) | EXTEND_MASK | ((count as u16) << 5);
                return Ok(());
            }
        }
    }
    // We have not been able to create a new range, nor extend the last one. Let's add this as an individual sample_id
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
        reserve(words, target)?;
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
        }
        result
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
    fn ranges_require_three_samples_and_merge_with_individuals() {
        let ids: Vec<u16> = [0, 1]
            .into_iter()
            .chain(4..90)
            .chain([94, 95, 99])
            .chain(110..200)
            .chain([205, 207])
            .collect();
        let words = group(&ids);
        assert_eq!(range_count(&words), 2);
        assert_eq!(range(&words, 0), (4, 89));
        assert_eq!(range(&words, 1), (110, 199));
        assert_eq!(individual_count(&words), 7);
        assert_eq!(GroupIds::new(&words).collect::<Vec<_>>(), ids);
        assert_eq!(words[0] & EXTEND_MASK, 0);

        let two = group(&[0, 1]);
        let three = group(&[0, 1, 2]);
        assert_eq!(range_count(&two), 0);
        assert_eq!(range(&three, 0), (0, 2));
    }

    #[test]
    fn distances_visit_all_range_members_and_singletons_once() {
        let ids: Vec<u16> = (0..90).chain([94, 95, 99]).chain(110..200).collect();
        let words = group(&ids);
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
        let words = group(&ids);
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
    fn ranged_groups_grow_with_at_most_32_spare_words() {
        let mut words = group(&(0..100).collect::<Vec<_>>());
        for id in (120..2048).step_by(2) {
            append_group(&mut words, id, 128, group_words(2048)).unwrap();
            assert!(words.capacity() <= words.len() + 32);
        }
        assert_eq!(range(&words, 0), (0, 99));
        assert_eq!(individual_count(&words), 964);
    }

    fn check_membership(words: &[u16], ids: &[u16]) {
        assert_eq!(GroupIds::new(words).collect::<Vec<_>>(), ids);
        assert_eq!(group_len(words), ids.len());
        assert_eq!(group_base(words), 4);
        for index in 0..range_count(words) {
            let (start, end) = range(words, index);
            assert!(end - start + 1 >= MIN_RANGE_SAMPLES as u16);
        }
        if words[0] & EXTEND_MASK != 0 {
            assert_eq!(range(words, range_count(words) - 1).1, *ids.last().unwrap());
        }
        let end = payload_end(words);
        assert_eq!(
            end,
            PAYLOAD_BIT + SAMPLE_ID_BITS * (2 * range_count(words) + individual_count(words))
        );
        let padding = words.len() * 16 - end;
        assert!(padding < 16);
        if padding != 0 {
            assert_eq!(words.last().unwrap() >> (16 - padding), 0);
        }
    }

    #[test]
    fn packed_shift_preserves_bits_at_every_offset_and_overlap() {
        for offset in 0..32 {
            for bits in 0..100 {
                for shift in [1, 6, 11, 16, 22, 31] {
                    let mut words: Vec<u16> = (0..12)
                        .map(|index| (index as u16).wrapping_mul(19873) ^ 0xa56c)
                        .collect();
                    let original = words.clone();
                    let mut expected = original.clone();
                    for bit in 0..bits {
                        let source = offset + bit;
                        let destination = source + shift;
                        let value = (original[source / 16] >> (source % 16)) & 1;
                        expected[destination / 16] = (expected[destination / 16]
                            & !(1 << (destination % 16)))
                            | (value << (destination % 16));
                    }
                    shift_packed_bits_right(&mut words, offset, bits, shift);
                    assert_eq!(
                        words, expected,
                        "offset={offset}, bits={bits}, shift={shift}"
                    );
                }
            }
        }
    }

    #[test]
    fn third_id_creates_range_and_extension_never_reserves() {
        let mut words = group(&[0, 1]).into_boxed_slice().into_vec();
        assert_eq!(range_count(&words), 0);
        let length = words.len();
        let capacity = words.capacity();
        let end = payload_end(&words);
        for id in 2..=2047 {
            append_group_with_reserve(&mut words, id, 128, group_words(MAX_SAMPLES), |_, _| {
                panic!("Range update must not reserve")
            })
            .unwrap();
            assert_eq!(range(&words, 0), (0, id));
            assert_eq!(group_len(&words), usize::from(id) + 1);
            assert_eq!(words.len(), length);
            assert_eq!(words.capacity(), capacity);
            assert_eq!(payload_end(&words), end);
        }
        check_membership(&words, &(0..2048).collect::<Vec<_>>());
        let snapshot = words.clone();
        assert!(append_group(&mut words, 2047, 128, group_words(MAX_SAMPLES)).is_err());
        assert_eq!(words, snapshot);
        assert!(new_group(0, 2048).is_err());
    }

    #[test]
    fn sparse_individuals_cover_padding_without_creating_ranges() {
        let mut words = new_group(4, 0).unwrap();
        let mut paddings = [false; 16];
        let mut ids = vec![0];
        for id in (2..2048).step_by(2) {
            append_group(&mut words, id, 128, group_words(MAX_SAMPLES)).unwrap();
            ids.push(id);
            assert_eq!(range_count(&words), 0);
            paddings[usize::from((words[1] & PADDING_MASK) >> 10)] = true;
            check_membership(&words, &ids);
        }
        assert!(paddings.iter().all(|seen| *seen));
    }

    #[test]
    fn promotion_shifts_earlier_individuals_and_preserves_ranges() {
        let mut words = group(&[0, 1, 3, 4, 6, 7, 9, 10]);
        let length = words.len();
        let capacity = words.capacity();
        append_group(&mut words, 11, 128, group_words(MAX_SAMPLES)).unwrap();
        assert_eq!(range(&words, 0), (9, 11));
        assert_eq!(individual_count(&words), 6);
        assert_eq!(words.len(), length);
        assert_eq!(words.capacity(), capacity);
        for id in [15, 16, 17, 18, 20, 21, 25, 26, 27] {
            append_group(&mut words, id, 128, group_words(MAX_SAMPLES)).unwrap();
        }
        assert_eq!(range(&words, 0), (9, 11));
        assert_eq!(range(&words, 1), (15, 18));
        assert_eq!(range(&words, 2), (25, 27));
        check_membership(
            &words,
            &[
                0, 1, 3, 4, 6, 7, 9, 10, 11, 15, 16, 17, 18, 20, 21, 25, 26, 27,
            ],
        );
    }

    #[test]
    fn dynamic_groups_match_plain_membership_after_every_append() {
        for seed in 0..16 {
            let mut ids = Vec::new();
            let mut words = None;
            for id in 0..2048u16 {
                let position = (usize::from(id) + seed) % 17;
                if position < 2
                    || (4..8).contains(&position)
                    || (9..12).contains(&position)
                    || position >= 13
                {
                    if let Some(words) = words.as_mut() {
                        append_group(words, id, 128, group_words(MAX_SAMPLES)).unwrap();
                    } else {
                        words = Some(new_group(4, id).unwrap());
                    }
                    ids.push(id);
                    check_membership(words.as_ref().unwrap(), &ids);
                }
            }
        }
    }

    #[test]
    fn failed_individual_reservation_preserves_group_and_capacity() {
        let mut words = group(&[0, 2, 4]).into_boxed_slice().into_vec();
        let snapshot = words.clone();
        let capacity = words.capacity();
        let failure =
            append_group_with_reserve(&mut words, 6, 128, group_words(MAX_SAMPLES), |_, _| {
                Err("Injected reservation failure".into())
            });
        assert!(failure.is_err());
        assert_eq!(words, snapshot);
        assert_eq!(words.capacity(), capacity);
        assert!(append_group(&mut words, 2048, 128, group_words(MAX_SAMPLES)).is_err());
        assert_eq!(words, snapshot);
        append_group(&mut words, 6, 128, group_words(MAX_SAMPLES)).unwrap();
        check_membership(&words, &[0, 2, 4, 6]);
    }

    #[test]
    fn polymorphic_groups_form_ranges_without_support_thresholds() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.set_expected_samples(20).unwrap();
        for sample in 0..20 {
            let base = if sample < 3 {
                b'A'
            } else if sample < 6 {
                b'G'
            } else if sample < 9 {
                b'R'
            } else {
                b'T'
            };
            let mut row = vec![0; sample];
            matrix
                .add_sample_words(pack_entry_u64(0, base).to_vec(), &mut row)
                .unwrap();
            for (earlier, distance) in row.into_iter().enumerate() {
                let old = if earlier < 3 {
                    b'A'
                } else if earlier < 6 {
                    b'G'
                } else if earlier < 9 {
                    b'R'
                } else {
                    b'T'
                };
                assert_eq!(
                    distance,
                    u32::from(base != old && base != b'R' && old != b'R')
                );
            }
        }
        let SiteStorage::Growing(sites) = &matrix.storage else {
            unreachable!()
        };
        let GrowingSite::Polymorphic(groups) = &sites[&0] else {
            unreachable!()
        };
        assert_eq!(groups.len(), 4);
        for (index, words) in groups.iter().enumerate() {
            assert_eq!(range_count(words), 1);
            assert_eq!(
                range(words, 0),
                if index < 3 {
                    (index as u16 * 3, index as u16 * 3 + 2)
                } else {
                    (9, 19)
                }
            );
        }
    }

    #[test]
    fn batch_preparation_preserves_dynamic_ranges_and_payloads() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.set_expected_samples(3).unwrap();
        for sample in 0..3 {
            matrix
                .add_sample_words(pack_entry_u64(0, b'A').to_vec(), &mut vec![0; sample])
                .unwrap();
        }
        let before = packed_capacity(&matrix);
        matrix.set_expected_samples(1000).unwrap();
        assert_eq!(packed_capacity(&matrix), before);
        for sample in 3..1000 {
            matrix
                .add_sample_words(pack_entry_u64(0, b'A').to_vec(), &mut vec![0; sample])
                .unwrap();
        }
        let SiteStorage::Growing(sites) = &matrix.storage else {
            unreachable!()
        };
        let GrowingSite::Monomorphic(words) = &sites[&0] else {
            unreachable!()
        };
        assert_eq!(range(&words, 0), (0, 999));
        assert_eq!(packed_capacity(&matrix), before);
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
            matrix.add_sample_words(packet, &mut row).unwrap();
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
            matrix.add_sample_words(packet, &mut row).unwrap();
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
        assert!(contains_ranges(&matrix));
        assert!(!matrix.is_sealed());
        let capacity = packed_capacity(&matrix);
        let mut row = vec![0; samples];
        if samples == MAX_SAMPLES {
            assert!(matrix.add_sample_words(Vec::new(), &mut row).is_err());
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
