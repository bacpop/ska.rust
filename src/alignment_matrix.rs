//! Sparse, sample-ID/base packed alignment storage for the WASM workflow.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::hash::Hash;
use std::io::{self, Write};

use hashbrown::hash_map::Entry;
use hashbrown::HashMap;
use num_traits::{AsPrimitive, PrimInt};

const GAP: u8 = 15;

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

const SAMPLE_ID_BITS: usize = 11;
const SAMPLE_ID_MASK: u16 = 0x07ff;
const BASE_MASK: u16 = 0x000f;

pub(crate) const MAX_SAMPLES: usize = 1 << SAMPLE_ID_BITS;

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

/// Number of words needed for a base/count header and `samples` packed 11-bit IDs.
#[inline]
fn group_words(samples: usize) -> usize {
    1 + (SAMPLE_ID_BITS * samples).div_ceil(16)
}

#[inline]
fn group_base(words: &[u16]) -> u8 {
    (words[0] & BASE_MASK) as u8
}

#[inline]
fn group_len(words: &[u16]) -> usize {
    usize::from(words[0] >> 4)
}

/// Read one sample ID from the bit stream following the group header.
#[inline]
fn group_id(words: &[u16], index: usize) -> u16 {
    debug_assert!(index < group_len(words));

    let bit = SAMPLE_ID_BITS * index;
    let word = 1 + bit / 16;
    let shift = bit % 16;

    let mut value = u32::from(words[word]);
    if shift + SAMPLE_ID_BITS > 16 {
        value |= u32::from(words[word + 1]) << 16;
    }

    ((value >> shift) as u16) & SAMPLE_ID_MASK
}

/// Reserve the run's final capacity when a group reaches the word threshold.
#[cold]
fn reserve_final_group(words: &mut Vec<u16>, final_words: usize) -> Result<(), String> {
    words
        .try_reserve_exact(final_words - words.len())
        .map_err(|error| format!("Could not extend site group: {error}"))
}

#[inline]
fn append_group(
    words: &mut Vec<u16>,
    sample_id: u16,
    threshold_words: usize,
    final_words: usize,
) -> Result<(), String> {
    if words.len() == threshold_words {
        reserve_final_group(words, final_words)?;
    }
    let count = group_len(words);
    let bit = SAMPLE_ID_BITS * count;
    let word = 1 + bit / 16;
    let shift = bit % 16;
    if shift == 0 || shift + SAMPLE_ID_BITS > 16 {
        words.push(0);
    }
    let value = u32::from(sample_id) << shift;
    words[word] |= value as u16;
    if shift + SAMPLE_ID_BITS > 16 {
        words[word + 1] |= (value >> 16) as u16;
    }
    words[0] = (words[0] & BASE_MASK) | (((count + 1) as u16) << 4);
    Ok(())
}

#[inline]
fn new_group(base: u8, sample_id: u16) -> Result<Vec<u16>, String> {
    let mut words = Vec::new();
    words
        .try_reserve_exact(2)
        .map_err(|error| format!("Could not allocate site group: {error}"))?;
    words.extend_from_slice(&[(1u16 << 4) | u16::from(base), sample_id]);
    Ok(words)
}

#[inline]
fn accumulate_group_distances(group: &[u16], distance_row: &mut [u32]) {
    // The validated packet contains at most u32::MAX unique split-kmer entries.
    // Each earlier sample can receive at most one increment per entry.
    for index in 0..group_len(group) {
        distance_row[usize::from(group_id(group, index))] += 1;
    }
}

/// A site while samples are still being added.
#[derive(Debug, Clone)]
enum GrowingSite {
    Monomorphic(Vec<u16>),
    Polymorphic(Vec<Vec<u16>>),
}

#[derive(Debug, Clone)]
struct FrozenGroup {
    words: Vec<u16>,
    cursor: u16,
}

impl FrozenGroup {
    fn take_sample(&mut self, sample_id: u16) -> Option<u8> {
        let count = group_len(&self.words);
        if usize::from(self.cursor) >= count
            || group_id(&self.words, usize::from(self.cursor)) != sample_id
        {
            return None;
        }

        self.cursor += 1;
        Some(group_base(&self.words))
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
        let final_words = group_words(total);
        let threshold_words = 1usize << (final_words - 1).ilog2();
        if self.sample_count != 0 {
            let SiteStorage::Growing(postings) = &mut self.storage else {
                return Err("The alignment is closed to new samples.".into());
            };
            let prepare = |words: &mut Vec<u16>| -> Result<(), String> {
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
                            append_group(existing, sample_id, self.growth_threshold_words, self.final_group_words)?;
                        }
                        GrowingSite::Monomorphic(existing) => {
                            let new_words = new_group(code, sample_id)?;

                            let mut groups = Vec::new(); // This should be equivalent to with_capacity
                            groups.try_reserve_exact(2).map_err(|error| {
                                format!("Could not allocate site groups: {error}")
                            })?;

                            // Get the distances
                            accumulate_group_distances(existing, distance_row);

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
                                    append_group(group, sample_id, self.growth_threshold_words, self.final_group_words)?;
                                    found = true;
                                } else {
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
                        .push(FrozenSite::Monomorphic(FrozenGroup { words, cursor: 0 }));
                }
                GrowingSite::Polymorphic(groups) => {
                    let start = frozen.polymorphic_groups.len();
                    let len = groups.len();
                    for words in groups {
                        frozen
                            .polymorphic_groups
                            .push(FrozenGroup { words, cursor: 0 });
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
                FrozenSite::Monomorphic(group) => group.cursor = 0,
                FrozenSite::Polymorphic { start, len } => {
                    for group in &mut frozen.polymorphic_groups[*start..(*start + *len)] {
                        group.cursor = 0;
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
