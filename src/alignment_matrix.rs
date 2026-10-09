//! Sparse, sample-ID/base packed alignment storage for the WASM workflow.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::hash::Hash;
use std::io::{self, Write};

use hashbrown::hash_map::Entry;
use hashbrown::HashMap;
use num_traits::PrimInt;

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

fn unpack_entry<IntT>(words: &[u32]) -> Result<(IntT, u8), String>
where
    IntT: PrimInt + num_traits::NumCast,
{
    let mut packed = 0u128;
    for (index, word) in words.iter().enumerate() {
        packed |= u128::from(*word) << (index * 32);
    }
    let shift = words.len() * 32 - 4;
    let code = ((packed >> shift) & 0xF) as u8;
    let key = num_traits::NumCast::from(packed & ((1u128 << shift) - 1))
        .ok_or_else(|| "Packed key does not fit its integer type.".to_string())?;
    Ok((key, code))
}

fn validate_words(words: &[u32], width: usize) -> Result<(), String> {
    if width != 2 && width != 4 {
        return Err("Packed entries must contain two or four words.".into());
    }
    if !words.len().is_multiple_of(width) {
        return Err("The packed sample contains an incomplete entry.".into());
    }
    Ok(())
}

/// Number of words needed for a base/count header and `samples` packed 11-bit IDs.
fn group_words(samples: usize) -> usize {
    1 + (SAMPLE_ID_BITS * samples).div_ceil(16)
}

fn group_base(words: &[u16]) -> u8 {
    (words[0] & BASE_MASK) as u8
}

fn group_len(words: &[u16]) -> usize {
    usize::from(words[0] >> 4)
}

/// Read one sample ID from the bit stream following the group header.
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

/// Reserve enough packed words for one more sample ID, growing geometrically up to the planned
/// sample count. An underestimated plan still grows to the required length.
fn reserve_group_append(words: &mut Vec<u16>, expected_samples: usize) -> Result<(), String> {
    let needed = group_words(group_len(words) + 1);
    if needed <= words.capacity() {
        return Ok(());
    }
    let target = words
        .capacity()
        .saturating_mul(2)
        .min(group_words(expected_samples))
        .max(needed);
    words
        .try_reserve_exact(target - words.len())
        .map_err(|error| format!("Could not extend site group: {error}"))
}

/// Append an ID after `reserve_group_append` has made the required capacity available.
fn append_group_reserved(words: &mut Vec<u16>, sample_id: u16) {
    let count = group_len(words);
    let bit = SAMPLE_ID_BITS * count;
    let word = 1 + bit / 16;
    let shift = bit % 16;

    words.resize(group_words(count + 1), 0);

    let value = u32::from(sample_id & SAMPLE_ID_MASK) << shift;
    words[word] |= value as u16;
    if shift + SAMPLE_ID_BITS > 16 {
        words[word + 1] |= (value >> 16) as u16;
    }

    words[0] = (words[0] & BASE_MASK) | (((count + 1) as u16) << 4);
}

fn new_group(base: u8, sample_id: u16) -> Result<Vec<u16>, String> {
    let mut words = Vec::new();
    words
        .try_reserve_exact(group_words(1))
        .map_err(|error| format!("Could not allocate site group: {error}"))?;
    words.push(u16::from(base & BASE_MASK as u8));
    // The new vector has capacity for its header and first ID.
    append_group_reserved(&mut words, sample_id);
    Ok(words)
}

fn increment_distance(distance: &mut u32) -> Result<(), String> {
    *distance = distance
        .checked_add(1)
        .ok_or_else(|| "SNP distance exceeds u32::MAX.".to_string())?;
    Ok(())
}

fn add_group_distances(
    group: &[u16],
    incoming_base: u8,
    distance_row: &mut [u32],
) -> Result<(), String> {
    let previous_base = group_base(group);
    if previous_base == GAP || previous_base == incoming_base {
        return Ok(());
    }

    for index in 0..group_len(group) {
        let sample_id = usize::from(group_id(group, index));
        increment_distance(&mut distance_row[sample_id])?;
    }
    Ok(())
}

/// A site while samples are still being added.
#[derive(Debug, Clone)]
enum GrowingSite {
    Monomorphic(Vec<u16>),
    Polymorphic(Vec<Vec<u16>>),
}

impl GrowingSite {
    #[cfg(test)]
    fn groups(&self) -> &[Vec<u16>] {
        match self {
            Self::Monomorphic(words) => std::slice::from_ref(words),
            Self::Polymorphic(groups) => groups,
        }
    }
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
    /// Samples the dataset will hold once the planned run is added; caps each site's capacity.
    expected_samples: usize,
}

impl<IntT> Default for AlignmentMatrix<IntT> {
    fn default() -> Self {
        Self {
            storage: SiteStorage::Growing(HashMap::new()),
            sample_count: 0,
            sealed: false,
            expected_samples: MAX_SAMPLES,
        }
    }
}

impl<IntT> AlignmentMatrix<IntT>
where
    IntT: Copy + Eq + Hash + PrimInt + num_traits::NumCast,
{
    pub(crate) fn is_sealed(&self) -> bool {
        self.sealed
    }

    /// Set the total the dataset will reach after the planned samples, at most [`MAX_SAMPLES`].
    pub(crate) fn set_expected_samples(&mut self, total: usize) {
        self.expected_samples = total.min(MAX_SAMPLES);
    }

    fn check_addition(&self, words: &[u32], width: usize) -> Result<(), String> {
        if self.sealed {
            return Err("The alignment is closed to new samples.".into());
        }
        if self.sample_count >= MAX_SAMPLES {
            return Err(format!(
                "Sparse alignment supports at most {MAX_SAMPLES} samples."
            ));
        }
        validate_words(words, width)
    }

    /// Store one sample and write its distances to earlier samples into `distance_row`,
    /// visiting each split k-mer once. A failure after storing has begun seals the matrix.
    pub(crate) fn add_sample_words(
        &mut self,
        words: Vec<u32>,
        words_per_entry: usize,
        distance_row: &mut [u32],
    ) -> Result<(), String> {
        self.check_addition(&words, words_per_entry)?;
        if distance_row.len() != self.sample_count {
            return Err("Distance row length does not match the existing sample count.".into());
        }
        let result = self.store_and_compare(&words, words_per_entry, distance_row);
        if result.is_err() {
            self.sealed = true;
        }
        result
    }

    fn store_and_compare(
        &mut self,
        words: &[u32],
        words_per_entry: usize,
        distance_row: &mut [u32],
    ) -> Result<(), String> {
        let SiteStorage::Growing(postings) = &mut self.storage else {
            return Err("The alignment is closed to new samples.".into());
        };
        let sample_id = self.sample_count as u16;
        let limit = self.expected_samples;
        distance_row.fill(0);
        for entry in words.chunks_exact(words_per_entry) {
            let (kmer, code) = unpack_entry::<IntT>(entry)?;
            // Reserving first keeps the insertion below fallible; with room to spare it is only a check.
            postings
                .try_reserve(1)
                .map_err(|error| format!("Could not extend site index: {error:?}"))?;
            match postings.entry(kmer) {
                Entry::Occupied(mut entry) => {
                    let site = entry.get_mut();
                    match site {
                        GrowingSite::Monomorphic(existing) if group_base(existing) == code => {
                            reserve_group_append(existing, limit)?;
                            append_group_reserved(existing, sample_id);
                        }
                        GrowingSite::Monomorphic(existing) => {
                            let new_words = new_group(code, sample_id)?;
                            let mut groups = Vec::new();
                            groups.try_reserve_exact(2).map_err(|error| {
                                format!("Could not allocate site groups: {error}")
                            })?;

                            add_group_distances(existing, code, distance_row)?;

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
                            let matching_group =
                                groups.iter().position(|group| group_base(group) == code);

                            if let Some(index) = matching_group {
                                reserve_group_append(&mut groups[index], limit)?;
                            } else {
                                groups.try_reserve(1).map_err(|error| {
                                    format!("Could not extend site groups: {error}")
                                })?;
                            }
                            let new_words = if matching_group.is_none() {
                                Some(new_group(code, sample_id)?)
                            } else {
                                None
                            };

                            for group in groups.iter() {
                                add_group_distances(group, code, distance_row)?;
                            }

                            if let Some(index) = matching_group {
                                append_group_reserved(&mut groups[index], sample_id);
                            } else {
                                groups.push(new_words.expect("new base group was prepared"));
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

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::GzDecoder;
    use std::io::Read;

    fn words(sample: &HashMap<u64, u8>) -> Vec<u32> {
        sample
            .iter()
            .flat_map(|(&key, &base)| pack_entry_u64(key, base))
            .collect()
    }

    fn add_words<IntT>(
        matrix: &mut AlignmentMatrix<IntT>,
        words: Vec<u32>,
        words_per_entry: usize,
        distances: &mut Vec<u32>,
    ) -> std::ops::Range<usize>
    where
        IntT: Copy + Eq + Hash + PrimInt + num_traits::NumCast,
    {
        let start = distances.len();
        let end = start + matrix.sample_count;
        distances.resize(end, 0);
        matrix
            .add_sample_words(words, words_per_entry, &mut distances[start..end])
            .unwrap();
        start..end
    }

    #[test]
    fn strict_triangle_helpers_cover_empty_and_all_pair_positions() {
        assert_eq!(strict_triangle_len(0), Some(0));
        assert_eq!(strict_triangle_len(1), Some(0));
        assert_eq!(strict_triangle_len(4), Some(6));
        assert_eq!(strict_triangle_index(1, 0), Some(0));
        assert_eq!(strict_triangle_index(2, 0), Some(1));
        assert_eq!(strict_triangle_index(2, 1), Some(2));
        assert_eq!(strict_triangle_index(3, 0), Some(3));
        assert_eq!(strict_triangle_index(3, 1), Some(4));
        assert_eq!(strict_triangle_index(3, 2), Some(5));
        assert_eq!(strict_triangle_index(2, 2), None);
        assert_eq!(strict_triangle_index(1, 2), None);
    }

    #[test]
    fn sparse_distance_rows_match_map_intersections() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut prior = Vec::<HashMap<u64, u8>>::new();
        let mut distances = Vec::new();
        for sample_id in 0..24 {
            let sample: HashMap<_, _> = (0..90u64)
                .filter(|key| (key * 17 + sample_id * 13) % 7 < 3)
                .map(|key| (key, decode_symbol(((key + sample_id) % 16) as u8)))
                .collect();
            let packed = words(&sample);
            let row = add_words(&mut matrix, packed, 2, &mut distances);
            for (j, other) in prior.iter().enumerate() {
                let expected = sample
                    .iter()
                    .filter(|(key, base)| {
                        other
                            .get(*key)
                            .is_some_and(|old| *old != b'-' && old != *base)
                    })
                    .count() as u32;
                assert_eq!(distances[row.start + j], expected);
            }
            prior.push(sample);
        }
    }

    #[test]
    fn invariant_sites_skip_comparisons_without_changing_rows() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut prior = Vec::<HashMap<u64, u8>>::new();
        let mut distances = Vec::new();
        for sample_id in 0..20u64 {
            let sample: HashMap<_, _> = (0..60u64)
                .filter(|key| (key + sample_id) % 9 != 0)
                .map(|key| {
                    let variant = key % 10 == 0 && sample_id % 3 == 0;
                    (key, if variant { b'T' } else { b'A' })
                })
                .collect();
            let row = add_words(&mut matrix, words(&sample), 2, &mut distances);
            for (j, other) in prior.iter().enumerate() {
                let expected = sample
                    .iter()
                    .filter(|(key, base)| other.get(*key).is_some_and(|old| old != *base))
                    .count() as u32;
                assert_eq!(distances[row.start + j], expected);
            }
            prior.push(sample);
        }
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        assert!(postings
            .iter()
            .all(|(key, site)| matches!(site, GrowingSite::Polymorphic(_)) == (key % 10 == 0)));
    }

    #[test]
    fn site_capacity_stops_at_the_expected_sample_count() {
        let shared = || pack_entry_u64(5, b'A').to_vec();
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut distances = Vec::new();
        matrix.set_expected_samples(6);
        for _ in 0..6 {
            add_words(&mut matrix, shared(), 2, &mut distances);
        }
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        // 1, then 4, then 6 rather than 8.
        let GrowingSite::Monomorphic(words) = &postings[&5] else {
            panic!()
        };
        assert_eq!(words.capacity(), 6);

        // An underestimate still grows, one call beyond the limit at a time.
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut distances = Vec::new();
        matrix.set_expected_samples(2);
        for _ in 0..5 {
            add_words(&mut matrix, shared(), 2, &mut distances);
        }
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        let GrowingSite::Monomorphic(words) = &postings[&5] else {
            panic!()
        };
        assert_eq!(group_len(words), 5);
    }

    #[test]
    fn zero_call_is_valid_and_sample_limit_is_checked() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut distances = Vec::new();
        add_words(
            &mut matrix,
            pack_entry_u64(9, b'A').to_vec(),
            2,
            &mut distances,
        );
        let second_row = add_words(
            &mut matrix,
            pack_entry_u64(9, b'C').to_vec(),
            2,
            &mut distances,
        );
        assert_eq!(&distances[second_row], &[1]);
        for _ in 2..MAX_SAMPLES - 1 {
            add_words(&mut matrix, Vec::new(), 2, &mut distances);
        }
        add_words(
            &mut matrix,
            pack_entry_u64(9, b'N').to_vec(),
            2,
            &mut distances,
        );
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        let GrowingSite::Polymorphic(groups) = &postings[&9] else {
            panic!()
        };
        assert_eq!(groups.len(), 3);
        assert_eq!(group_len(&groups[0]), 1);
        assert_eq!(group_id(&groups[0], 0), 0);
        assert_eq!(group_base(&groups[0]), 0);
        assert_eq!(group_id(&groups[1], 0), 1);
        assert_eq!(group_base(&groups[1]), 1);
        assert_eq!(group_id(&groups[2], 0), (MAX_SAMPLES - 1) as u16);
        assert_eq!(group_base(&groups[2]), 14);
        assert!(matrix.add_sample_words(Vec::new(), 2, &mut []).is_err());
    }

    #[test]
    fn large_keys_ambiguity_and_empty_samples_work() {
        let mut matrix = AlignmentMatrix::<u128>::default();
        let mut distances = Vec::new();
        let key = (1u128 << 100) + 27;
        add_words(
            &mut matrix,
            pack_entry_u128(key, b'R').to_vec(),
            4,
            &mut distances,
        );
        add_words(&mut matrix, Vec::new(), 4, &mut distances);
        add_words(
            &mut matrix,
            pack_entry_u128(key, b'R').to_vec(),
            4,
            &mut distances,
        );
        let fourth_row = add_words(
            &mut matrix,
            pack_entry_u128(key, b'A').to_vec(),
            4,
            &mut distances,
        );
        assert_eq!(&distances[fourth_row], &[1, 0, 1]);
        assert!(matrix
            .add_sample_words(vec![1], 4, &mut [0, 0, 0, 0])
            .is_err());
        assert!(matrix
            .add_sample_words(Vec::new(), 0, &mut [0, 0, 0, 0])
            .is_err());
        let mut distance = u32::MAX;
        assert!(increment_distance(&mut distance).is_err());
    }

    #[test]
    fn failed_addition_after_a_partial_insert_seals_the_matrix() {
        let mut matrix = AlignmentMatrix::<u8>::default();
        let mut invalid_sample = pack_entry_u64(3, b'A').to_vec();
        invalid_sample.extend(pack_entry_u64(300, b'C'));
        assert!(matrix.add_sample_words(invalid_sample, 2, &mut []).is_err());
        assert!(matrix.is_sealed());
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        assert_eq!(postings.len(), 1);
        let GrowingSite::Monomorphic(group) = &postings[&3] else {
            panic!()
        };
        assert_eq!(group_len(group), 1);
        assert!(matrix.add_sample_words(Vec::new(), 2, &mut [0, 0]).is_err());
    }

    #[test]
    fn freeze_moves_buffers_and_export_recovers_gaps_and_symbols() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut distances = Vec::new();
        let first: HashMap<_, _> = (0..16)
            .map(|code| (code as u64, decode_symbol(code)))
            .collect();
        add_words(&mut matrix, words(&first), 2, &mut distances);
        add_words(
            &mut matrix,
            pack_entry_u64(0, b'G').to_vec(),
            2,
            &mut distances,
        );
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        let mut pointers: Vec<_> = postings
            .values()
            .flat_map(|site| site.groups().iter().map(|group| group.as_ptr() as usize))
            .collect();
        pointers.sort_unstable();
        matrix.freeze().unwrap();
        let SiteStorage::Frozen(frozen) = &matrix.storage else {
            panic!()
        };
        let mut frozen_pointers = Vec::new();
        for site in &frozen.sites {
            match site {
                FrozenSite::Monomorphic(group) => {
                    frozen_pointers.push(group.words.as_ptr() as usize);
                }
                FrozenSite::Polymorphic { start, len } => frozen_pointers.extend(
                    frozen.polymorphic_groups[*start..(*start + *len)]
                        .iter()
                        .map(|group| group.words.as_ptr() as usize),
                ),
            }
        }
        frozen_pointers.sort_unstable();
        assert_eq!(frozen_pointers, pointers);
        let mut bytes = Vec::new();
        GzDecoder::new(
            matrix
                .export_fasta_gzip(&["one".into(), "two".into()], &|_| {})
                .unwrap()
                .as_slice(),
        )
        .read_to_end(&mut bytes)
        .unwrap();
        let text = String::from_utf8(bytes).unwrap();
        let lines: Vec<_> = text.lines().collect();
        assert_eq!((lines[0], lines[2]), ("> one", "> two"));
        assert_eq!(lines[1].len(), 16);
        let mut columns: Vec<_> = lines[1].bytes().zip(lines[3].bytes()).collect();
        columns.sort_unstable();
        let mut expected: Vec<_> = (0..16)
            .map(|code| (decode_symbol(code), if code == 0 { b'G' } else { b'-' }))
            .collect();
        expected.sort_unstable();
        assert_eq!(columns, expected);
        assert!(matches!(matrix.storage, SiteStorage::Released));
        assert!(matrix.add_sample_words(Vec::new(), 2, &mut []).is_err());
    }

    #[test]
    fn failed_writer_keeps_frozen_data_and_retry_resets_cursors() {
        struct FailAfter(usize);
        impl Write for FailAfter {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                if self.0 == 0 {
                    return Err(io::Error::other("simulated write failure"));
                }
                let n = self.0.min(bytes.len());
                self.0 -= n;
                Ok(n)
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut distances = Vec::new();
        add_words(
            &mut matrix,
            pack_entry_u64(1, b'A').to_vec(),
            2,
            &mut distances,
        );
        add_words(
            &mut matrix,
            pack_entry_u64(1, b'T').to_vec(),
            2,
            &mut distances,
        );
        matrix.freeze().unwrap();
        assert!(matrix
            .write_fasta(&["a".into(), "b".into()], &mut FailAfter(6), &|_| {})
            .is_err());
        assert!(matrix.is_sealed());
        let mut plain = Vec::new();
        matrix
            .write_fasta(&["a".into(), "b".into()], &mut plain, &|_| {})
            .unwrap();
        assert_eq!(plain, b"> a\nA\n> b\nT\n");
    }

    #[test]
    fn packed_group_header_and_ids_cover_the_full_sample_range() {
        let mut words = new_group(14, 0).unwrap();
        for sample_id in 1..MAX_SAMPLES as u16 {
            reserve_group_append(&mut words, MAX_SAMPLES).unwrap();
            append_group_reserved(&mut words, sample_id);
        }

        assert_eq!(group_base(&words), 14);
        assert_eq!(group_len(&words), MAX_SAMPLES);
        assert_eq!(words.len(), group_words(MAX_SAMPLES));
        for sample_id in 0..MAX_SAMPLES {
            assert_eq!(group_id(&words, sample_id), sample_id as u16);
        }
    }

    #[test]
    fn freeze_reservation_failure_preserves_growing_data_for_retry() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut distances = Vec::new();
        add_words(
            &mut matrix,
            pack_entry_u64(7, b'A').to_vec(),
            2,
            &mut distances,
        );
        add_words(
            &mut matrix,
            pack_entry_u64(7, b'C').to_vec(),
            2,
            &mut distances,
        );

        assert!(matrix
            .freeze_with(|_, _| Err("injected reservation failure".into()))
            .is_err());
        assert!(matrix.is_sealed());
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!("growing data must remain available after reservation failure")
        };
        let GrowingSite::Polymorphic(groups) = &postings[&7] else {
            panic!()
        };
        assert_eq!(groups.len(), 2);
        assert_eq!(group_len(&groups[0]), 1);
        assert_eq!(group_len(&groups[1]), 1);

        matrix.freeze().unwrap();
        assert!(matches!(matrix.storage, SiteStorage::Frozen(_)));
    }
}
