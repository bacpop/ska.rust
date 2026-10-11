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
// Profile 1/64 of entries in contiguous blocks. Batching keeps the imported
// JavaScript clock cost small relative to the Rust work being measured.
const AGGREGATION_PROFILE_SAMPLE_INTERVAL: usize = 1 << 18;
const AGGREGATION_PROFILE_BLOCK_ENTRIES: usize = 1 << 12;

pub(crate) const MAX_SAMPLES: usize = 1 << SAMPLE_ID_BITS;

/// Batched timings for sampled blocks in the hot loop that adds one sample to
/// the sparse matrix. Unpack and lookup are measured in prepasses; the real
/// block is timed inclusively, without subtracting incompatible replay timings.
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct AggregationPhaseTimings {
    total_entries: usize,
    sampled_entries: usize,
    sampled_update_entries: usize,
    sampled_lookup_hits: usize,
    sampled_unpack_ms: f64,
    sampled_hash_lookup_ms: f64,
    sampled_block_work_ms: f64,
    sampled_resize_block_ms: f64,
    sampled_resize_blocks: usize,
    fine: FineAggregationTimings,
}

/// Independent batched replays, never exclusive phases of the production loop.
#[derive(Debug, Default, Clone, Copy)]
struct FineAggregationTimings {
    sampled_entries: usize,
    skipped_blocks: usize,
    operation_entries: [usize; 4],
    operation_ms: [f64; 4],
    distance_ms: f64,
    distance_groups: usize,
    distance_ids: u64,
    preparation_ms: f64,
    diagnostic_ms: f64,
}

impl AggregationPhaseTimings {
    /// Values returned to the browser in this order: total Rust call, estimated
    /// unpack, estimated hash lookup, estimated inclusive block, sampled entries,
    /// total entries, entries in non-resize blocks, sampling interval, block
    /// size, sampled lookup hits, total resize-block time, resize-block count.
    /// The versioned suffix contains independent operation replay estimates,
    /// their workload counts, and observed diagnostic overhead. Keep this layout
    /// in sync with aggregationProfileIndex in Mapper.ts.
    #[cfg(target_family = "wasm")]
    pub(crate) fn wasm_values(self, rust_core_ms: f64) -> [f64; 28] {
        let sample_scale = if self.sampled_entries == 0 {
            0.0
        } else {
            self.total_entries as f64 / self.sampled_entries as f64
        };
        let update_scale = if self.sampled_update_entries == 0 {
            0.0
        } else {
            self.total_entries as f64 / self.sampled_update_entries as f64
        };
        let fine_scale = if self.fine.sampled_entries == 0 {
            0.0
        } else {
            self.total_entries as f64 / self.fine.sampled_entries as f64
        };
        [
            rust_core_ms,
            self.sampled_unpack_ms * sample_scale,
            self.sampled_hash_lookup_ms * sample_scale,
            self.sampled_block_work_ms * update_scale,
            self.sampled_entries as f64,
            self.total_entries as f64,
            self.sampled_update_entries as f64,
            AGGREGATION_PROFILE_SAMPLE_INTERVAL as f64,
            AGGREGATION_PROFILE_BLOCK_ENTRIES as f64,
            self.sampled_lookup_hits as f64,
            self.sampled_resize_block_ms,
            self.sampled_resize_blocks as f64,
            3.0, // Profile schema version; legacy profiles have no version field.
            self.fine.operation_ms[0] * fine_scale,
            self.fine.operation_ms[1] * fine_scale,
            self.fine.operation_ms[2] * fine_scale,
            self.fine.operation_ms[3] * fine_scale,
            self.fine.distance_ms * fine_scale,
            self.fine.operation_entries[0] as f64,
            self.fine.operation_entries[1] as f64,
            self.fine.operation_entries[2] as f64,
            self.fine.operation_entries[3] as f64,
            self.fine.distance_groups as f64,
            self.fine.distance_ids as f64,
            self.fine.preparation_ms,
            self.fine.diagnostic_ms,
            self.fine.sampled_entries as f64,
            self.fine.skipped_blocks as f64,
        ]
    }
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
    //// OLD
    // for index in 0..group_len(group) {
    //     distance_row[usize::from(group_id(group, index))] += 1;
    // }

    // Let's try to do this, perhaps we could get some SIMD optimisations?
    let n_sampleids   = group_len(group);
    let n_full_chunks = (n_sampleids / 16) * SAMPLE_ID_BITS; // this, unless we go over one full chunk, will be zero NOTE THE ORDER!
    let thedata = &group[1..];
    let remainder = &thedata[n_full_chunks..];

    for wordchunk in thedata[..n_full_chunks].chunks_exact(SAMPLE_ID_BITS) {
        for i in 0..16 {
            let bit_offset  = i * SAMPLE_ID_BITS;
            let word        = bit_offset / 16; // truncated div., i.e. the lower possible value
            let shift       = bit_offset % 16; // remainder of the 16 bits

            let mut value = (wordchunk[word] as u32) >> shift;

            if shift > 5 { // NOTE: this needs change if modifying SAMPLE_ID_BITS
                value |= (wordchunk[word + 1] as u32) << (16 - shift); // the other bits are in the next one!
            }

            // NOTE: this needs change if modifying SAMPLE_ID_BITS
            let sample_id = (value & 0x7ff) as usize; // 0x7ff = 0b11111111111, i.e. this gets the first SAMPLE_ID_BITS bits
            distance_row[sample_id] += 1;
        }
    }
    // Now, if we miss anyone
    for i in 0..(n_sampleids % 16) {
        let bit_offset  = i * SAMPLE_ID_BITS;
        let word        = bit_offset / 16;
        let shift       = bit_offset % 16;

        let mut value = remainder[word] as u32 >> shift;
        if shift > 5 {// NOTE: this needs change if modifying SAMPLE_ID_BITS
            value |= u32::from(remainder[word + 1]) << (16 - shift);
        }

        let sample_id = (value & 0x7ff) as usize;// NOTE: this needs change if modifying SAMPLE_ID_BITS
        distance_row[sample_id] += 1;
    }

}

/// A site while samples are still being added.
#[derive(Debug, Clone)]
enum GrowingSite {
    Monomorphic(Vec<u16>),
    Polymorphic(Vec<Vec<u16>>),
}

/// Shared by the production loop and the diagnostic replay. A no-op distance
/// callback lets the replay measure site mutation without distance arithmetic.
#[inline]
fn update_existing_site<F>(
    site: &mut GrowingSite,
    code: u8,
    sample_id: u16,
    threshold_words: usize,
    final_words: usize,
    mut accumulate: F,
) -> Result<(), String>
where
    F: FnMut(&[u16]),
{
    match site {
        GrowingSite::Monomorphic(existing) if group_base(existing) == code => {
            append_group(existing, sample_id, threshold_words, final_words)?;
        }
        GrowingSite::Monomorphic(existing) => {
            let new_words = new_group(code, sample_id)?;
            let mut groups = Vec::new();
            groups
                .try_reserve_exact(2)
                .map_err(|error| format!("Could not allocate site groups: {error}"))?;
            accumulate(existing);
            let old_site = std::mem::replace(site, GrowingSite::Monomorphic(Vec::new()));
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
                if group_base(group) == code {
                    append_group(group, sample_id, threshold_words, final_words)?;
                    found = true;
                } else {
                    accumulate(group);
                }
            }
            if !found {
                let new_words = new_group(code, sample_id)?;
                groups
                    .try_reserve(1)
                    .map_err(|error| format!("Could not extend site groups: {error}"))?;
                groups.push(new_words);
            }
        }
    }
    Ok(())
}

/// Preserve spare capacity so replayed append allocations follow the same
/// thresholds as the real site. Clone allocation/copying is timed separately.
fn clone_profile_site(site: &GrowingSite) -> Result<GrowingSite, String> {
    let clone_words = |source: &Vec<u16>| -> Result<Vec<u16>, String> {
        let mut words = Vec::new();
        words
            .try_reserve_exact(source.capacity())
            .map_err(|error| format!("Could not allocate profile group: {error}"))?;
        words.extend_from_slice(source);
        Ok(words)
    };
    match site {
        GrowingSite::Monomorphic(words) => Ok(GrowingSite::Monomorphic(clone_words(words)?)),
        GrowingSite::Polymorphic(source) => {
            let mut groups = Vec::new();
            groups
                .try_reserve_exact(source.capacity())
                .map_err(|error| format!("Could not allocate profile groups: {error}"))?;
            for words in source {
                groups.push(clone_words(words)?);
            }
            Ok(GrowingSite::Polymorphic(groups))
        }
    }
}

/// Replay sampled work in batches without modifying postings or real distances.
/// No per-entry clocks, unsafe pointers, or production-loop reordering are used.
/// Costs reflect scratch/cache behaviour and exclude hash-map insertion; they
/// must not be presented as an additive decomposition of the Rust core time.
fn profile_site_operations<IntT, F>(
    postings: &HashMap<IntT, GrowingSite>,
    entries: &[(IntT, u8)],
    sample_id: u16,
    threshold_words: usize,
    final_words: usize,
    now_ms: &mut F,
) -> Result<FineAggregationTimings, String>
where
    IntT: Eq + Hash,
    F: FnMut() -> f64,
{
    let preparation_started_at = now_ms();
    let mut buckets: [Vec<(GrowingSite, u8)>; 3] = std::array::from_fn(|_| Vec::new());
    for bucket in &mut buckets {
        bucket
            .try_reserve_exact(entries.len())
            .map_err(|error| format!("Could not allocate profile bucket: {error}"))?;
    }
    let mut new_bases = Vec::new();
    new_bases
        .try_reserve_exact(entries.len())
        .map_err(|error| format!("Could not allocate profile bases: {error}"))?;
    // There are at most 16 encoded base groups per site. References remain
    // valid because this entire replay precedes any mutation of the real map.
    let mut distance_groups: Vec<&[u16]> = Vec::new();
    distance_groups
        .try_reserve_exact(entries.len() * 16)
        .map_err(|error| format!("Could not allocate profile distances: {error}"))?;
    let mut timings = FineAggregationTimings {
        sampled_entries: entries.len(),
        ..FineAggregationTimings::default()
    };
    for (kmer, code) in entries {
        if let Some(site) = postings.get(kmer) {
            let bucket = match site {
                GrowingSite::Monomorphic(words) => {
                    if group_base(words) == *code {
                        0
                    } else {
                        distance_groups.push(words);
                        1
                    }
                }
                GrowingSite::Polymorphic(groups) => {
                    for words in groups {
                        if group_base(words) != *code {
                            distance_groups.push(words);
                        }
                    }
                    2
                }
            };
            buckets[bucket].push((clone_profile_site(site)?, *code));
            timings.operation_entries[bucket] += 1;
        } else {
            new_bases.push(*code);
            timings.operation_entries[3] += 1;
        }
    }
    timings.distance_groups = distance_groups.len();
    timings.distance_ids = distance_groups.iter().map(|group| group_len(group) as u64).sum();
    let mut scratch_distances = Vec::new();
    scratch_distances
        .try_reserve_exact(usize::from(sample_id))
        .map_err(|error| format!("Could not allocate profile distance row: {error}"))?;
    scratch_distances.resize(usize::from(sample_id), 0u32);
    let mut new_groups = Vec::new();
    new_groups
        .try_reserve_exact(new_bases.len())
        .map_err(|error| format!("Could not allocate profile new sites: {error}"))?;
    timings.preparation_ms = now_ms() - preparation_started_at;

    for (index, bucket) in buckets.iter_mut().enumerate() {
        if bucket.is_empty() {
            continue;
        }
        let started_at = now_ms();
        for (site, code) in bucket.iter_mut() {
            update_existing_site(site, *code, sample_id, threshold_words, final_words, |_| {})?;
        }
        std::hint::black_box(&*bucket);
        timings.operation_ms[index] = now_ms() - started_at;
    }
    if !new_bases.is_empty() {
        let started_at = now_ms();
        for code in &new_bases {
            new_groups.push(new_group(*code, sample_id)?);
        }
        std::hint::black_box(&new_groups);
        timings.operation_ms[3] = now_ms() - started_at;
    }
    if !distance_groups.is_empty() {
        let started_at = now_ms();
        for group in &distance_groups {
            accumulate_group_distances(group, &mut scratch_distances);
        }
        std::hint::black_box(&scratch_distances);
        timings.distance_ms = now_ms() - started_at;
    }
    // Scratch destruction is excluded from phase clocks. The caller measures
    // the whole replay, including preparation and destruction, as overhead.
    Ok(timings)
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
    pub(crate) fn add_sample_words<F>(
        &mut self,
        words: Vec<u32>,
        distance_row: &mut [u32],
        now_ms: &mut F,
    ) -> Result<AggregationPhaseTimings, String>
    where
        F: FnMut() -> f64,
    {
        if self.sealed {
            return Err("The alignment is closed to new samples.".into());
        }
        if distance_row.len() != self.sample_count {
            return Err("Distance row length does not match the existing sample count.".into());
        }
        let result = self.store_and_compare(&words, distance_row, now_ms);
        if result.is_err() {
            self.sealed = true;
        }
        result
    }

    fn store_and_compare<F>(
        &mut self,
        words: &[u32],
        distance_row: &mut [u32],
        now_ms: &mut F,
    ) -> Result<AggregationPhaseTimings, String>
    where
        F: FnMut() -> f64,
    {
        let SiteStorage::Growing(postings) = &mut self.storage else {
            return Err("The alignment is closed to new samples.".into());
        };
        let sample_id = self.sample_count as u16;
        let words_per_entry = std::mem::size_of::<IntT>() / std::mem::size_of::<u32>();
        let entries = words.len() / words_per_entry;
        let mut timings = AggregationPhaseTimings {
            total_entries: entries,
            ..AggregationPhaseTimings::default()
        };

        let mut profile_entries = Vec::new();
        let profile_enabled = profile_entries
            .try_reserve_exact(AGGREGATION_PROFILE_BLOCK_ENTRIES)
            .is_ok();

        let mut block_start = 0;
        while block_start < entries {
            let sampled_block = profile_enabled
                && block_start % AGGREGATION_PROFILE_SAMPLE_INTERVAL == 0;
            let entries_until_next_sample = if !profile_enabled {
                entries - block_start
            } else if sampled_block {
                AGGREGATION_PROFILE_BLOCK_ENTRIES
            } else {
                AGGREGATION_PROFILE_SAMPLE_INTERVAL
                    - (block_start % AGGREGATION_PROFILE_SAMPLE_INTERVAL)
            };
            let block_end = block_start
                .saturating_add(entries_until_next_sample)
                .min(entries);
            let block_entries = block_end - block_start;
            let first_word = block_start * words_per_entry;
            let end_word = block_end * words_per_entry;
            let block_words = &words[first_word..end_word];

            if sampled_block {
                timings.sampled_entries += block_entries;
                profile_entries.clear();

                let started_at = now_ms();
                for entry in block_words.chunks_exact(words_per_entry) {
                    profile_entries.push(unpack_entry::<IntT>(entry));
                }
                let sampled_unpack_ms = now_ms() - started_at;

                let started_at = now_ms();
                for (kmer, _) in &profile_entries {
                    if postings.contains_key(kmer) {
                        timings.sampled_lookup_hits += 1;
                    }
                }
                let sampled_hash_lookup_ms = now_ms() - started_at;
                timings.sampled_unpack_ms += sampled_unpack_ms;
                timings.sampled_hash_lookup_ms += sampled_hash_lookup_ms;

                let diagnostic_started_at = now_ms();
                let fine_result = profile_site_operations(
                    postings,
                    &profile_entries,
                    sample_id,
                    self.growth_threshold_words,
                    self.final_group_words,
                    now_ms,
                );
                timings.fine.diagnostic_ms += now_ms() - diagnostic_started_at;
                match fine_result {
                    Ok(fine) => {
                        timings.fine.sampled_entries += fine.sampled_entries;
                        timings.fine.preparation_ms += fine.preparation_ms;
                        for index in 0..4 {
                            timings.fine.operation_entries[index] += fine.operation_entries[index];
                            timings.fine.operation_ms[index] += fine.operation_ms[index];
                        }
                        timings.fine.distance_ms += fine.distance_ms;
                        timings.fine.distance_groups += fine.distance_groups;
                        timings.fine.distance_ids += fine.distance_ids;
                    }
                    Err(_) => {
                        // Diagnostic allocation failures must never abort the
                        // actual alignment or report partial replay estimates.
                        timings.fine.skipped_blocks += 1;
                    }
                }
            }

            let sampled_capacity_before = sampled_block.then(|| postings.capacity());
            let sampled_block_started_at = sampled_block.then(|| now_ms());
            for entry in block_words.chunks_exact(words_per_entry) {
                let (kmer, code) = unpack_entry::<IntT>(entry);
                let site_entry = postings.entry(kmer);
                match site_entry {
                    Entry::Occupied(mut entry) => {
                        update_existing_site(
                            entry.get_mut(),
                            code,
                            sample_id,
                            self.growth_threshold_words,
                            self.final_group_words,
                            |group| accumulate_group_distances(group, distance_row),
                        )?;
                    }
                    Entry::Vacant(site) => {
                        site.insert(GrowingSite::Monomorphic(new_group(code, sample_id)?));
                    }
                }
            }

            if let Some(started_at) = sampled_block_started_at {
                let sampled_block_ms = now_ms() - started_at;
                let map_capacity_grew = sampled_capacity_before
                    .is_some_and(|capacity| postings.capacity() > capacity);
                if map_capacity_grew {
                    // Exclude the whole sampled block from the scaled estimate:
                    // the resize cost cannot be isolated without per-entry clocks.
                    timings.sampled_resize_block_ms += sampled_block_ms;
                    timings.sampled_resize_blocks += 1;
                } else {
                    timings.sampled_update_entries += block_entries;
                    timings.sampled_block_work_ms += sampled_block_ms;
                }
            }

            block_start = block_end;
        }
        self.sample_count += 1;
        Ok(timings)
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
