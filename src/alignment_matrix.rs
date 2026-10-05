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

pub(crate) const MAX_SAMPLES: usize = 4096;

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
    if words.len() % width != 0 {
        return Err("The packed sample contains an incomplete entry.".into());
    }
    Ok(())
}

/// Room for one more call, doubling as `Vec` would but never past `limit` calls, so a site's last
/// growth goes straight to its final size.
fn reserve_call(calls: &mut Vec<u16>, limit: usize) -> Result<(), String> {
    if calls.len() < calls.capacity() {
        return Ok(());
    }
    let target = (calls.capacity() * 2)
        .max(4)
        .min(limit)
        .max(calls.len() + 1);
    calls
        .try_reserve_exact(target - calls.len())
        .map_err(|error| format!("Could not extend site calls: {error}"))
}

fn increment_distance(distance: &mut u32) -> Result<(), String> {
    *distance = distance
        .checked_add(1)
        .ok_or_else(|| "SNP distance exceeds u32::MAX.".to_string())?;
    Ok(())
}

#[derive(Debug, Clone)]
struct SparseSite {
    calls: Vec<u16>,
    cursor: usize,
}

/// A site while samples are still being added.
#[derive(Debug, Clone)]
struct GrowingSite {
    calls: Vec<u16>,
    /// Some stored call differs from the first one's base; set once and never cleared.
    polymorphic: bool,
}

#[derive(Debug, Clone)]
enum SiteStorage<IntT> {
    Growing(HashMap<IntT, GrowingSite>),
    Frozen(Vec<SparseSite>),
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
            return Err("Sparse alignment supports at most 4096 samples.".into());
        }
        validate_words(words, width)
    }

    /// Store one sample and return its distances to the earlier ones, visiting each split k-mer once.
    /// A failure after storing has begun seals the matrix, as it may hold part of the sample.
    pub(crate) fn add_sample_words(
        &mut self,
        words: Vec<u32>,
        words_per_entry: usize,
    ) -> Result<Vec<u32>, String> {
        self.check_addition(&words, words_per_entry)?;
        let result = self.store_and_compare(&words, words_per_entry);
        if result.is_err() {
            self.sealed = true;
        }
        result
    }

    fn store_and_compare(
        &mut self,
        words: &[u32],
        words_per_entry: usize,
    ) -> Result<Vec<u32>, String> {
        let SiteStorage::Growing(postings) = &mut self.storage else {
            return Err("The alignment is closed to new samples.".into());
        };
        let sample_id = self.sample_count as u16;
        let limit = self.expected_samples;
        let mut row = Vec::<u32>::new();
        row.try_reserve_exact(self.sample_count + 1)
            .map_err(|error| format!("Could not allocate distance row: {error}"))?;
        row.resize(self.sample_count + 1, 0);
        for entry in words.chunks_exact(words_per_entry) {
            let (kmer, code) = unpack_entry::<IntT>(entry)?;
            let call = (sample_id << 4) | u16::from(code);
            // Reserving first keeps the insertion below fallible; with room to spare it is only a check.
            postings
                .try_reserve(1)
                .map_err(|error| format!("Could not extend site index: {error:?}"))?;
            match postings.entry(kmer) {
                Entry::Occupied(mut site) => {
                    let site = site.get_mut();
                    if site.calls.last().is_some_and(|last| last >> 4 == sample_id) {
                        return Err("A sample contains a duplicate split k-mer.".into());
                    }
                    // A site whose calls all share the new base adds no distance, so only
                    // polymorphic sites (or one becoming so) are walked.
                    if site.polymorphic || (site.calls[0] & 0xF) as u8 != code {
                        for &previous in &site.calls {
                            let previous_base = (previous & 0xF) as u8;
                            if previous_base != GAP && previous_base != code {
                                increment_distance(&mut row[usize::from(previous >> 4)])?;
                            }
                        }
                        site.polymorphic = true;
                    }
                    reserve_call(&mut site.calls, limit)?;
                    site.calls.push(call);
                }
                Entry::Vacant(site) => {
                    let mut calls = Vec::new();
                    calls
                        .try_reserve_exact(1)
                        .map_err(|error| format!("Could not allocate site calls: {error}"))?;
                    calls.push(call);
                    site.insert(GrowingSite {
                        calls,
                        polymorphic: false,
                    });
                }
            }
        }
        self.sample_count += 1;
        Ok(row)
    }

    /// Move call vectors without copying their payloads, then discard k-mer keys.
    fn freeze(&mut self) -> Result<(), String> {
        self.sealed = true;
        if let SiteStorage::Growing(postings) = &self.storage {
            let mut sites = Vec::new();
            sites
                .try_reserve_exact(postings.len())
                .map_err(|error| format!("Could not allocate sparse site array: {error}"))?;
            let SiteStorage::Growing(postings) =
                std::mem::replace(&mut self.storage, SiteStorage::Released)
            else {
                unreachable!()
            };
            for (_key, site) in postings {
                sites.push(SparseSite {
                    calls: site.calls,
                    cursor: 0,
                });
            }
            self.storage = SiteStorage::Frozen(sites);
        }
        if matches!(self.storage, SiteStorage::Released) {
            return Err("Alignment sites have already been exported and released.".into());
        }
        Ok(())
    }

    /// Reset cursors on every attempt, permitting retries after writer errors.
    fn write_fasta<W: Write>(
        &mut self,
        names: &[String],
        writer: &mut W,
        on_sample: &dyn Fn(usize),
    ) -> io::Result<()> {
        let SiteStorage::Frozen(sites) = &mut self.storage else {
            return Err(io::Error::other("The sparse array is not available."));
        };
        if names.len() != self.sample_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Sample names do not match the matrix.",
            ));
        }
        for site in sites.iter_mut() {
            site.cursor = 0;
        }
        let mut buffer = [0u8; 8192];
        for (sample, name) in names.iter().enumerate() {
            writer.write_all(b"> ")?;
            writer.write_all(name.as_bytes())?;
            writer.write_all(b"\n")?;
            let mut buffered = 0;
            for site in sites.iter_mut() {
                let code = match site.calls.get(site.cursor) {
                    Some(&call) if usize::from(call >> 4) == sample => {
                        site.cursor += 1;
                        (call & 0xF) as u8
                    }
                    _ => GAP,
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

    #[test]
    fn sparse_distance_rows_match_map_intersections() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut prior = Vec::<HashMap<u64, u8>>::new();
        for sample_id in 0..24 {
            let sample: HashMap<_, _> = (0..90u64)
                .filter(|key| (key * 17 + sample_id * 13) % 7 < 3)
                .map(|key| (key, decode_symbol(((key + sample_id) % 16) as u8)))
                .collect();
            let packed = words(&sample);
            let row = matrix.add_sample_words(packed, 2).unwrap();
            for (j, other) in prior.iter().enumerate() {
                let expected = sample
                    .iter()
                    .filter(|(key, base)| {
                        other
                            .get(*key)
                            .is_some_and(|old| *old != b'-' && old != *base)
                    })
                    .count() as u32;
                assert_eq!(row[j], expected);
            }
            assert_eq!(row.last(), Some(&0));
            prior.push(sample);
        }
    }

    #[test]
    fn invariant_sites_skip_comparisons_without_changing_rows() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut prior = Vec::<HashMap<u64, u8>>::new();
        for sample_id in 0..20u64 {
            let sample: HashMap<_, _> = (0..60u64)
                .filter(|key| (key + sample_id) % 9 != 0)
                .map(|key| {
                    let variant = key % 10 == 0 && sample_id % 3 == 0;
                    (key, if variant { b'T' } else { b'A' })
                })
                .collect();
            let row = matrix.add_sample_words(words(&sample), 2).unwrap();
            for (j, other) in prior.iter().enumerate() {
                let expected = sample
                    .iter()
                    .filter(|(key, base)| other.get(*key).is_some_and(|old| old != *base))
                    .count() as u32;
                assert_eq!(row[j], expected);
            }
            prior.push(sample);
        }
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        assert!(postings
            .iter()
            .all(|(key, site)| site.polymorphic == (key % 10 == 0)));
    }

    #[test]
    fn site_capacity_stops_at_the_expected_sample_count() {
        let shared = || pack_entry_u64(5, b'A').to_vec();
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.set_expected_samples(6);
        for _ in 0..6 {
            matrix.add_sample_words(shared(), 2).unwrap();
        }
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        // 1, then 4, then 6 rather than 8.
        assert_eq!(postings[&5].calls.capacity(), 6);

        // An underestimate still grows, one call beyond the limit at a time.
        let mut matrix = AlignmentMatrix::<u64>::default();
        matrix.set_expected_samples(2);
        for _ in 0..5 {
            matrix.add_sample_words(shared(), 2).unwrap();
        }
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        assert_eq!(postings[&5].calls.len(), 5);
    }

    #[test]
    fn zero_call_is_valid_and_sample_limit_is_checked() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        assert_eq!(
            matrix
                .add_sample_words(pack_entry_u64(9, b'A').to_vec(), 2)
                .unwrap(),
            vec![0]
        );
        assert_eq!(
            matrix
                .add_sample_words(pack_entry_u64(9, b'C').to_vec(), 2)
                .unwrap(),
            vec![1, 0]
        );
        for _ in 2..MAX_SAMPLES - 1 {
            matrix.add_sample_words(Vec::new(), 2).unwrap();
        }
        matrix
            .add_sample_words(pack_entry_u64(9, b'N').to_vec(), 2)
            .unwrap();
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        assert_eq!(postings[&9].calls, vec![0, (1 << 4) | 1, (4095 << 4) | 14]);
        assert!(postings[&9].polymorphic);
        assert!(matrix.add_sample_words(Vec::new(), 2).is_err());
    }

    #[test]
    fn large_keys_ambiguity_and_empty_samples_work() {
        let mut matrix = AlignmentMatrix::<u128>::default();
        let key = (1u128 << 100) + 27;
        assert_eq!(
            matrix
                .add_sample_words(pack_entry_u128(key, b'R').to_vec(), 4)
                .unwrap(),
            vec![0]
        );
        assert_eq!(matrix.add_sample_words(Vec::new(), 4).unwrap(), vec![0, 0]);
        assert_eq!(
            matrix
                .add_sample_words(pack_entry_u128(key, b'R').to_vec(), 4)
                .unwrap(),
            vec![0, 0, 0]
        );
        assert_eq!(
            matrix
                .add_sample_words(pack_entry_u128(key, b'A').to_vec(), 4)
                .unwrap(),
            vec![1, 0, 1, 0]
        );
        assert!(matrix.add_sample_words(vec![1], 4).is_err());
        assert!(matrix.add_sample_words(Vec::new(), 0).is_err());
        let mut distance = u32::MAX;
        assert!(increment_distance(&mut distance).is_err());
    }

    #[test]
    fn failed_addition_seals_the_matrix() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let mut duplicate = pack_entry_u64(3, b'A').to_vec();
        duplicate.extend(pack_entry_u64(3, b'C'));
        assert!(matrix.add_sample_words(duplicate, 2).is_err());
        assert!(matrix.is_sealed());
        assert!(matrix.add_sample_words(Vec::new(), 2).is_err());
    }

    #[test]
    fn freeze_moves_buffers_and_export_recovers_gaps_and_symbols() {
        let mut matrix = AlignmentMatrix::<u64>::default();
        let first: HashMap<_, _> = (0..16)
            .map(|code| (code as u64, decode_symbol(code)))
            .collect();
        matrix.add_sample_words(words(&first), 2).unwrap();
        matrix
            .add_sample_words(pack_entry_u64(0, b'G').to_vec(), 2)
            .unwrap();
        let SiteStorage::Growing(postings) = &matrix.storage else {
            panic!()
        };
        let pointers: Vec<_> = postings.values().map(|site| site.calls.as_ptr()).collect();
        matrix.freeze().unwrap();
        let SiteStorage::Frozen(sites) = &matrix.storage else {
            panic!()
        };
        assert!(sites
            .iter()
            .all(|site| pointers.contains(&site.calls.as_ptr())));
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
        assert!(matrix.add_sample_words(Vec::new(), 2).is_err());
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
        matrix
            .add_sample_words(pack_entry_u64(1, b'A').to_vec(), 2)
            .unwrap();
        matrix
            .add_sample_words(pack_entry_u64(1, b'T').to_vec(), 2)
            .unwrap();
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
}
