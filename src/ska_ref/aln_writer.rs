//! Generate a linear pseudoalignment.
//!
//! An empty sequence with the same length is set as the sequence
//! Matches are written as the first half, plus the middle base, e.g. for 7-mer
//! ACTGCTA then ACTG would be written. The first half is copied over from the
//! reference. Matches are skipped over until the next `k / 2 + 1` position
//! which is where more sequence can then be written. If bases are missing
//! or at the end of the contig these are caught up by adding the skipped
//! matches and the front half of the k-mer.
//!
//! All sequence is stored in memory. Calling [`AlnWriter::get_seq()`] finalises
//! and writes to the end of the last contig.

use crate::ska_dict::bit_encoding::is_ambiguous;

/// Stores indexes which keep track of the writing, the output sequence and a
/// reference to the reference sequences which are written in to flanks of
/// each match.
#[derive(Clone)]
pub struct AlnWriter<'a> {
    /// Next position where it is valid to write a match.
    next_pos: usize,
    /// The current chromsome.
    curr_chrom: usize,
    /// The latest position where a base has been mapped (may not have been written).
    last_mapped: usize,
    /// The latest position in `seq_out` that has been written to.
    last_written: usize,
    /// An offset which is added to convert chromosome positions to concatenated position.
    chrom_offset: usize,
    /// The reference sequences.
    ref_seq: &'a Vec<Vec<u8>>,
    /// The output alignment.
    seq_out: Vec<u8>,
    /// The size of the flanking region for each split k-mer match.
    half_split_len: usize,
    /// Whether the finalise function has been run, filling to the end of the
    /// final contig.
    finalised: bool,
    /// Repeats to mask, which happens at finalisation
    repeat_regions: &'a Vec<usize>,
    /// Whether to treat all ambiguous bases as Ns
    mask_ambig: bool,
    /// Buffer for middle bases, which are written at finalisation.
    _middle_out: Vec<(u8, usize)>,
}

impl<'a> AlnWriter<'a> {
    /// Create a new [`AlnWriter`] taking the reference sequence mapped against
    /// and the k-mer size used for the mapping
    pub fn new(
        ref_seq: &'a Vec<Vec<u8>>,
        k: usize,
        repeat_regions: &'a Vec<usize>,
        mask_ambig: bool,
    ) -> Self {
        let total_size = ref_seq.iter().map(|x| x.len()).sum();
        let half_split_len = (k - 1) / 2;
        Self {
            next_pos: half_split_len,
            curr_chrom: 0,
            last_mapped: 0,
            last_written: 0,
            chrom_offset: 0,
            ref_seq,
            seq_out: vec![b'-'; total_size],
            half_split_len,
            finalised: false,
            repeat_regions,
            mask_ambig,
            _middle_out: Vec::new(),
        }
    }

    /// Get the total length of the concatenated sequence output
    pub fn total_size(&self) -> usize {
        self.seq_out.len()
    }

    // Fill fwd bases, and any skipped over.
    // e.g. with split 7-mers perfectly matching
    // CCGA AAGT
    //    1   23
    // 1: Last written
    // 2: Last mapped
    // 3: Current mapped position
    // Don't write anything
    // CCGA[TT]-----AAGT
    //    1 A2   B     3
    // Want to write between
    // A = last written + 1
    // B = A + half-k + (last mapped - last written)
    // maximum is provided as the base about to be written next, this will not
    // write over this position
    fn fill_fwd_bases(&mut self, maximum: usize) {
        // bases at the end of last valid match not yet written
        if self.last_written > 0 {
            let last_match_overhang =
                (self.last_mapped + self.half_split_len).saturating_sub(self.last_written);
            let start = self.last_written + 1;
            let end = usize::min(start + last_match_overhang, maximum);
            if end > start {
                self.seq_out[(start + self.chrom_offset)..(end + self.chrom_offset)]
                    .copy_from_slice(&self.ref_seq[self.curr_chrom][start..end]);
                self.last_written = end;
            }
        }
    }

    // Fill up to the end of the contig and update indices
    fn fill_contig(&mut self) {
        let chrom_length = self.ref_seq[self.curr_chrom].len();
        self.fill_fwd_bases(chrom_length);
        self.chrom_offset += chrom_length;
        self.curr_chrom += 1;
        self.next_pos = self.half_split_len;
    }

    /// Write the first half of a split k-mer (by copying from the reference)
    /// and its middle base, also filling any skipped bases as appropiate.
    ///
    /// Call in a loop over mapped positions, usually excluding '-' base which
    /// haven't actually been mapped
    pub fn write_split_kmer(&mut self, mapped_pos: usize, mapped_chrom: usize, base: u8) {
        while mapped_chrom > self.curr_chrom {
            self.fill_contig();
        }
        // Middle bases may clash with the flanks in complex repeats, which then
        // are copied from reference. Deal with these in `finalise`.
        self._middle_out.push((
            if is_ambiguous(base) && self.mask_ambig {
                b'N'
            } else {
                base
            },
            mapped_pos + self.chrom_offset,
        ));

        if mapped_pos < self.next_pos {
            self.last_mapped = mapped_pos;
        } else {
            // Write bases between last match and this one
            if mapped_pos > self.next_pos {
                self.fill_fwd_bases(mapped_pos - self.half_split_len);
            }

            // First half of split k-mer
            let start = mapped_pos - self.half_split_len;
            let end = mapped_pos;
            self.seq_out[(start + self.chrom_offset)..(end + self.chrom_offset)]
                .copy_from_slice(&self.ref_seq[self.curr_chrom][start..end]);

            // update indices
            self.next_pos = mapped_pos + self.half_split_len + 1;
            self.last_mapped = mapped_pos;
            self.last_written = mapped_pos;
        }
    }

    /// Fills to the end of the final contig
    /// Should only be called after the last call to [`AlnWriter::write_split_kmer()`].
    pub fn finalise(&mut self) {
        if !self.finalised {
            while self.curr_chrom < self.ref_seq.len() {
                self.fill_contig();
            }
            // Make sure any ambiguous bases are correct
            for (middle_base, middle_pos) in &self._middle_out {
                self.seq_out[*middle_pos] = *middle_base;
            }
            // Mask repeats
            for repeat_idx in self.repeat_regions {
                if self.seq_out[*repeat_idx] != b'-' {
                    self.seq_out[*repeat_idx] = b'N';
                }
            }
            self.finalised = true;
        }
    }

    /// Retrieve the written sequence. Calls [`AlnWriter::finalise()`] if not already called.
    pub fn get_seq(&'a mut self) -> &'a [u8] {
        self.finalise();
        self.seq_out.as_slice()
    }
}
