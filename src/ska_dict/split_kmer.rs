//! Create split k-mers from sequences by rolling through input
//!
//! The [`SplitKmer`] struct stores a reference to then input sequence, the
//! necessary values to extract k-mers, and the current split k-mer:middle base
//! (as well as its reverse complement)
//!
//! Use [`SplitKmer::new()`] to initialise with sequence input from [`needletail`],
//! and [`SplitKmer::get_next_kmer()`] to advance through the sequence.
//!
//! (This could be made into a rust iterator, but I didn't do this when I wrote
//! it as I didn't know how yet.)

use crate::ska_dict::bit_encoding::*;
use std::borrow::Cow;
use std::cmp::Ordering;

/// Struct to generate all split k-mers from an input sequence
///
/// Holds reference to input sequence, current encoded k-mer, and other
/// information (k-mer size, masks etc.)
#[derive(Debug)]
pub struct SplitKmer<'a, IntT> {
    /// K-mer size
    k: usize,
    /// Mask to extract upper k-mer
    upper_mask: IntT,
    /// Mask to extract lower k-mer
    lower_mask: IntT,
    /// Reference to input sequence
    seq: Cow<'a, [u8]>,
    /// Size of seq
    seq_len: usize,
    /// Reference to sequence quality scores
    qual: Option<&'a [u8]>,
    /// Minimum quality score to allow in a split k-mers
    min_qual: u8,
    /// Current index in input sequence
    index: usize,
    /// Current upper part of split k-mer
    upper: IntT,
    /// Current lower part of split k-mer
    lower: IntT,
    /// Current middle base
    middle_base: u8,
    /// Whether reverse complements are being used
    rc: bool,
    /// Current upper part of split k-mer, reverse complemented
    rc_upper: IntT,
    /// Current lower part of split k-mer, reverse complemented
    rc_lower: IntT,
    /// Current middle base, reverse complemented
    rc_middle_base: u8,
}

impl<'a, IntT: for<'b> RevComp<'b>> SplitKmer<'a, IntT> {
    /// Quality score is at least minimum.
    #[inline(always)]
    fn valid_qual(idx: usize, qual: Option<&'a [u8]>, min_qual: u8) -> bool {
        match qual {
            Some(qual_seq) => qual_seq[idx] - 33 > min_qual, // ASCII encoding starts from b'!' = 33
            None => true,
        }
    }

    /// Build a new split k-mer at the given index.
    ///
    /// Called when initialised, or after skipping unknown bases. Returns
    /// [`None`] if the end of the input has been reached.
    fn build(
        seq: &[u8],
        seq_len: usize,
        qual: Option<&'a [u8]>,
        k: usize,
        idx: &mut usize,
        min_qual: u8,
    ) -> Option<(IntT, IntT, u8)> {
        if *idx + k >= seq_len {
            return None;
        }
        let mut upper = IntT::zero_init();
        let mut lower = IntT::zero_init();
        let mut middle_base: u8 = 0;
        let middle_idx = (k + 1) / 2 - 1;
        let mut i = 0;
        while i < k {
            if valid_base(seq[i + *idx]) && Self::valid_qual(i + *idx, qual, min_qual) {
                // Checks for N or n
                let next_base = encode_base(seq[i + *idx]);
                match i.cmp(&middle_idx) {
                    Ordering::Greater => {
                        lower <<= 2;
                        lower |= IntT::from_encoded_base(next_base);
                    }
                    Ordering::Less => {
                        upper <<= 2;
                        upper |= IntT::from_encoded_base(next_base) << (middle_idx * 2);
                    }
                    Ordering::Equal => {
                        middle_base = next_base;
                    }
                }
                i += 1;
            } else {
                // Start again, skipping over N
                *idx += i + 1;
                if *idx + k >= seq_len {
                    return None;
                }
                upper = IntT::zero_init();
                lower = IntT::zero_init();
                middle_base = 0;
                i = 0;
            }
        }
        *idx += k - 1;
        Some((upper, lower, middle_base))
    }

    /// Update the stored reverse complement using the stored split-k and middle base
    fn update_rc(&mut self) {
        self.rc_upper = self.lower.rev_comp(self.k - 1) & self.upper_mask;
        self.rc_middle_base = rc_base(self.middle_base);
        self.rc_lower = self.upper.rev_comp(self.k - 1) & self.lower_mask;
    }

    /// Move forward to the next valid split k-mer
    ///
    /// Usually the next base, but if an N skips over it.
    /// If end of sequence encountered then returns `false`.
    fn roll_fwd(&mut self) -> bool {
        let mut success = false;
        self.index += 1;
        if self.index >= self.seq_len {
            return success;
        }
        let base = self.seq[self.index];
        if !valid_base(base) || !Self::valid_qual(self.index, self.qual, self.min_qual) {
            let new_kmer = Self::build(
                &self.seq,
                self.seq_len,
                self.qual,
                self.k,
                &mut self.index,
                self.min_qual,
            );
            if let Some(kmer_tuple) = new_kmer {
                (self.upper, self.lower, self.middle_base) = kmer_tuple;
                if self.rc {
                    self.update_rc();
                }
                success = true;
            }
        } else {
            let half_k: usize = (self.k - 1) / 2;
            self.upper = (self.upper << 2
                | (IntT::from_encoded_base(self.middle_base) << (half_k * 2)))
                & self.upper_mask;
            self.middle_base = (self.lower >> (2 * (half_k - 1))).as_u8();
            let new_base = encode_base(base);
            self.lower = ((self.lower << 2) | (IntT::from_encoded_base(new_base))) & self.lower_mask;
            if self.rc {
                self.rc_lower = (self.rc_lower >> 2
                    | ((IntT::from_encoded_base(self.rc_middle_base)) << (2 * (half_k - 1))))
                    & self.lower_mask;
                self.rc_middle_base = rc_base(self.middle_base);
                self.rc_upper = (self.rc_upper >> 2
                    | (IntT::from_encoded_base(rc_base(new_base))) << (2 * ((half_k * 2) - 1)))
                    & self.upper_mask;
            }
            success = true;
        }
        success
    }

    /// Create a [`SplitKmer`] iterator given reference to sequence input.
    ///
    /// Sequence, length and quality come from [`needletail`].
    ///
    /// Returns [`None`] if no valid split k-mers found in input (e.g. too short,
    /// no sequence, too many Ns).
    pub fn new(
        seq: Cow<'a, [u8]>,
        seq_len: usize,
        qual: Option<&'a [u8]>,
        k: usize,
        rc: bool,
        min_qual: u8,
    ) -> Option<Self> {
        let mut index = 0;
        let first_kmer = Self::build(&seq, seq_len, qual, k, &mut index, min_qual);
        if let Some((upper, lower, middle_base)) = first_kmer {
            let (lower_mask, upper_mask) = IntT::generate_masks(k);
            let mut split_kmer = Self {
                k,
                upper_mask,
                lower_mask,
                seq_len,
                seq,
                qual,
                min_qual,
                upper,
                lower,
                middle_base,
                rc,
                rc_upper: IntT::zero_init(),
                rc_lower: IntT::zero_init(),
                rc_middle_base: 0,
                index,
            };
            if rc {
                split_kmer.update_rc();
            }
            Some(split_kmer)
        } else {
            None
        }
    }

    /// Get the current split k-mer
    ///
    /// Returns split k-mer, middle base, and whether the reverse complement
    pub fn get_curr_kmer(&self) -> (IntT, u8, bool) {
        let split_kmer = self.upper | self.lower;
        // Some of the most useful prints for debugging left as comments here
        //let (upper, lower) = decode_kmer(self.k, split_kmer, self.upper_mask, self.lower_mask);
        //println!("{} {} {}", upper, lower, self.middle_base);
        if self.rc {
            let rc_split_kmer = self.rc_upper | self.rc_lower;
            //let (upper_rc, lower_rc) = decode_kmer(self.k, rc_split_kmer, self.upper_mask, self.lower_mask);
            //println!("{} {} {}", upper_rc, lower_rc, self.rc_middle_base);
            if split_kmer > rc_split_kmer {
                return (rc_split_kmer, self.rc_middle_base, true);
            }
        }
        (split_kmer, self.middle_base, false)
    }

    /// Get the next split k-mer in the sequence
    ///
    /// Returns split k-mer, middle base, and whether the reverse complement
    /// or [`None`] if no next k-mer
    pub fn get_next_kmer(&mut self) -> Option<(IntT, u8, bool)> {
        let next = self.roll_fwd();
        match next {
            true => Some(self.get_curr_kmer()),
            false => None,
        }
    }

    /// Get the index in the sequence of the current middle base
    pub fn get_middle_pos(&self) -> usize {
        let middle_idx = (self.k + 1) / 2 - 1;
        self.index - middle_idx
    }
}
