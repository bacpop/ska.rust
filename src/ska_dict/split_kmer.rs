use crate::ska_dict::bit_encoding::*;
use std::borrow::Cow;

pub struct SplitKmer<'a> {
    k: usize,
    upper_mask: u64,
    lower_mask: u64,
    seq: Cow<'a, [u8]>,
    seq_len: usize,
    index: usize,
    upper: u64,
    lower: u64,
    middle_base: u8,
    rc: bool,
    rc_upper: u64,
    rc_lower: u64,
    rc_middle_base: u8,
}

impl<'a> SplitKmer<'a> {
    fn build(seq: &[u8], seq_len: usize, k: usize, idx: &mut usize) -> Option<(u64, u64, u8)> {
        if *idx + k >= seq_len {
            return None;
        }
        let mut upper: u64 = 0;
        let mut lower: u64 = 0;
        let mut middle_base: u8 = 0;
        let middle_idx = (k + 1) / 2 - 1;
        let mut i = 0;
        while i < k {
            if valid_base(seq[i + *idx]) {
                // Checks for N or n
                let next_base = encode_base(seq[i + *idx]);
                if i < middle_idx {
                    upper = upper << 2;
                    upper |= (next_base as u64) << (middle_idx * 2);
                } else if i > middle_idx {
                    lower = lower << 2;
                    lower |= next_base as u64;
                } else {
                    middle_base = next_base;
                }
                i += 1;
            } else {
                // Start again, skipping over N
                *idx += i + 1;
                if *idx + k >= seq_len {
                    return None;
                }
                upper = 0;
                lower = 0;
                middle_base = 0;
                i = 0;
            }
        }
        *idx += k - 1;
        return Some((upper, lower, middle_base));
    }

    fn update_rc(&mut self) {
        self.rc_upper = revcomp64_v2(self.lower, 30) & self.upper_mask;
        self.rc_middle_base = rc_base(self.middle_base);
        self.rc_lower = revcomp64_v2(self.upper, 30) & self.lower_mask;
    }

    fn roll_fwd(&mut self) -> bool {
        let mut success = false;
        self.index += 1;
        if self.index >= self.seq_len {
            return success;
        }
        let base = self.seq[self.index];
        if !valid_base(base) {
            let new_kmer = Self::build(&*self.seq, self.seq_len, self.k, &mut self.index);
            if new_kmer.is_some() {
                (self.upper, self.lower, self.middle_base) = new_kmer.unwrap();
                if self.rc {
                    self.update_rc();
                }
                success = true;
            }
        } else {
            let half_k: usize = (self.k - 1) / 2;
            self.upper = (self.upper << 2 | ((self.middle_base as u64) << (half_k * 2))) & self.upper_mask;
            self.middle_base = (self.lower >> (2 * (half_k - 1))) as u8;
            let new_base = encode_base(base);
            self.lower = (self.lower << 2 | (new_base as u64)) & self.lower_mask;
            if self.rc {
                self.rc_lower =
                    (self.rc_lower >> 2 | ((self.rc_middle_base as u64) << (2 * (half_k - 1)))) & self.lower_mask;
                self.rc_middle_base = rc_base(self.middle_base);
                self.rc_upper =
                    (self.rc_upper >> 2 | (rc_base(new_base) as u64) << (2*((half_k * 2) - 1))) & self.upper_mask;
            }
            success = true;
        }
        return success;
    }

    pub fn new(seq: Cow<'a, [u8]>, seq_len: usize, k: usize, rc: bool) -> Option<Self> {
        let (mut index, rc_upper, rc_lower, rc_middle_base) = (0, 0, 0, 0);
        let first_kmer = Self::build(&*seq, seq_len, k, &mut index);
        if first_kmer.is_some() {
            let (upper, lower, middle_base) = first_kmer.unwrap();
            let (upper_mask, lower_mask) = generate_masks(k);
            let mut split_kmer = Self {
                k,
                upper_mask,
                lower_mask,
                seq_len,
                seq,
                upper,
                lower,
                middle_base,
                rc,
                rc_upper,
                rc_lower,
                rc_middle_base,
                index,
            };
            if rc {
                split_kmer.update_rc();
            }
            return Some(split_kmer);
        } else {
            return None;
        }
    }

    pub fn get_curr_kmer(&self) -> (u64, u8) {
        let split_kmer = self.upper | self.lower;
        if self.rc {
            let rc_split_kmer = self.rc_upper | self.rc_lower;
            if split_kmer > rc_split_kmer {
                return (rc_split_kmer, decode_base(self.rc_middle_base));
            }
        }
        return (split_kmer, decode_base(self.middle_base));
    }

    pub fn get_next_kmer(&mut self) -> Option<(u64, u8)> {
        let next = self.roll_fwd();
        match next {
            true => Some(self.get_curr_kmer()),
            false => None,
        }
    }
}
