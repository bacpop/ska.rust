
use std::borrow::Cow;
use crate::ska_dict::bit_encoding::*;

pub struct SplitKmer<'a> {
    k: usize,
    seq: Cow<'a, [u8]>,
    seq_len: usize,
    start: u64, // upper bits
    end: u64, // lower bits
    middle_base: u8,
    rc: bool,
    rc_start: u64,
    rc_end: u64,
    rc_middle_base: u8,
    index: usize
}

impl<'a> SplitKmer<'a> {
    fn build(seq: &[u8], seq_len: usize, k: usize, idx: &mut usize) -> Option<(u64, u64, u8)> {
        if *idx + k >= seq_len {
            return None;
        }
        let mut start: u64 = 0;
        let mut end: u64 = 0;
        let mut middle_base: u8 = 0;
        let middle_idx = (k + 1) / 2;
        let mut i = 0;
        while i < k {
            if valid_base(seq[i + *idx]) { // Checks for N or n
                let next_base = encode_base(seq[i + *idx]);
                if i < middle_idx {
                    end = end << 2;
                    end |= next_base as u64;
                } else if i > middle_idx {
                    start = start << 2;
                    start |= (next_base as u64) << ((middle_idx + 1) * 2);
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
                start = 0;
                end = 0;
                middle_base = 0;
                i = 0;
            }
        }
        *idx += k + 1;
        return Some((start, end, middle_base));
    }

    fn update_rc(&mut self) {
        self.rc_start = revcomp64_v2(self.end, 15);
        self.rc_middle_base = rc_base(self.middle_base);
        self.rc_end = revcomp64_v2(self.start, 15);
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
                (self.start, self.end, self.middle_base) = new_kmer.unwrap();
                if self.rc {
                    self.update_rc();
                }
                success = true;
            }
        } else {
            self.start = (self.start << 2 | ((self.middle_base as u64) << 31)) & UPPER_MASK;
            self.middle_base = (self.end >> 30) as u8;
            let new_base = encode_base(base);
            self.end = (self.end << 2 | (new_base as u64)) & LOWER_MASK;
            if self.rc {
                self.rc_start = (self.rc_start >> 2 | (self.rc_middle_base as u64) << 62) & LOWER_MASK;
                self.rc_middle_base = rc_base(self.middle_base);
                self.rc_end = (self.rc_end << 2 | (rc_base(new_base) as u64)) & UPPER_MASK;
            }
            success = true;
        }
        return success;
    }

    pub fn new(seq: Cow<'a, [u8]>, seq_len: usize, rc: bool) -> Option<Self> {
        let (mut index, rc_start, rc_end, rc_middle_base) = (0, 0, 0, 0);
        let k = 31;
        let first_kmer = Self::build(&*seq, seq_len, k, &mut index);
        if first_kmer.is_some() {
            let (start, end, middle_base) = first_kmer.unwrap();
            let mut split_kmer = Self {k, seq_len, seq, start, end, middle_base, rc, rc_start, rc_end, rc_middle_base, index};
            if rc {
                split_kmer.update_rc();
            }
            return Some(split_kmer);
        } else {
            return None;
        }
    }

    pub fn get_curr_kmer(&self) -> (u64, u8) {
        let split_kmer = self.start | self.end;
        if self.rc {
            let rc_split_kmer = self.rc_start | self.rc_end;
            if split_kmer > rc_split_kmer {
                return (rc_split_kmer, self.rc_middle_base);
            }
        }
        return (split_kmer, self.middle_base);
    }

    pub fn get_next_kmer(&mut self) -> Option<(u64, u8)> {
        let next = self.roll_fwd();
        match next {
            true => Some(self.get_curr_kmer()),
            false => None
        }
    }
}
