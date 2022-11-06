
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
    fn build(seq: Cow<'a, [u8]>, seq_len: usize, k: usize, idx: &mut usize) -> Option<(u64, u64, u8)> {
        if *idx + k >= seq_len {
            return None;
        }
        let mut start: u64 = 0;
        let mut end: u64 = 0;
        let mut middle_base: u8 = 0;
        let middle_idx = (k + 1) / 2;
        for i in 0..k {
            if seq[i + *idx] & 0xF != 14 { // Checks for N or n
                let next_base = encode_base(seq[i + *idx]) as u64;
                if i < middle_idx {
                    end = end << 2;
                    end |= next_base;
                } else if i > middle_idx {
                    start = start << 2;
                    start |= next_base << ((middle_idx + 1) * 2);
                } else {
                    middle_base = next_base;
                }
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
        return Some(start, end, middle_base);
    }

    fn update_rc(&mut self) {
        self.rc_start = revcomp64_v2(self.end, 15);
        self.rc_middle_base = rc_base(self.middle_base);
        self.rc_end = revcomp64_v2(self.start, 15);
    }

    fn roll_fwd(&mut self) -> bool {
        let success = false;
        self.index += 1;
        if self.index >= self.seq_len {
            return success;
        }
        let base = self.seq[self.index];
        if !valid_base(base) {
            let new_kmer = Self::build(self.seq, 31, &self.idx);
            if new_kmer.is_some() {
                (self.start, self.end, self.middle_base) = new_kmer.unwrap_unchecked();
                if self.rc {
                    self.update_rc();
                }
                success = true;
            }
        } else {
            self.start = (self.start << 2 | (self.middle_base << 31)) & upper_mask;
            self.middle_base = self.end >> 30;
            let new_base = encode_base(base);
            self.end = (self.end << 2 | new_base) & lower_mask;
            if self.rc {
                self.rc_start = (self.rc_start >> 2 | self.rc_middle_base << 62) & lower_mask;
                self.rc_middle_base = rc_base(self.middle_base);
                self.rc_end = (self.rc_end << 2 | rc_base(new_base)) & upper_mask;
            }
            success = true;
        }
        return success;
    }

    pub fn new(seq: &u8, seq_len: usize, rc: bool) -> Self {
        let (mut idx, rc_start, rc_end, rc_middle_base) = 0;
        let k = 31;
        let (start, end, middle_base) = Self::build(seq, k, &idx).expect("Sequence too short or too many Ns");
        let mut split_kmer = Self {k, seq, start, end, middle_base, rc, rc_start, rc_end, rc_middle_base, idx};
        if rc {
            split_kmer.update_rc();
        }
        return split_kmer;
    }

    pub fn get_curr_kmer(&self) -> (u64, u8) {
        let split_kmer = self.start | self.end;
        if self.rc {
            let rc_split_kmer = self.rc_start | self.rc_end;
            if split_kmer > rc_split_kmer {
                return (self.rc_split_kmer, self.rc_middle_base);
            }
        }
        return (self.split_kmer, self.middle_base);
    }

    pub fn get_next_kmer(&mut self) -> Option<(u64, u8)> {
        let next = self.roll_fwd();
        match next {
            true => Some(self.get_curr_kmer()),
            false => None
        }
    }
}
