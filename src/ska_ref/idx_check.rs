//! Keep track of current chromosome during iteration over a concatenated
//! alignment.
//!
//! Simple counter with state to avoid doing range lookups (e.g. with an
//! interval tree) at every position. Used when converting alignment to VCF.

/// Builds the contig change indices and keeps track of the current chromosome.
///
/// NB: not robust to leaps in the passed index, as chromosome is only incremented
/// by one (so works for linear iteration over pos only)
#[derive(Debug)]
pub struct IdxCheck {
    end_coor: Vec<usize>,
}

impl IdxCheck {
    pub fn new(ref_seq: &[Vec<u8>]) -> Self {
        let mut end_coor = Vec::new();

        let mut cum_pos = 0;
        for chrom in ref_seq {
            cum_pos += chrom.len();
            end_coor.push(cum_pos);
        }

        Self { end_coor }
    }

    pub fn iter(&self) -> IdxCheckIter<'_> {
        IdxCheckIter {
            end_coor: &self.end_coor,
            current_chr: 0,
            idx: 0,
        }
    }
}

pub struct IdxCheckIter<'a> {
    end_coor: &'a Vec<usize>,
    current_chr: usize,
    idx: usize,
}

impl<'a> Iterator for IdxCheckIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<(usize, usize)> {
        if self.idx >= self.end_coor[self.current_chr] {
            self.current_chr += 1;
        }
        if self.current_chr < self.end_coor.len() {
            let mut pos = self.idx;
            if self.current_chr > 0 {
                pos -= self.end_coor[self.current_chr - 1];
            }
            self.idx += 1;
            Some((self.current_chr, pos))
        } else {
            None
        }
    }
}
