//! Keep track of current chromosome during iteration over a concatenated
//! alignment.
//!
//! Simple counter with state to avoid doing range lookups (e.g. with an
//! interval tree) at every position. Used when converting alignment to VCF.

/// Builds the contig change indices and keeps track of the current chromosome.
#[derive(Debug)]
pub struct IdxCheck {
    /// The end coordinates (last base) of each chromosome
    end_coor: Vec<usize>,
}

impl IdxCheck {
    /// Create indicies from Vec of Vec repr of a reference sequence
    pub fn new(ref_seq: &[Vec<u8>]) -> Self {
        let mut end_coor = Vec::new();

        let mut cum_pos = 0;
        for chrom in ref_seq {
            cum_pos += chrom.len();
            end_coor.push(cum_pos);
        }

        Self { end_coor }
    }

    /// Iterate over index to return (chromosome, position) tuples
    pub fn iter(&self) -> IdxCheckIter<'_> {
        IdxCheckIter {
            end_coor: &self.end_coor,
            current_chr: 0,
            idx: 0,
        }
    }
}

/// Carries state which keeps track of chromosome during iteration
///
/// Iterator as separate class so [`IdxCheck`] not modified during iteration
pub struct IdxCheckIter<'a> {
    /// Ref to end coordinates
    end_coor: &'a Vec<usize>,
    /// Current chromosome
    current_chr: usize,
    /// Current absolute position
    idx: usize,
}

impl Iterator for IdxCheckIter<'_> {
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
