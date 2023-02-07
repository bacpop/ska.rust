#[derive(Debug)]
pub struct IdxCheck {
    end_coor: Vec<usize>,
    current_chr: usize,
}

impl IdxCheck {
    pub fn new(ref_seq: &[Vec<u8>]) -> Self {
        let mut end_coor = Vec::new();
        let current_chr = 0;

        let mut cum_pos = 0;
        for chrom in ref_seq {
            cum_pos += chrom.len();
            end_coor.push(cum_pos);
        }

        Self {
            end_coor,
            current_chr,
        }
    }

    /// Assumes idx is incremented, no skips or interesting intervals
    pub fn idx_to_coor(&mut self, idx: usize) -> (usize, usize) {
        if idx >= self.end_coor[self.current_chr] {
            self.current_chr += 1;
        }
        let mut pos = idx;
        if self.current_chr > 0 {
            pos -= self.end_coor[self.current_chr - 1];
        }
        (self.current_chr, pos)
    }
}
