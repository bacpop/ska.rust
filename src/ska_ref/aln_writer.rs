#[derive(Clone)]
pub struct AlnWriter<'a> {
    next_pos: usize,
    curr_chrom: usize,
    last_mapped: usize,
    last_written: usize,
    chrom_offset: usize,
    ref_seq: &'a Vec<Vec<u8>>,
    seq_out: Vec<u8>,
    half_split_len: usize,
}

impl<'a> AlnWriter<'a> {
    pub fn new(ref_seq: &'a Vec<Vec<u8>>, k: usize) -> Self {
        let (curr_chrom, last_mapped, last_written, chrom_offset) = (0, 0, 0, 0);
        let total_size = ref_seq.iter().map(|x| x.len()).sum();
        let seq_out = vec![b'-'; total_size];
        let half_split_len = (k - 1) / 2;
        let next_pos = half_split_len;
        Self {
            next_pos,
            curr_chrom,
            last_mapped,
            last_written,
            chrom_offset,
            ref_seq,
            seq_out,
            half_split_len,
        }
    }

    // Fill to end of chromosome. Also called at the end
    // AC[GGC]----[AT  |   CC]
    //  1   2          3
    // 1: Last written -> last map-pos + half-k
    // 2: Last mapped + half-k
    // 3: Mapped position
    fn fill_bases(&mut self) {
        // bases at the end of last valid match not yet written
        if self.last_written > 0 {
            let last_match_overhang =
                (self.last_mapped + self.half_split_len).saturating_sub(self.last_written);
            if last_match_overhang > 0 {
                let start = self.last_written;
                let end = self.last_written + last_match_overhang;
                self.seq_out[(start + self.chrom_offset)..=(end + self.chrom_offset)]
                    .copy_from_slice(&self.ref_seq[self.curr_chrom][start..=end]);
            }
        }
    }

    // Fill up to the end of the contig and update indices
    fn fill_contig(&mut self) {
        self.fill_bases();
        let chrom_length = self.ref_seq[self.curr_chrom].len();
        self.chrom_offset += chrom_length;
        self.curr_chrom += 1;
        self.next_pos = self.half_split_len;
    }

    pub fn write_split_kmer(&mut self, mapped_pos: usize, mapped_chrom: usize, base: u8) {
        if mapped_chrom > self.curr_chrom {
            self.fill_contig();
        }

        if mapped_pos < self.next_pos {
            self.last_mapped = mapped_pos;
        } else {
            // Write bases between last match and this one
            if mapped_pos > self.next_pos {
                self.fill_bases();
            }

            // First half of split k-mer
            let start = mapped_pos - self.half_split_len;
            let end = mapped_pos;
            self.seq_out[(start + self.chrom_offset)..(end + self.chrom_offset)]
                .copy_from_slice(&self.ref_seq[self.curr_chrom][start..end]);
            // Middle base
            self.seq_out[mapped_pos + self.chrom_offset] = base;
            // Second half of split k-mer
            let start = mapped_pos + 1;
            let end = mapped_pos + self.half_split_len;
            self.seq_out[(start + self.chrom_offset)..(end + self.chrom_offset)]
                .copy_from_slice(&self.ref_seq[self.curr_chrom][start..end]);

            // update indices
            self.next_pos = mapped_pos + 2 * self.half_split_len;
            self.last_mapped = mapped_pos;
            self.last_written = mapped_pos + self.half_split_len;
        }
    }

    pub fn get_seq(&'a mut self) -> &'a [u8] {
        self.fill_contig();
        self.seq_out.as_slice()
    }
}
