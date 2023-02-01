
use std::{error::Error, fmt};

pub struct AlnWriter<'a> {
    next_pos: usize,
    curr_chrom: usize,
    last_mapped: usize,
    last_written: usize,
    ref_seq: &'a Vec<Vec<u8>>,
    total_size: usize,
    seq_out: Vec<u8>,
    half_split_len: usize
}

// TODO two things to try (and time)
// 1) Copy the whole seq then only edit the middle base and missing bases
// 2) Parallelise this, and write out at the end
impl<'a> AlnWriter<'a> {
    pub fn new(ref_seq: &'a Vec<Vec<u8>>, total_size: usize, k: usize) -> Self {
        let (curr_chrom, last_mapped, last_written) = (0, 0, 0);
        let seq_out = vec![b'-'; total_size];
        let half_split_len = (k - 1) / 2;
        let next_pos = half_split_len;
        Self {next_pos, curr_chrom, last_mapped, last_written, ref_seq, total_size, seq_out, half_split_len}
    }

    // Fill to end of chromosome. Also called at the end
    // AC[GGC]----[AT  |   C`]
    //  1   2          3
    // 1: Last written -> last map-pos + half-k
    // 2: Last mapped + half-k
    // 3: Mapped position
    fn fill_to(&mut self, end_position: usize, contig_end: bool) {
        // bases at the end of last valid match not yet written
        if self.last_written > 0 {
            let last_match_overhang = (self.last_mapped + self.half_split_len).saturating_sub(self.last_written);
            if last_match_overhang > 0 {
                self.seq_out.extend_from_slice(
                    &self.ref_seq[self.curr_chrom][(self.last_written + 1)..(self.last_written + last_match_overhang)]
                );
            }
        }
        println!("{}", String::from_utf8(self.seq_out.clone()).unwrap());
        // 2-start of split k-mer: missing bases
        let mut missing_bases = end_position.saturating_sub(self.last_mapped + self.half_split_len);
        // This accounts for the fact that at the start of the sequence no bases have
        // been written yet
        if !contig_end && self.last_written > 0 {
            missing_bases = missing_bases.saturating_sub(self.half_split_len);
        }
        if missing_bases > 0 {
            self.seq_out.extend(vec![b'-'; missing_bases]);
        }
        println!("{}", String::from_utf8(self.seq_out.clone()).unwrap());
    }

    // Fill up to the end of the contig and update indices
    fn fill_contig(&mut self) {
        self.fill_to(self.ref_seq[self.curr_chrom].len(), true);
        self.curr_chrom += 1;
        self.next_pos = self.half_split_len;
    }

    pub fn write_split_kmer(&mut self, mapped_pos: usize, mapped_chrom: usize, base: u8) {
        println!("{} {}", mapped_pos, base);
        if mapped_chrom > self.curr_chrom {
            self.fill_contig();
        }
        if base == b'-' {
            return;
        }
        if mapped_pos < self.next_pos {
            self.last_mapped = mapped_pos;
            return;
        }

        // Write bases between last match and this one
        if mapped_pos > self.next_pos {
            self.fill_to(mapped_pos, false);
        }

        // First half of split k-mer
        self.seq_out.extend_from_slice(
            &self.ref_seq[self.curr_chrom][(mapped_pos - self.half_split_len)..mapped_pos],
        );
        // Middle base
        self.seq_out.push(base);
        // Second half of split k-mer
        self.seq_out.extend_from_slice(
            &self.ref_seq[self.curr_chrom][(mapped_pos + 1)..=(mapped_pos + self.half_split_len)],
        );
        println!("{}", String::from_utf8(self.seq_out.clone()).unwrap());

        // update indices
        self.next_pos = mapped_pos + 2 * self.half_split_len + 1;
        self.last_mapped = mapped_pos;
        self.last_written = mapped_pos + self.half_split_len;

        return;
    }

    pub fn get_seq(&'a mut self) -> Result<&'a[u8], Box<dyn self::Error + 'a>> {
        println!("{}", String::from_utf8(self.seq_out.clone()).unwrap());
        self.fill_contig();
        println!("{}", String::from_utf8(self.seq_out.clone()).unwrap());
        if self.seq_out.len() != self.total_size {
            Result::Err(Box::new(self))
        } else {
            Result::Ok(self.seq_out.as_slice())
        }
    }
}

impl Error for &mut AlnWriter<'_> {}

impl fmt::Display for AlnWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8(self.seq_out.clone()).unwrap())
    }
}

impl fmt::Debug for AlnWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Output length {} not identical to chromosome length {}",
        self.seq_out.len(),
        self.total_size)
    }
}