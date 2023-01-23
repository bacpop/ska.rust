
use std::{error::Error, fmt};

#[derive(Debug)]
pub struct AlnWriter<'a> {
    next_pos: usize,
    curr_chrom: usize,
    last_mapped: usize,
    last_written: usize,
    ref_seq: &'a Vec<Vec<u8>>,
    total_size: usize,
    seq_out: Vec<u8>,
    k: usize,
    half_split_len: usize
}

impl<'a> AlnWriter<'a> {
    pub fn new(ref_seq: &'a Vec<Vec<u8>>, total_size: usize, k: usize) -> Self {
        let (next_pos, curr_chrom, last_mapped, last_written) = (0, 0, 0, 0);
        let seq_out = Vec::new();
        seq_out.reserve(total_size);
        let half_split_len = (k - 1) / 2;
        Self {next_pos, curr_chrom, last_mapped, last_written, ref_seq, total_size, seq_out, k, half_split_len}
    }

    // TODO: fill to end of chromosome. Also called at the end
    fn fill_contig(&mut self) {
        seq.extend_from_slice(
            &self.seq[curr_chrom][next_pos..(next_pos + half_split_len)],
        );
        seq.extend(vec![
            b'-';
            self.seq[curr_chrom].len() - (next_pos + half_split_len) // This shouldn't overflow
        ]);
        curr_chrom += 1;
        next_pos = 0;
    }

    pub fn write_split_kmer(&mut self, mapped_pos: usize, mapped_chrom: usize, base: u8) {
        if mapped_chrom > self.curr_chrom {
            self.fill_contig();
        }
        if base != b'-' {
            if mapped_pos < self.next_pos {
                self.last_mapped = mapped_pos + self.half_split_len;
            } else {
                // Write bases
                // TODO easiest to write this out on a piece of paper
                // probably want to calculate lengths first, then write to seq
                // to make things clearer
                // TODO bases at the end of last valid match not yet written
                last_match_overhang = 
                seq.extend_from_slice(
                    &self.seq[curr_chrom][seq.len()..(last_mapped + half_split_len)],
                );
                // TODO missing bases
                seq.extend(vec![b'-'; *map_pos - next_pos - half_split_len]);
                // TODO first half of split k-mer
                seq.extend_from_slice(
                    &self.seq[curr_chrom][(*map_pos - half_split_len)..*map_pos],
                );
                // TODO base
                seq.push(*base);
                // TODO second half of split k-mer
                seq.extend_from_slice(
                    &self.seq[curr_chrom][(*map_pos + 1)..=(*map_pos + half_split_len)],
                );
                // TODO update indices
                next_pos = *map_pos + self.k;
            }
        }
    }

    pub fn get_seq(&mut self) -> Result<&'a[u8], Box<dyn self::Error>> {
        self.fill_contig();
        if self.ref_seq.len() != self.total_size {
            Result::Err(Box::new(self))
        } else {
            Result::Ok(self.seq_out.as_slice())
        }
    }
}

impl Error for &mut AlnWriter<'_> {}

impl fmt::Display for AlnWriter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Internal error: output length {} not identical to chromosome length {}",
        self.ref_seq.len(),
        self.total_size)
    }
}
