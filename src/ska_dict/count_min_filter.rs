use ahash::RandomState;

#[derive(Debug, Clone)]
pub struct CountMin {
    width: usize,
    height: usize,
    hash_factory: Vec<RandomState>,
    mask: u64,
    counts: Vec<u16>,
    min_count: u16,
}

impl CountMin {
    pub fn empty(width: usize, height: usize, min_count: u16) -> Self {
        // Consistent with consts used, but ensures a power of two
        let width_bits: usize = f64::floor(f64::log2(width as f64)) as usize;
        let width = 1 << (width_bits + 1);
        let mask = width as u64 - 1;

        // Reserve for these gets call by the vec! macro used in init
        let hash_factory = Vec::new();
        let counts = Vec::new();

        Self {
            width,
            height,
            hash_factory,
            mask,
            counts,
            min_count,
        }
    }

    pub fn is_init(&self) -> bool {
        self.counts.len() > 0
    }

    pub fn init(&mut self) {
        if self.counts.len() == 0 {
            self.hash_factory = vec![RandomState::new(); self.height];
            self.counts = vec![0; self.width * self.height];
        }
    }

    pub fn reset(&mut self) {
        self.counts.fill(0);
    }

    pub fn filter(&mut self, kmer: u64) -> bool {
        let mut count = 0;
        for hash_it in self.hash_factory.iter().enumerate() {
            let (row_idx, hash) = hash_it;
            let table_idx = row_idx * self.width + ((hash.hash_one(kmer) & self.mask) as usize);
            if self.counts[table_idx] < u16::MAX {
                self.counts[table_idx] += 1;
            }
            count = u16::min(count, self.counts[table_idx])
        }
        return count >= self.min_count;
    }
}
