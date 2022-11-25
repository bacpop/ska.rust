// Easy encoding from ASCII are the 3rd and 2nd bits
// This encodes as A: 00; C: 01; T: 10; G: 11
// EOR w/ 10          10     10     10     10
// gives rc           10     11     00     01
const LETTER_CODE: [u8; 4] = [b'A', b'C', b'T', b'G'];

pub fn generate_masks(k: usize) -> (u64, u64) {
    let half_size: usize = (k - 1) / 2;
    let mut lower_mask = 0;
    for _idx in 0..(half_size * 2) {
        lower_mask <<= 1;
        lower_mask += 1;
    }
    let upper_mask = lower_mask << (half_size * 2);
    return (lower_mask, upper_mask);
}

#[inline(always)]
pub fn encode_base(base: u8) -> u8 {
    (base >> 1) & 0x3
}

#[inline(always)]
pub fn decode_base(bitbase: u8) -> u8 {
    LETTER_CODE[bitbase as usize]
}

#[inline(always)]
pub fn rc_base(base: u8) -> u8 {
    base ^ 2
}

// Checks for N or n
#[inline(always)]
pub fn valid_base(base: u8) -> bool {
    base & 0xF != 14
}

pub fn decode_kmer(k: usize, kmer: u64, upper_mask: u64, lower_mask: u64) -> (String, String) {
    let half_k: usize = (k - 1) / 2;
    let mut upper_bits = (kmer & upper_mask) >> (half_k * 2);
    let mut upper_kmer = String::with_capacity(half_k);
    for _idx in 0..half_k {
        let base = decode_base((upper_bits & 0x3) as u8);
        upper_kmer.push(base as char);
        upper_bits = upper_bits >> 2;
    }
    upper_kmer = upper_kmer.chars().rev().collect::<String>();

    let mut lower_bits = kmer & lower_mask;
    let mut lower_kmer = String::with_capacity(half_k);
    for _idx in 0..half_k {
        let base = decode_base((lower_bits & 0x3) as u8);
        lower_kmer.push(base as char);
        lower_bits = lower_bits >> 2;
    }
    lower_kmer = lower_kmer.chars().rev().collect::<String>();
    (upper_kmer, lower_kmer)
}

// https://www.biostars.org/p/113640/
#[inline(always)]
pub fn revcomp64_v2(mut res: u64, k_size: usize) -> u64 {
    res = (res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2;
    res = (res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4;
    res = (res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8;
    res = (res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16;
    res = (res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32;
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return res >> (2 * (32 - k_size));
}

// A 	Adenine
// C 	Cytosine
// G 	Guanine
// T (or U) 	Thymine (or Uracil)
// R 	A or G
// Y 	C or T
// S 	G or C
// W 	A or T
// K 	G or T
// M 	A or C
// B 	C or G or T
// D 	A or G or T
// H 	A or C or T
// V 	A or C or G
// N 	any base
// . or - 	gap

// A + A -> A   C + A -> M   T + A -> W   G + A -> R
// A + C -> M   C + C -> C   T + C -> Y   G + C -> S
// A + G -> R   C + G -> S   T + G -> K   G + G -> G
// A + T -> W   C + T -> Y   T + T -> T   G + T -> K
// A + R -> R   C + R -> V   T + R -> D   G + R -> R
// A + Y -> H   C + Y -> Y   T + Y -> Y   G + Y -> B
// A + S -> V   C + S -> S   T + S -> B   G + S -> S
// A + W -> W   C + W -> H   T + W -> W   G + W -> D
// A + K -> D   C + K -> B   T + K -> K   G + K -> K
// A + M -> M   C + M -> M   T + M -> H   G + M -> V
// A + B -> N   C + B -> B   T + B -> B   G + B -> B
// A + D -> D   C + D -> N   T + D -> D   G + D -> D
// A + H -> H   C + H -> H   T + H -> H   G + H -> N
// A + V -> V   C + V -> V   T + V -> N   G + V -> V
// A + N -> N   C + N -> N   T + N -> N   G + N -> N

pub const IUPAC: [u8; 1024] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0-15
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16-31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32-47
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48-63
    0, b'A', b'N', b'M', b'D', 0, 0, b'R', b'H', 0, 0, b'D', 0, b'M', b'N', 0, // 64-79
    0, 0, b'R', b'V', b'W', 0, b'V', b'W', 0, b'H', 0, 0, 0, 0, 0, 0, // 80-95
    0, b'A', b'N', b'M', b'D', 0, 0, b'R', b'H', 0, 0, b'D', 0, b'M', b'N', 0, // 96-111
    0, 0, b'R', b'V', b'W', 0, b'V', b'W', 0, b'H', 0, 0, 0, 0, 0, 0, // 112-127
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128-143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 144-159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 160-175
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 176-191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 192-207
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 208-223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 224-239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 240-255
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0-15
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16-31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32-47
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48-63
    0, b'M', b'B', b'C', b'N', 0, 0, b'S', b'H', 0, 0, b'B', 0, b'M', b'N', 0, // 64-79
    0, 0, b'V', b'S', b'Y', 0, b'V', b'H', 0, b'Y', 0, 0, 0, 0, 0, 0, // 80-95
    0, b'M', b'B', b'C', b'N', 0, 0, b'S', b'H', 0, 0, b'B', 0, b'M', b'N', 0, // 96-111
    0, 0, b'V', b'S', b'Y', 0, b'V', b'H', 0, b'Y', 0, 0, 0, 0, 0, 0, // 112-127
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128-143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 144-159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 160-175
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 176-191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 192-207
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 208-223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 224-239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //240-255
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0-15
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16-31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32-47
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48-63
    0, b'W', b'B', b'Y', b'D', 0, 0, b'K', b'H', 0, 0, b'K', 0, b'H', b'N', 0, // 64-79
    0, 0, b'D', b'B', b'T', 0, b'N', b'W', 0, b'Y', 0, 0, 0, 0, 0, 0, // 80-95
    0, b'W', b'B', b'Y', b'D', 0, 0, b'K', b'H', 0, 0, b'K', 0, b'H', b'N', 0, // 96-111
    0, 0, b'D', b'B', b'T', 0, b'N', b'W', 0, b'Y', 0, 0, 0, 0, 0, 0, // 112-127
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128-143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 144-159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 160-175
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 176-191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 192-207
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 208-223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 224-239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //240-255
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0-15
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 16-31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32-47
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 48-63
    0, b'R', b'B', b'S', b'D', 0, 0, b'G', b'N', 0, 0, b'K', 0, b'V', b'N', 0, // 64-79
    0, 0, b'R', b'S', b'K', 0, b'V', b'D', 0, b'B', 0, 0, 0, 0, 0, 0, // 80-95
    0, b'R', b'B', b'S', b'D', 0, 0, b'G', b'N', 0, 0, b'K', 0, b'V', b'N', 0, // 96-111
    0, 0, b'R', b'S', b'K', 0, b'V', b'D', 0, b'B', 0, 0, 0, 0, 0, 0, // 112-127
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 128-143
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 144-159
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 160-175
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 176-191
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 192-207
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 208-223
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 224-239
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]; //240-255
