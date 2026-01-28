//! Encode and decode DNA bases into bit representation.
//!
//! The main functions are those to encode and decode bases, and get their
//! reverse complement. Each bases uses two bits, and are packed into 64-bit
//! integers.
//!
//! Easy encoding from ASCII are the 3rd and 2nd bits (works with both
//! upper and lowercase)
//!
//! ```none
//! This encodes as A: 00; C: 01; T: 10; G: 11
//!
//! EOR w/ 10          10     10     10     10
//! gives rc           10     11     00     01
//! ```
//! (Same as used in GATB library)
//!
//! Ns/n are checked cheaply in a similar way
//!
//! There are also lookup tables to support ambiguity using IUPAC codes.

#[cfg(not(target_arch = "wasm32"))]
use ahash::RandomState;

use num_traits::{PrimInt, Unsigned};
#[cfg(not(target_arch = "wasm32"))]
use std::hash::{BuildHasher, Hasher};

/// Table from bits 0-3 to ASCII (use [`decode_base()`] not this table).
const LETTER_CODE: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Encode an ASCII char to bits 0-3.
#[inline(always)]
pub fn encode_base(base: u8) -> u8 {
    (base >> 1) & 0x3
}

/// Decode bits 0-3 to ASCII.
#[inline(always)]
pub fn decode_base(bitbase: u8) -> u8 {
    LETTER_CODE[bitbase as usize]
}

/// Reverse complement an encoded base.
#[inline(always)]
pub fn rc_base(base: u8) -> u8 {
    base ^ 2
}

/// Checks for N or n with ASCII input.
#[inline(always)]
pub fn valid_base(base: u8) -> bool {
    base & 0xF != 14
}

/// Checks for A, C, G, T/U or gap with ASCII input
#[inline(always)]
pub fn is_ambiguous(mut base: u8) -> bool {
    base |= 0x20; // to lower
    !matches!(base, b'a' | b'c' | b'g' | b't' | b'u' | b'-')
}

/// Convert an ASCII base into a probability vector
/// [p(A), p(C), p(T), p(G)]
pub fn base_to_prob(base: u8) -> [f64; 4] {
    match base {
        //      (A    C    T    G  )
        b'A' => [1.0, 0.0, 0.0, 0.0],
        b'C' => [0.0, 1.0, 0.0, 0.0],
        b'G' => [0.0, 0.0, 0.0, 1.0],
        b'T' | b'U' => [0.0, 0.0, 1.0, 0.0],
        b'R' => [0.5, 0.0, 0.0, 0.5],
        b'Y' => [0.0, 0.5, 0.5, 0.0],
        b'S' => [0.0, 0.5, 0.0, 0.5],
        b'W' => [0.5, 0.0, 0.5, 0.0],
        b'K' => [0.0, 0.0, 0.5, 0.5],
        b'M' => [0.5, 0.5, 0.0, 0.0],
        b'B' => [0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        b'D' => [1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0 / 3.0],
        b'H' => [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
        b'V' => [1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0 / 3.0],
        b'N' => [0.0, 0.0, 0.0, 0.0], // This gives more expected behaviour than 0.25 in each entry
        _ => [0.0, 0.0, 0.0, 0.0],
    }
}

/// Trait to support both `u64` and `u128` representation of split k-mers
pub trait UInt<'a>:
    PrimInt
    + Unsigned
    + std::fmt::Display
    + std::fmt::Debug
    + std::marker::Send
    + std::marker::Sync
    + std::hash::Hash
    + std::ops::Shl<usize, Output = Self>
    + std::ops::BitAnd
    + std::ops::ShrAssign<i32>
    + std::ops::ShlAssign<i32>
    + std::ops::BitOrAssign
    + std::borrow::BorrowMut<Self>
    + serde::Serialize
    + serde::Deserialize<'a>
{
    /// Reverse complement of an encoded and packed split k-mer (64-bits).
    ///
    /// Neat trick from <https://www.biostars.org/p/113640/>.
    ///
    /// In the [`crate::ska_dict::split_kmer::SplitKmer`] class reverse complement is given by a rolling
    /// method, but this is used when an entire split k-mer needs to be built
    /// e.g. on construction or after skipping an N.
    fn rev_comp(self, k_size: usize) -> Self;
    /// Get the lowest (furthest right) base, encoded
    fn lsb_u8(self) -> u8;
    /// Convert to u8 as primitive
    fn as_u8(self) -> u8;
    /// Generate bit masks which can be applied to the packed k-mer representation
    /// too extract the upper and lower parts of the split k-mer (as bits).
    fn generate_masks(k: usize) -> (Self, Self);
    /// Generate a mask for the skalo algorithm
    fn skalo_mask(k: usize) -> Self;
    /// Encodes a kmer from a string to UInt
    fn encode_kmer(kmer: &[u8]) -> Self {
        kmer.iter().fold(Self::zero_init(), |result, nt| {
            (result << 2) | (Self::from_encoded_base(encode_base(*nt)))
        })
    }
    /// Encodes a kmer in str form to UInt
    fn encode_kmer_str(kmer: &str) -> Self {
        Self::encode_kmer(kmer.as_bytes())
    }
    /// Combines two kmers togethers
    fn combine_kmers(encoded_kmer1: Self, encoded_kmer2: Self) -> Self {
        let last_nucleotide_mask: Self = Self::from_encoded_base(0b11); // Mask for 2 bits

        // shift the first k-mer left by 2 bits to make space for the new nucleotide
        let shifted_kmer1 = encoded_kmer1 << 2;

        // extract the last nucleotide from the second k-mer
        let last_nucleotide = encoded_kmer2 & last_nucleotide_mask;

        // combine the two k-mers into a (k+1)-mer encoding
        shifted_kmer1 | last_nucleotide
    }
    /// Get last nucleotides from a kmer
    fn get_last_nucl(encoded_kmer: Self) -> char {
        // mask the last 2 bits to get the encoded nucleotide
        let last_bits = Self::as_u8(encoded_kmer & Self::from_encoded_base(0b11));
        // decode the nucleotide based on the 2-bit pattern
        decode_base(last_bits) as char
    }
    /// Decodes kmer string for use in skalo code
    fn skalo_decode_kmer(encoded: Self, k: usize) -> String {
        let mut kmer = String::with_capacity(k);

        let mask: Self = Self::skalo_mask(k);
        let mut value = encoded & mask;

        for _ in 0..k {
            let nucleotide =
                decode_base(Self::as_u8(value & Self::from_encoded_base(0b11))) as char;
            kmer.insert(0, nucleotide);
            value >>= 2;
        }
        kmer
    }
    /// Set to zero
    fn zero_init() -> Self {
        Self::from_encoded_base(0)
    }
    /// Convert from u8, encoded bases 0-3
    fn from_encoded_base(encoded_base: u8) -> Self;
    /// Number of bits in the representation
    fn n_bits() -> u32;
    /// Generate a `u64` hash
    #[cfg(not(target_arch = "wasm32"))]
    fn hash_val(self, hash_fn: &RandomState) -> u64;
}

impl UInt<'_> for u64 {
    #[inline(always)]
    fn rev_comp(mut self, k_size: usize) -> Self {
        // This part reverses the bases by shuffling them using an on/off pattern
        // of bits
        self = (self >> 2 & 0x3333333333333333) | (self & 0x3333333333333333) << 2;
        self = (self >> 4 & 0x0F0F0F0F0F0F0F0F) | (self & 0x0F0F0F0F0F0F0F0F) << 4;
        self = (self >> 8 & 0x00FF00FF00FF00FF) | (self & 0x00FF00FF00FF00FF) << 8;
        self = (self >> 16 & 0x0000FFFF0000FFFF) | (self & 0x0000FFFF0000FFFF) << 16;
        self = (self >> 32 & 0x00000000FFFFFFFF) | (self & 0x00000000FFFFFFFF) << 32;
        // This reverse complements
        self ^= 0xAAAAAAAAAAAAAAAA;

        // Shifts so LSB is at the bottom
        self >> (2 * (32 - k_size))
    }

    #[inline(always)]
    fn lsb_u8(self) -> u8 {
        (self & 0x3) as u8
    }

    #[inline(always)]
    fn as_u8(self) -> u8 {
        self as u8
    }

    #[inline(always)]
    fn generate_masks(k: usize) -> (Self, Self) {
        let half_size: usize = (k - 1) / 2;
        let lower_mask: Self = (1 << (half_size * 2)) - 1;
        let upper_mask: Self = lower_mask << (half_size * 2);
        (lower_mask, upper_mask)
    }

    #[inline(always)]
    fn skalo_mask(k: usize) -> Self {
        (1 << (k * 2)) - 1
    }

    #[inline(always)]
    fn from_encoded_base(encoded_base: u8) -> Self {
        encoded_base as Self
    }

    #[inline(always)]
    fn n_bits() -> u32 {
        Self::BITS
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[inline(always)]
    fn hash_val(self, hash_fn: &RandomState) -> u64 {
        let mut hasher = hash_fn.build_hasher();
        hasher.write_u64(self);
        hasher.finish()
    }
}

impl UInt<'_> for u128 {
    #[inline(always)]
    fn rev_comp(mut self, k_size: usize) -> Self {
        // This part reverses the bases by shuffling them using an on/off pattern
        // of bits
        self = (self >> 2 & 0x33333333333333333333333333333333)
            | (self & 0x33333333333333333333333333333333) << 2;
        self = (self >> 4 & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F)
            | (self & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F) << 4;
        self = (self >> 8 & 0x00FF00FF00FF00FF00FF00FF00FF00FF)
            | (self & 0x00FF00FF00FF00FF00FF00FF00FF00FF) << 8;
        self = (self >> 16 & 0x0000FFFF0000FFFF0000FFFF0000FFFF)
            | (self & 0x0000FFFF0000FFFF0000FFFF0000FFFF) << 16;
        self = (self >> 32 & 0x00000000FFFFFFFF00000000FFFFFFFF)
            | (self & 0x00000000FFFFFFFF00000000FFFFFFFF) << 32;
        self = (self >> 64 & 0x0000000000000000FFFFFFFFFFFFFFFF)
            | (self & 0x0000000000000000FFFFFFFFFFFFFFFF) << 64;
        // This reverse complements
        self ^= 0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA;

        // Shifts so LSB is at the bottom
        self >> (2 * (64 - k_size))
    }

    #[inline(always)]
    fn lsb_u8(self) -> u8 {
        (self & 0x3) as u8
    }

    #[inline(always)]
    fn as_u8(self) -> u8 {
        self as u8
    }

    #[inline(always)]
    fn generate_masks(k: usize) -> (Self, Self) {
        let half_size: usize = (k - 1) / 2;
        let lower_mask: Self = (1 << (half_size * 2)) - 1;
        let upper_mask: Self = lower_mask << (half_size * 2);
        (lower_mask, upper_mask)
    }

    #[inline(always)]
    fn skalo_mask(k: usize) -> Self {
        (1 << (k * 2)) - 1
    }

    #[inline(always)]
    fn from_encoded_base(encoded_base: u8) -> Self {
        encoded_base as Self
    }

    #[inline(always)]
    fn n_bits() -> u32 {
        Self::BITS
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[inline(always)]
    fn hash_val(self, hash_fn: &RandomState) -> u64 {
        let mut hasher = hash_fn.build_hasher();
        hasher.write_u128(self);
        hasher.finish()
    }
}

/// Decodes an encoded and packed split k-mer (64-bits) into strings for upper
/// and lower parts.
pub fn decode_kmer<IntT>(
    k: usize,
    kmer: IntT,
    upper_mask: IntT,
    lower_mask: IntT,
) -> (String, String)
where
    IntT: for<'a> UInt<'a>,
{
    let half_k: usize = (k - 1) / 2;
    let mut upper_bits = (kmer & upper_mask) >> (half_k * 2);
    let mut upper_kmer = String::with_capacity(half_k);
    for _idx in 0..half_k {
        let base = decode_base(upper_bits.lsb_u8());
        upper_kmer.push(base as char);
        upper_bits >>= 2;
    }
    upper_kmer = upper_kmer.chars().rev().collect::<String>();

    let mut lower_bits = kmer & lower_mask;
    let mut lower_kmer = String::with_capacity(half_k);
    for _idx in 0..half_k {
        let base = decode_base(lower_bits.lsb_u8());
        lower_kmer.push(base as char);
        lower_bits >>= 2;
    }
    lower_kmer = lower_kmer.chars().rev().collect::<String>();
    (upper_kmer, lower_kmer)
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

/// Lookup table to return an ambiguity code by adding a new base to an existing code.
///
/// Table is indexed as `[existing_base, new_base]` where `new_base` is two-bit encoded
/// and `existing_base` is ASCII/`u8`. Returns ASCII/`u8` IUPAC code.
///
/// Table is flattened and row-major i.e. strides are `(1, 256)`.
///
/// # Examples
///
/// ```
/// use ska::ska_dict::bit_encoding::{encode_base, IUPAC};
///
/// // A + Y -> H
/// let new_base = encode_base(b'A');
/// let existing_base = b'Y';
/// // Returns b'H'
/// let updated_base = IUPAC[new_base as usize * 256 + existing_base as usize];
/// ```
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

// IUPAC UInt
// A -> T
// C -> G
// G -> C
// T -> A
// U -> A
// R -> Y
// Y -> R
// S -> S
// W -> W
// K -> M
// M -> K
// B -> V
// D -> H
// H -> D
// V -> B
// N -> N
// - -> -

/// Lookup table which gives reverse complement of a single IUPAC code (ASCII/`u8`).
pub const RC_IUPAC: [u8; 256] = [
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 0-15
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 16-31
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 32-47
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 48-63
    b'-', b'T', b'V', b'G', b'H', b'-', b'-', b'C', b'D', b'-', b'-', b'M', b'-', b'K', b'N',
    b'-', // 64-79
    b'-', b'-', b'Y', b'S', b'A', b'-', b'B', b'W', b'-', b'R', b'-', b'-', b'-', b'-', b'-',
    b'-', // 80-95
    b'-', b'T', b'V', b'G', b'H', b'-', b'-', b'C', b'D', b'-', b'-', b'M', b'-', b'K', b'N',
    b'-', // 96-111
    b'-', b'-', b'Y', b'S', b'A', b'-', b'B', b'W', b'-', b'R', b'-', b'-', b'-', b'-', b'-',
    b'-', // 112-127
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 128-143
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 144-159
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 160-175
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 176-191
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 192-207
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 208-223
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 224-239
    b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-', b'-',
    b'-', // 240-255
];

#[cfg(test)]
mod tests {
    use super::*;

    use pretty_assertions::assert_eq;

    fn overlap(b1: &[f64; 4], b2: &[f64; 4]) -> f64 {
        b1.iter().zip(b2).map(|(p1, p2)| *p1 * p2).sum()
    }

    #[test]
    fn test_base_to_prob() {
        let a = base_to_prob(b'A');
        let c = base_to_prob(b'C');
        let g = base_to_prob(b'G');
        let t = base_to_prob(b'T');
        let u = base_to_prob(b'U');
        let r = base_to_prob(b'R');
        let y = base_to_prob(b'Y');
        let s = base_to_prob(b'S');
        let w = base_to_prob(b'W');
        let k = base_to_prob(b'K');
        let m = base_to_prob(b'M');
        let b = base_to_prob(b'B');
        let d = base_to_prob(b'D');
        let h = base_to_prob(b'H');
        let v = base_to_prob(b'V');
        let n = base_to_prob(b'N');
        let empty = base_to_prob(b'-');

        assert_eq!(overlap(&a, &a), 1.0);
        assert_eq!(overlap(&a, &c), 0.0);
        assert_eq!(overlap(&t, &u), 1.0);
        assert_eq!(overlap(&g, &u), 0.0);

        assert_eq!(overlap(&r, &a), 0.5);
        assert_eq!(overlap(&r, &y), 0.0);
        assert_eq!(overlap(&s, &g), 0.5);
        assert_eq!(overlap(&w, &w), 0.5);
        assert_eq!(overlap(&m, &y), 0.25);

        assert_eq!(overlap(&k, &b), 1.0 / 3.0);
        assert_eq!(overlap(&d, &h), 2.0 / 9.0);
        assert_eq!(overlap(&v, &n), 0.0);

        assert_eq!(overlap(&n, &empty), 0.0);
    }
}
