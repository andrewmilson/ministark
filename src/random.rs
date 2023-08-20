use crate::hash::Digest;
use crate::hash::ElementHashFn;
use crate::hash::HashFn;
use alloc::vec::Vec;
use ark_ff::Field;
use rand::Rng;
use rand::RngCore;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::BTreeSet;
use std::fmt::Debug;
use std::marker::PhantomData;

// TODO: alternative approach
// trait Seedable<T>: Sync + Debug {
//     fn reseed(&mut self, v: &T);
// }
// trait PublicCoin:
//     Seedable<Self::Field> + Seedable<&[Self::Field]> +
// Seedable<FriRemainder<Self::Field>> + Seedable {
//     type Field: Field;
// }
// Seedable<Self::Fp> + Seedable<Self::Fp> + Seedable<Self::Fq> +
// Seedable<FriRemainder<Self::Fq>>

/// `PublicCoin` trait adapted from Winterfell
pub trait PublicCoin: Sized + Send + Sync + Debug {
    type Digest: Digest;
    type Field: Field;

    fn new(digest: Self::Digest) -> Self;

    fn reseed_with_digest(&mut self, val: &Self::Digest);

    fn reseed_with_field_elements(&mut self, vals: &[Self::Field]);

    fn reseed_with_field_element_vector(&mut self, vector: &[Self::Field]) {
        self.reseed_with_field_elements(vector);
    }

    fn reseed_with_int(&mut self, val: u64);

    fn draw(&mut self) -> Self::Field;

    /// Draws a maximum of n unique queries in the range `[0, domain_size)`
    fn draw_queries(&mut self, max_n: usize, domain_size: usize) -> BTreeSet<usize>;

    fn grind_proof_of_work(&self, proof_of_work_bits: u8) -> Option<u64> {
        #[cfg(not(feature = "parallel"))]
        return (1..u64::MAX).find(|&nonce| self.verify_proof_of_work(proof_of_work_bits, nonce));
        #[cfg(feature = "parallel")]
        return (1..u64::MAX)
            .into_par_iter()
            .find_any(|&nonce| self.verify_proof_of_work(proof_of_work_bits, nonce));
    }

    fn verify_proof_of_work(&self, proof_of_work_bits: u8, nonce: u64) -> bool;

    fn security_level_bits() -> u32;
}

pub struct PublicCoinImpl<F: Field, H: HashFn> {
    pub seed: H::Digest,
    counter: u64,
    bytes: Vec<u8>,
    _phantom: PhantomData<F>,
}

impl<F: Field, H: ElementHashFn<F>> PublicCoinImpl<F, H> {
    fn reseed_with_field_element(&mut self, val: &F) {
        let val_digest = H::hash_elements([*val]);
        self.seed = H::merge(&self.seed, &val_digest);
        self.counter = 0;
        self.bytes = Vec::new();
    }
}

impl<F: Field, H: HashFn> Debug for PublicCoinImpl<F, H> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PublicCoinImpl")
            .field("seed", &self.seed)
            .field("counter", &self.counter)
            .field("bytes", &self.bytes)
            .finish()
    }
}

impl<F: Field, H: HashFn> PublicCoinImpl<F, H> {
    /// Updates the state by incrementing the counter and returns hash(seed ||
    /// counter)
    fn gen_next(&mut self) -> H::Digest {
        self.counter += 1;
        self.bytes = Vec::new();
        H::merge_with_int(&self.seed, self.counter)
    }
}

impl<F: Field, H: ElementHashFn<F>> PublicCoin for PublicCoinImpl<F, H> {
    type Digest = H::Digest;
    type Field = F;

    fn new(digest: H::Digest) -> Self {
        Self {
            seed: digest,
            counter: 0,
            bytes: Vec::new(),
            _phantom: PhantomData,
        }
    }

    fn reseed_with_digest(&mut self, val: &H::Digest) {
        self.seed = H::merge(&self.seed, val);
        self.counter = 0;
        self.bytes = Vec::new();
    }

    fn reseed_with_field_elements(&mut self, vals: &[Self::Field]) {
        for val in vals {
            self.reseed_with_field_element(val);
        }
    }

    fn reseed_with_int(&mut self, val: u64) {
        self.seed = H::merge_with_int(&self.seed, val);
        self.counter = 0;
        self.bytes = Vec::new();
    }

    fn verify_proof_of_work(&self, proof_of_work_bits: u8, nonce: u64) -> bool {
        let digest = H::merge_with_int(&self.seed, nonce);
        leading_zeros(&digest.as_bytes()) >= u32::from(proof_of_work_bits)
    }

    fn draw(&mut self) -> F {
        F::rand(self)
    }

    fn draw_queries(&mut self, max_n: usize, domain_size: usize) -> BTreeSet<usize> {
        (0..max_n).map(|_| self.gen_range(0..domain_size)).collect()
    }

    fn security_level_bits() -> u32 {
        H::COLLISION_RESISTANCE
    }
}

impl<F: Field, H: HashFn> Iterator for PublicCoinImpl<F, H> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.bytes.is_empty() {
            self.bytes = self.gen_next().as_bytes().to_vec();
        }
        self.bytes.pop()
    }
}

impl<F: Field, H: HashFn> RngCore for PublicCoinImpl<F, H> {
    fn next_u32(&mut self) -> u32 {
        let mut bytes = [0; 4];
        self.fill_bytes(&mut bytes);
        u32::from_be_bytes(bytes)
    }

    fn next_u64(&mut self) -> u64 {
        let mut bytes = [0; 8];
        self.fill_bytes(&mut bytes);
        u64::from_be_bytes(bytes)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        dest.iter_mut().for_each(|v| *v = self.next().unwrap());
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

pub fn leading_zeros(hash: &[u8]) -> u32 {
    let mut zeros = 0;
    for byte in hash {
        let leading_zeros = byte.leading_zeros();
        zeros += leading_zeros;

        if leading_zeros != 8 {
            break;
        }
    }
    zeros
}

pub fn draw_multiple<P: PublicCoin>(public_coin: &mut P, n: usize) -> Vec<P::Field> {
    (0..n).map(|_| public_coin.draw()).collect()
}
