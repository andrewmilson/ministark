use alloc::vec::Vec;
use ark_ff::Field;
use digest::Digest;
use digest::Output;
use rand::Rng;
use rand::RngCore;
use std::collections::BTreeSet;
use std::fmt::Debug;
use std::marker::PhantomData;

/// PublicCoin trait adapted from Winterfell
pub trait PublicCoin: Sync + Debug {
    type Digest: Digest;
    type Field: Field;

    fn new(digest: Output<Self::Digest>) -> Self;

    fn reseed_with_hash(&mut self, val: &Output<Self::Digest>);

    fn reseed_with_field_element(&mut self, val: &Self::Field);

    fn reseed_with_field_elements(&mut self, vals: &[Self::Field]) {
        for val in vals {
            self.reseed_with_field_element(val);
        }
    }

    fn reseed_with_int(&mut self, val: u64);

    fn draw(&mut self) -> Self::Field;

    /// Draws a maximum of n unique queries in the range `[0, domain_size)`
    fn draw_queries(&mut self, max_n: usize, domain_size: usize) -> BTreeSet<usize>;

    fn verify_proof_of_work(&self, proof_of_work_bits: u8, nonce: u64) -> bool;
}

pub struct PublicCoinImpl<D: Digest, F: Field> {
    pub seed: Output<D>,
    counter: usize,
    bytes: Vec<u8>,
    _phantom: PhantomData<F>,
}

impl<D: Digest, F: Field> Debug for PublicCoinImpl<D, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PublicCoinImpl")
            .field("seed", &self.seed)
            .field("counter", &self.counter)
            .field("bytes", &self.bytes)
            .field("_phantom", &self._phantom)
            .finish()
    }
}

impl<D: Digest, F: Field> PublicCoinImpl<D, F> {
    /// Updates the state by incrementing the counter and returns hash(seed ||
    /// counter)
    fn gen_next(&mut self) -> Output<D> {
        self.counter += 1;
        self.bytes = Vec::new();
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(self.counter.to_be_bytes());
        hasher.finalize()
    }
}

impl<D: Digest, F: Field> PublicCoin for PublicCoinImpl<D, F> {
    type Digest = D;
    type Field = F;

    fn new(digest: Output<D>) -> Self {
        Self {
            seed: digest,
            counter: 0,
            bytes: Vec::new(),
            _phantom: PhantomData,
        }
    }

    fn reseed_with_hash(&mut self, val: &Output<D>) {
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(val);
        self.seed = hasher.finalize();
        self.counter = 0;
        self.bytes = Vec::new();
    }

    fn reseed_with_field_element(&mut self, val: &Self::Field) {
        let mut val_bytes = Vec::new();
        val.serialize_uncompressed(&mut val_bytes).unwrap();
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(val_bytes);
        self.seed = hasher.finalize();
        self.counter = 0;
        self.bytes = Vec::new();
    }

    fn reseed_with_int(&mut self, val: u64) {
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(val.to_be_bytes());
        self.seed = hasher.finalize();
        self.counter = 0;
        self.bytes = Vec::new();
    }

    fn verify_proof_of_work(&self, proof_of_work_bits: u8, nonce: u64) -> bool {
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(nonce.to_be_bytes());
        leading_zeros(&hasher.finalize()) >= u32::from(proof_of_work_bits)
    }

    fn draw(&mut self) -> F {
        F::rand(self)
    }

    fn draw_queries(&mut self, max_n: usize, domain_size: usize) -> BTreeSet<usize> {
        (0..max_n).map(|_| self.gen_range(0..domain_size)).collect()
    }
}

impl<D: Digest, F: Field> Iterator for PublicCoinImpl<D, F> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.bytes.is_empty() {
            self.bytes = self.gen_next().to_vec();
        }
        self.bytes.pop()
    }
}

impl<D: Digest, F: Field> RngCore for PublicCoinImpl<D, F> {
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
