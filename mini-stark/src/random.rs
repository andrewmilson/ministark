use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use digest::Digest;
use digest::Output;
use rand_chacha::rand_core::SeedableRng;
use rand_chacha::ChaCha20Rng;

pub struct PublicCoin<D: Digest> {
    pub seed: Output<D>,
    counter: usize,
}

impl<D: Digest> PublicCoin<D> {
    pub fn new(seed: &[u8]) -> Self {
        let mut hasher = D::new();
        hasher.update(seed);
        PublicCoin {
            seed: hasher.finalize(),
            counter: 0,
        }
    }

    pub fn reseed(&mut self, item: &impl CanonicalSerialize) {
        let mut data = Vec::new();
        item.serialize_compressed(&mut data).unwrap();
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(data);
        self.seed = hasher.finalize();
        self.counter = 0;
    }

    pub fn draw<F: Field>(&mut self) -> F {
        F::rand(&mut self.draw_rng())
    }

    // TODO: make this generic
    pub fn draw_rng(&mut self) -> ChaCha20Rng {
        let mut seed: [u8; 32] = Default::default();
        seed.copy_from_slice(&self.next()[0..32]);
        ChaCha20Rng::from_seed(seed)
    }

    /// Updates the state by incrementing the counter and returns hash(seed ||
    /// counter)
    fn next(&mut self) -> Output<D> {
        self.counter += 1;
        let mut hasher = D::new();
        hasher.update(&self.seed);
        hasher.update(self.counter.to_be_bytes());
        hasher.finalize()
    }
}
