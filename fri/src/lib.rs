#![feature(int_log)]

use algebra::StarkFelt;
use rand::Rng;
use std::marker::PhantomData;

pub struct Fri<E> {
    expansion_factor: usize,
    _phantom: PhantomData<E>,
}

impl<E: StarkFelt> Fri<E> {
    pub fn new(expansion_factor: usize) -> Fri<E> {
        assert!(expansion_factor.is_power_of_two(), "not a power of two");
        Fri {
            expansion_factor,
            _phantom: PhantomData,
        }
    }

    // fn commit() {}

    pub fn prove(&self, codeword: Vec<E>, rng: impl Rng) {
        assert!(codeword.len().is_power_of_two(), "length not power of two");
        let num_rounds = (codeword.len() / self.expansion_factor).ilog2();
        // let mut codewords
    }
}
