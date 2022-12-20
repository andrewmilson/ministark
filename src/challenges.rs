use crate::constraints::VerifierChallenge;
use alloc::vec::Vec;
use ark_std::rand::Rng;
use core::ops::Deref;
use core::ops::Index;
use gpu_poly::GpuField;

#[derive(Default)]
pub struct Challenges<F: GpuField>(Vec<F>);

impl<F: GpuField> Challenges<F> {
    pub fn new<R: Rng + ?Sized>(rng: &mut R, num_challenges: usize) -> Self {
        Challenges((0..num_challenges).map(|_| F::rand(rng)).collect())
    }
}

impl<F: GpuField> Deref for Challenges<F> {
    type Target = Vec<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: GpuField, C: VerifierChallenge> Index<C> for Challenges<F> {
    type Output = F;

    fn index(&self, challenge: C) -> &Self::Output {
        &self.0[challenge.index()]
    }
}
