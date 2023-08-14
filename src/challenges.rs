use crate::constraints::VerifierChallenge;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use core::ops::Deref;
use core::ops::Index;

#[derive(Default, Clone, Debug, CanonicalDeserialize, CanonicalSerialize)]
pub struct Challenges<F: Field>(Vec<F>);

impl<F: Field> Challenges<F> {
    pub fn new(challenges: Vec<F>) -> Self {
        Self(challenges)
    }
}

impl<F: Field> Deref for Challenges<F> {
    type Target = Vec<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: Field, C: VerifierChallenge> Index<C> for Challenges<F> {
    type Output = F;

    fn index(&self, challenge: C) -> &Self::Output {
        &self.0[challenge.index()]
    }
}
