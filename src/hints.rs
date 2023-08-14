use crate::constraints::Hint;
use alloc::vec::Vec;
use ark_ff::Field;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use core::ops::Deref;
use core::ops::Index;

#[derive(Default, Debug, Clone, CanonicalDeserialize, CanonicalSerialize)]
pub struct Hints<F: Field>(Vec<F>);

impl<F: Field> Hints<F> {
    pub fn new(mut hints: Vec<(usize, F)>) -> Self {
        hints.sort();
        for [(a, _), (b, _)] in hints.array_windows() {
            assert!(a != b, "multiple hints exist at index {a}");
        }
        for (expected, (actual, _)) in hints.iter().enumerate() {
            assert!(expected == *actual, "missing hint at index {expected}");
        }
        Self(hints.into_iter().map(|(_, value)| value).collect())
    }
}

impl<F: Field> Deref for Hints<F> {
    type Target = Vec<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: Field, H: Hint> Index<H> for Hints<F> {
    type Output = F;

    fn index(&self, hint: H) -> &Self::Output {
        &self.0[hint.index()]
    }
}
