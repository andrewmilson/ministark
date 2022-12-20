use crate::constraints::Hint;
use alloc::vec::Vec;
use core::ops::Deref;
use core::ops::Index;
use gpu_poly::GpuField;

#[derive(Default)]
pub struct Hints<F: GpuField>(Vec<F>);

impl<F: GpuField> Hints<F> {
    pub fn new(mut hints: Vec<(usize, F)>) -> Self {
        hints.sort();
        for [(a, _), (b, _)] in hints.array_windows() {
            assert!(a != b, "multiple hints exist at index {a}");
        }
        for (expected, (actual, _)) in (0..hints.len()).zip(&hints) {
            assert!(expected == *actual, "missing hint at index {expected}")
        }
        Hints(hints.into_iter().map(|(_, value)| value).collect())
    }
}

impl<F: GpuField> Deref for Hints<F> {
    type Target = Vec<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: GpuField, H: Hint> Index<H> for Hints<F> {
    type Output = F;

    fn index(&self, hint: H) -> &Self::Output {
        &self.0[hint.index()]
    }
}
