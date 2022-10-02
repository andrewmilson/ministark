// TODO: figure this shit out

use crate::plan::Fft;
use crate::plan::Planner;
use crate::GpuField;
use ark_poly::domain::radix2::Elements;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ark_serialize::Valid;
use std::fmt::Debug;
use std::hash::Hash;

pub struct GpuDomain<F: GpuField> {
    /// GPU compute object.
    fft: Fft<F>,
    /// Radix2 domain.
    pub base_domain: Radix2EvaluationDomain<F>,
}

impl<F: GpuField> Debug for GpuDomain<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        unreachable!()
    }
}

impl<F: GpuField> Hash for GpuDomain<F> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.base_domain.hash(state);
    }
}

impl<F: GpuField> PartialEq for GpuDomain<F> {
    fn eq(&self, other: &Self) -> bool {
        self.base_domain == other.base_domain
    }
}

impl<F: GpuField> Eq for GpuDomain<F> {}

impl<F: GpuField> CanonicalSerialize for &'static GpuDomain<F> {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        unreachable!()
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        unreachable!()
    }
}

impl<F: GpuField> Valid for &'static GpuDomain<F> {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        unreachable!()
    }
}

impl<F: GpuField> CanonicalDeserialize for &'static GpuDomain<F> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        _reader: R,
        _compress: ark_serialize::Compress,
        _validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        unreachable!()
    }
}

impl<F: GpuField> EvaluationDomain<F> for &'static GpuDomain<F> {
    type Elements = Elements<F>;

    fn new(num_coeffs: usize) -> Option<Self> {
        let base_domain = Radix2EvaluationDomain::<F>::new(num_coeffs)?;
        let planner = Planner::default();
        Some(Box::leak(Box::new(GpuDomain {
            fft: planner.plan_fft(base_domain.size()),
            base_domain,
        })))
    }

    fn get_coset(&self, offset: F) -> Option<Self> {
        self
    }

    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        todo!()
    }

    fn size(&self) -> usize {
        todo!()
    }

    fn log_size_of_group(&self) -> u64 {
        todo!()
    }

    fn size_inv(&self) -> F {
        todo!()
    }

    fn group_gen(&self) -> F {
        todo!()
    }

    fn group_gen_inv(&self) -> F {
        todo!()
    }

    fn coset_offset(&self) -> F {
        todo!()
    }

    fn coset_offset_inv(&self) -> F {
        todo!()
    }

    fn coset_offset_pow_size(&self) -> F {
        todo!()
    }

    fn fft_in_place<T: ark_poly::domain::DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        todo!()
    }

    fn ifft_in_place<T: ark_poly::domain::DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        todo!()
    }

    fn elements(&self) -> Elements<F> {
        self.base_domain.elements()
    }
}
