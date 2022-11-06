use ark_ff::BigInt;
use ark_ff::Field;
use ark_ff::Fp3;
use ark_ff::Fp3Config;
use ark_ff::FpConfig;
use ark_ff::One;
use ark_ff::SqrtPrecomputation;
use ark_ff::Zero;
use ark_ff_optimized::fp64::Fp;
use ark_ff_optimized::fp64::FpParams;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalDeserializeWithFlags;
use ark_serialize::CanonicalSerialize;
use ark_serialize::CanonicalSerializeWithFlags;
use ark_serialize::Valid;
use gpu_poly::GpuField;
use ministark::wrap_field;
use std::cmp::Ordering;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::DivAssign;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Neg;
use std::ops::Sub;
use std::ops::SubAssign;
use zeroize::Zeroize;

const TRACE: ark_ff::BigInt<3> = BigInt!("1461501636310055817916238417282618014431694553085");

pub struct Fq3Config;

impl Fp3Config for Fq3Config {
    type Fp = Fp;
    const NONRESIDUE: Fp = /* =2 */ ark_ff::Fp(BigInt([8589934590]), PhantomData);
    const TWO_ADICITY: u32 = FpParams::TWO_ADICITY;
    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = &TRACE.divide_by_2_round_down().0;
    const QUADRATIC_NONRESIDUE_TO_T: Fp3<Fq3Config> = Fp3::<Fq3Config>::new(
        /* =16140901060737761281 */ ark_ff::Fp(BigInt([2305843009213693952]), PhantomData),
        Fp::ZERO,
        Fp::ZERO,
    );

    // NOTE: these are used for pairings which I don't need so they are left empty
    const FROBENIUS_COEFF_FP3_C1: &'static [Fp] = &[];
    const FROBENIUS_COEFF_FP3_C2: &'static [Fp] = &[];
}

wrap_field!(Fq3; Fp3<Fq3Config>);

impl MulAssign<Fp> for Fq3 {
    fn mul_assign(&mut self, rhs: Fp) {
        self.0.mul_assign_by_base_field(&rhs)
    }
}

impl From<Fp> for Fq3 {
    fn from(value: Fp) -> Self {
        Fq3(Fp3::<Fq3Config>::from_base_prime_field(value))
    }
}

impl GpuField for Fq3 {
    type FftField = Fp;

    fn field_name() -> String {
        todo!()
    }
}
