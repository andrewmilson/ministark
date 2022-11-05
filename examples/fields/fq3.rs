use ark_ff::BigInt;
use ark_ff::Field;
use ark_ff::Fp3;
use ark_ff::Fp3Config;
use ark_ff::FpConfig;
use ark_ff_optimized::fp64::Fp;
use ark_ff_optimized::fp64::FpParams;
use std::marker::PhantomData;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::MulAssign;
use std::ops::Sub;
use std::ops::SubAssign;

const TRACE: ark_ff::BigInt<3> = BigInt!("1461501636310055817916238417282618014431694553085");

pub struct Fq3Config;

type Fq3 = Fp3<Fq3Config>;

impl Fp3Config for Fq3Config {
    type Fp = Fp;
    const NONRESIDUE: Fp = /* =2 */ ark_ff::Fp(BigInt([8589934590]), PhantomData);
    const TWO_ADICITY: u32 = FpParams::TWO_ADICITY;
    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = &TRACE.divide_by_2_round_down().0;
    const QUADRATIC_NONRESIDUE_TO_T: Fq3 = Fq3::new(
        /* =16140901060737761281 */ ark_ff::Fp(BigInt([2305843009213693952]), PhantomData),
        Fp::ZERO,
        Fp::ZERO,
    );

    // NOTE: these are used for pairings which I don't need so they are left empty
    const FROBENIUS_COEFF_FP3_C1: &'static [Fp] = &[];
    const FROBENIUS_COEFF_FP3_C2: &'static [Fp] = &[];
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WrappedFq3(Fq3);

impl ark_ff::Zero for WrappedFq3 {
    fn zero() -> Self {
        WrappedFq3(Fq3::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl MulAssign<Fp> for WrappedFq3 {
    fn mul_assign(&mut self, rhs: Fp) {
        self.0 *= Fq3::from_base_prime_field(rhs);
    }
}

impl Add<WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn add(self, rhs: WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 + rhs.0)
    }
}

impl Sub<WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn sub(self, rhs: WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 - rhs.0)
    }
}

impl AddAssign<WrappedFq3> for WrappedFq3 {
    fn add_assign(&mut self, rhs: WrappedFq3) {
        self.0 += rhs.0
    }
}

impl SubAssign<WrappedFq3> for WrappedFq3 {
    fn sub_assign(&mut self, rhs: WrappedFq3) {
        self.0 -= rhs.0
    }
}
