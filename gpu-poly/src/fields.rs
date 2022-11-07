use crate::GpuFftField;
use crate::GpuField;
use crate::GpuMul;
use ark_ff::BigInt;
use ark_ff::Field;
use ark_ff::Fp3;
use ark_ff::Fp3Config;
use ark_ff::FpConfig;
use std::ops::Add;
use std::ops::Mul;
use std::ops::MulAssign;

pub mod p18446744069414584321 {
    use super::*;
    use ark_ff_optimized::fp64;
    pub use fp64::Fp;
    pub use fp64::FpParams;
    use std::marker::PhantomData;

    impl GpuField for fp64::Fp {
        type FftField = Self;

        fn field_name() -> String {
            "p18446744069414584321_fp".to_string()
        }
    }

    impl GpuMul<fp64::Fp, fp64::Fp> for fp64::Fp {}

    impl GpuFftField for fp64::Fp {}

    const TRACE: ark_ff::BigInt<3> = BigInt!("1461501636310055817916238417282618014431694553085");

    pub struct Fq3Config;

    impl Fp3Config for Fq3Config {
        type Fp = Fp;
        const NONRESIDUE: Fp = /* =2 */ ark_ff::Fp(BigInt([8589934590]), PhantomData);
        const TWO_ADICITY: u32 = FpParams::TWO_ADICITY;
        const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = &TRACE.divide_by_2_round_down().0;
        const QUADRATIC_NONRESIDUE_TO_T: Fp3<Fq3Config> = Fp3::<Fq3Config>::new(
            /* =16140901060737761281 */
            ark_ff::Fp(BigInt([2305843009213693952]), PhantomData),
            Fp::ZERO,
            Fp::ZERO,
        );

        // NOTE: these are used for pairings which I don't need so they are left empty
        const FROBENIUS_COEFF_FP3_C1: &'static [Fp] = &[];
        const FROBENIUS_COEFF_FP3_C2: &'static [Fp] = &[];
    }

    wrap_field!(Fq3; Fp3<Fq3Config>);

    impl MulAssign<&Fp> for Fq3 {
        fn mul_assign(&mut self, rhs: &Fp) {
            self.0.mul_assign_by_base_field(rhs)
        }
    }

    impl MulAssign<Fp> for Fq3 {
        fn mul_assign(&mut self, rhs: Fp) {
            self.0.mul_assign_by_base_field(&rhs)
        }
    }

    impl Add<&Fp> for Fq3 {
        type Output = Fq3;

        fn add(self, rhs: &Fp) -> Self::Output {
            self + Fq3::from(*rhs)
        }
    }

    impl Mul<&Fp> for Fq3 {
        type Output = Fq3;

        fn mul(self, rhs: &Fp) -> Self::Output {
            let mut tmp = self;
            tmp.0.mul_assign_by_base_field(rhs);
            tmp
        }
    }

    impl Mul<Fp> for Fq3 {
        type Output = Fq3;

        fn mul(self, rhs: Fp) -> Self::Output {
            let mut tmp = self;
            tmp.0.mul_assign_by_base_field(&rhs);
            tmp
        }
    }

    impl From<Fp> for Fq3 {
        fn from(value: Fp) -> Self {
            Fq3(Fp3::<Fq3Config>::from_base_prime_field(value))
        }
    }

    impl GpuMul<Fp, Fq3> for Fq3 {}

    impl GpuField for Fq3 {
        type FftField = Fp;

        fn field_name() -> String {
            "p18446744069414584321_fq3".to_string()
        }
    }
}
