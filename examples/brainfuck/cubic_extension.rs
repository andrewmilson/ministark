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

#[derive(Default, Clone, Copy, Debug, PartialEq, Eq, Hash, Zeroize)]
pub struct WrappedFq3(Fq3);

impl Display for WrappedFq3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl ark_ff::UniformRand for WrappedFq3 {
    fn rand<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
        WrappedFq3(Fq3::rand(rng))
    }
}

impl Field for WrappedFq3 {
    type BasePrimeField = <Fq3 as Field>::BasePrimeField;
    type BasePrimeFieldIter = <Fq3 as Field>::BasePrimeFieldIter;
    const SQRT_PRECOMP: Option<SqrtPrecomputation<Self>> =
        Some(SqrtPrecomputation::TonelliShanks {
            two_adicity: Fq3Config::TWO_ADICITY,
            quadratic_nonresidue_to_trace: WrappedFq3(Fp3Config::QUADRATIC_NONRESIDUE_TO_T),
            trace_of_modulus_minus_one_div_two: Fq3Config::TRACE_MINUS_ONE_DIV_TWO,
        });

    const ZERO: Self = WrappedFq3(Fq3::ZERO);
    const ONE: Self = WrappedFq3(Fq3::ONE);

    fn extension_degree() -> u64 {
        Fq3::extension_degree()
    }

    fn to_base_prime_field_elements(&self) -> Self::BasePrimeFieldIter {
        self.0.to_base_prime_field_elements()
    }

    fn from_base_prime_field_elems(elems: &[Self::BasePrimeField]) -> Option<Self> {
        Fq3::from_base_prime_field_elems(elems).map(WrappedFq3)
    }

    fn from_base_prime_field(elem: Self::BasePrimeField) -> Self {
        WrappedFq3(Fq3::from_base_prime_field(elem))
    }

    fn double(&self) -> Self {
        WrappedFq3(self.0.double())
    }

    fn double_in_place(&mut self) -> &mut Self {
        self.0.double_in_place();
        self
    }

    fn neg_in_place(&mut self) -> &mut Self {
        self.0.neg_in_place();
        self
    }

    fn from_random_bytes_with_flags<F: ark_serialize::Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        Fq3::from_random_bytes_with_flags(bytes)
            .map(|(element, flags)| (WrappedFq3(element), flags))
    }

    fn legendre(&self) -> ark_ff::LegendreSymbol {
        self.0.legendre()
    }

    fn square(&self) -> Self {
        WrappedFq3(self.0.square())
    }

    fn square_in_place(&mut self) -> &mut Self {
        self.0.square_in_place();
        self
    }

    fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(WrappedFq3)
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        match self.0.inverse_in_place() {
            Some(_) => Some(self),
            None => None,
        }
    }

    fn frobenius_map(&mut self, power: usize) {
        self.0.frobenius_map(power)
    }
}

impl Ord for WrappedFq3 {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for WrappedFq3 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl CanonicalDeserializeWithFlags for WrappedFq3 {
    fn deserialize_with_flags<R: ark_serialize::Read, F: ark_serialize::Flags>(
        reader: R,
    ) -> Result<(Self, F), ark_serialize::SerializationError> {
        Fq3::deserialize_with_flags(reader).map(|(element, flags)| (WrappedFq3(element), flags))
    }
}

impl CanonicalSerializeWithFlags for WrappedFq3 {
    fn serialize_with_flags<W: ark_serialize::Write, F: ark_serialize::Flags>(
        &self,
        writer: W,
        flags: F,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.0.serialize_with_flags(writer, flags)
    }

    fn serialized_size_with_flags<F: ark_serialize::Flags>(&self) -> usize {
        self.0.serialized_size_with_flags::<F>()
    }
}

impl CanonicalSerialize for WrappedFq3 {
    fn serialize_with_mode<W: ark_serialize::Write>(
        &self,
        writer: W,
        compress: ark_serialize::Compress,
    ) -> Result<(), ark_serialize::SerializationError> {
        self.0.serialize_with_mode(writer, compress)
    }

    fn serialized_size(&self, compress: ark_serialize::Compress) -> usize {
        self.0.serialized_size(compress)
    }
}

impl CanonicalDeserialize for WrappedFq3 {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        reader: R,
        compress: ark_serialize::Compress,
        validate: ark_serialize::Validate,
    ) -> Result<Self, ark_serialize::SerializationError> {
        Fq3::deserialize_with_mode(reader, compress, validate).map(WrappedFq3)
    }
}

impl Valid for WrappedFq3 {
    fn check(&self) -> Result<(), ark_serialize::SerializationError> {
        self.0.check()
    }
}

impl Zero for WrappedFq3 {
    fn zero() -> Self {
        WrappedFq3(Fq3::zero())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl One for WrappedFq3 {
    fn one() -> Self {
        WrappedFq3(Fq3::one())
    }
}

impl Neg for WrappedFq3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        WrappedFq3(self.0.neg())
    }
}

// =====
// =====
// =====
// =====
// =====

impl core::iter::Product<Self> for WrappedFq3 {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), core::ops::Mul::mul)
    }
}

impl<'a> core::iter::Product<&'a Self> for WrappedFq3 {
    fn product<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), core::ops::Mul::mul)
    }
}

impl core::iter::Sum<Self> for WrappedFq3 {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), core::ops::Add::add)
    }
}

impl<'a> core::iter::Sum<&'a Self> for WrappedFq3 {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), core::ops::Add::add)
    }
}

impl Mul<WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn mul(self, rhs: WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 * rhs.0)
    }
}

impl Mul<&WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn mul(self, rhs: &WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 * rhs.0)
    }
}

impl Mul<&mut WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn mul(self, rhs: &mut WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 * rhs.0)
    }
}

impl Div<WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn div(self, rhs: WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 / rhs.0)
    }
}

impl Div<&WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn div(self, rhs: &WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 / rhs.0)
    }
}

impl Div<&mut WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn div(self, rhs: &mut WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 / rhs.0)
    }
}

impl Add<WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn add(self, rhs: WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 + rhs.0)
    }
}

impl Add<&WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn add(self, rhs: &WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 + rhs.0)
    }
}

impl Add<&mut WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn add(self, rhs: &mut WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 + rhs.0)
    }
}

impl Sub<WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn sub(self, rhs: WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 - rhs.0)
    }
}

impl Sub<&WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn sub(self, rhs: &WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 - rhs.0)
    }
}

impl Sub<&mut WrappedFq3> for WrappedFq3 {
    type Output = Self;

    fn sub(self, rhs: &mut WrappedFq3) -> Self::Output {
        WrappedFq3(self.0 - rhs.0)
    }
}

impl SubAssign<WrappedFq3> for WrappedFq3 {
    fn sub_assign(&mut self, rhs: WrappedFq3) {
        self.0 -= rhs.0;
    }
}

impl SubAssign<&WrappedFq3> for WrappedFq3 {
    fn sub_assign(&mut self, rhs: &WrappedFq3) {
        self.0 -= &rhs.0;
    }
}

impl SubAssign<&mut WrappedFq3> for WrappedFq3 {
    fn sub_assign(&mut self, rhs: &mut WrappedFq3) {
        self.0 -= &mut rhs.0;
    }
}

impl AddAssign<WrappedFq3> for WrappedFq3 {
    fn add_assign(&mut self, rhs: WrappedFq3) {
        self.0 += rhs.0;
    }
}

impl AddAssign<&WrappedFq3> for WrappedFq3 {
    fn add_assign(&mut self, rhs: &WrappedFq3) {
        self.0 += &rhs.0;
    }
}

impl AddAssign<&mut WrappedFq3> for WrappedFq3 {
    fn add_assign(&mut self, rhs: &mut WrappedFq3) {
        self.0 += &mut rhs.0;
    }
}

impl DivAssign<WrappedFq3> for WrappedFq3 {
    fn div_assign(&mut self, rhs: WrappedFq3) {
        self.0 /= rhs.0;
    }
}

impl DivAssign<&WrappedFq3> for WrappedFq3 {
    fn div_assign(&mut self, rhs: &WrappedFq3) {
        self.0 /= &rhs.0;
    }
}

impl DivAssign<&mut WrappedFq3> for WrappedFq3 {
    fn div_assign(&mut self, rhs: &mut WrappedFq3) {
        self.0 /= &mut rhs.0;
    }
}

impl MulAssign<WrappedFq3> for WrappedFq3 {
    fn mul_assign(&mut self, rhs: WrappedFq3) {
        self.0 *= rhs.0;
    }
}

impl MulAssign<&WrappedFq3> for WrappedFq3 {
    fn mul_assign(&mut self, rhs: &WrappedFq3) {
        self.0 *= &rhs.0;
    }
}

impl MulAssign<&mut WrappedFq3> for WrappedFq3 {
    fn mul_assign(&mut self, rhs: &mut WrappedFq3) {
        self.0 *= &mut rhs.0;
    }
}

impl From<u128> for WrappedFq3 {
    fn from(value: u128) -> Self {
        WrappedFq3(Fq3::from(value))
    }
}

impl From<u64> for WrappedFq3 {
    fn from(value: u64) -> Self {
        WrappedFq3(Fq3::from(value))
    }
}

impl From<u32> for WrappedFq3 {
    fn from(value: u32) -> Self {
        WrappedFq3(Fq3::from(value))
    }
}

impl From<u16> for WrappedFq3 {
    fn from(value: u16) -> Self {
        WrappedFq3(Fq3::from(value))
    }
}

impl From<u8> for WrappedFq3 {
    fn from(value: u8) -> Self {
        WrappedFq3(Fq3::from(value))
    }
}

impl From<bool> for WrappedFq3 {
    fn from(value: bool) -> Self {
        WrappedFq3(Fq3::from(value))
    }
}

// =====
// =====
// =====
// =====
// =====

impl MulAssign<Fp> for WrappedFq3 {
    fn mul_assign(&mut self, rhs: Fp) {
        self.0 *= Fq3::from_base_prime_field(rhs);
    }
}

impl From<Fp> for WrappedFq3 {
    fn from(value: Fp) -> Self {
        WrappedFq3(Fq3::from_base_prime_field(value))
    }
}

impl GpuField for WrappedFq3 {
    type FftField = Fp;

    fn field_name() -> String {
        todo!()
    }
}
