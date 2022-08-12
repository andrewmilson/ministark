use super::AffinePoint;
use super::CurveConfig;
use super::WeierstrassCurveConfig;
use crate::fields::fp_u256::BaseFelt;
use crate::fields::fp_u256::U256;

/// secp256k curve
struct Curve;

impl CurveConfig for Curve {
    type BaseFelt = BaseFelt;
}

impl WeierstrassCurveConfig for Curve {
    const A: BaseFelt = BaseFelt::ZERO;
    const B: BaseFelt = BaseFelt::new(U256 { high: 0, low: 7 });
}

/// Generator point of the abelian group used in Bitcoin
const GENERATOR: AffinePoint<Curve> = AffinePoint::new_unchecked(
    BaseFelt::new(U256 {
        high: 0x79BE667EF9DCBBAC55A06295CE870B07,
        low: 0x029BFCDB2DCE28D959F2815B16F81798,
    }),
    BaseFelt::new(U256 {
        high: 0x483ADA7726A3C4655DA4FBFC0E1108A8,
        low: 0xFD17B448A68554199C47D08FFB10D4B8,
    }),
);

const ORDER: U256 = U256 {
    high: 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE,
    low: 0xBAAEDCE6AF48A03BBFD25E8CD0364141,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ec::AffinePoint;

    #[test]
    fn generator_has_correct_order() {
        assert_eq!(GENERATOR * ORDER, AffinePoint::<Curve>::identity());
    }
}
