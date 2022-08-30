use algebra::fp_u256::U256;
use algebra::Felt;
use num_traits::One;
use num_traits::Zero;
use std::fmt::Debug;
use std::fmt::Display;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;

pub mod schoofs;
pub mod secp256k;
/// Config of an elliptic curve
pub trait CurveConfig: Sized {
    type BaseFelt: Felt;
}

/// Config for curves in Weierstrass Form i.e. curves of `y^2 = x^3 + ax + b`
pub trait WeierstrassCurveConfig: CurveConfig {
    /// Coefficient `a` of the curve equation.
    const A: Self::BaseFelt;

    /// Coefficient `b` of the curve equation.
    const B: Self::BaseFelt;
}

/// Affine point on a Weierstrass curve
pub struct AffinePoint<C: WeierstrassCurveConfig> {
    pub x: C::BaseFelt,
    pub y: C::BaseFelt,
    pub infinity: bool,
}

impl<C: WeierstrassCurveConfig> Debug for AffinePoint<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AffinePoint")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("infinity", &self.infinity)
            .finish()
    }
}

impl<C: WeierstrassCurveConfig> Clone for AffinePoint<C> {
    fn clone(&self) -> Self {
        Self {
            x: self.x,
            y: self.y,
            infinity: self.infinity,
        }
    }
}

impl<C: WeierstrassCurveConfig> Copy for AffinePoint<C> {}

impl<C: WeierstrassCurveConfig> PartialEq<AffinePoint<C>> for AffinePoint<C> {
    fn eq(&self, other: &AffinePoint<C>) -> bool {
        self.x == other.x && self.y == other.y && self.infinity == other.infinity
    }
}

impl<C: WeierstrassCurveConfig> Display for AffinePoint<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.infinity {
            true => write!(f, "infinity"),
            false => write!(f, "({}, {})", self.x, self.y),
        }
    }
}

impl<C: WeierstrassCurveConfig> AffinePoint<C> {
    pub const fn new_unchecked(x: C::BaseFelt, y: C::BaseFelt) -> Self {
        AffinePoint {
            x,
            y,
            infinity: false,
        }
    }

    pub fn identity() -> Self {
        // TODO: make constant
        Self {
            x: C::BaseFelt::zero(),
            y: C::BaseFelt::zero(),
            infinity: true,
        }
    }

    pub fn is_on_curve(&self) -> bool {
        if self.infinity {
            true
        } else {
            self.y.square() == self.x * self.x * self.x + C::A * self.x + C::B
        }
    }
}

impl<C: WeierstrassCurveConfig> Zero for AffinePoint<C> {
    /// Returns the point at infinity
    fn zero() -> Self {
        Self::identity()
    }

    /// Checks if `self` is the point at infinity.
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
}

impl<C: WeierstrassCurveConfig> Add for AffinePoint<C> {
    type Output = AffinePoint<C>;

    fn add(self, rhs: AffinePoint<C>) -> Self::Output {
        let mut copy = self;
        copy += rhs;
        copy
    }
}

impl<C: WeierstrassCurveConfig> AddAssign for AffinePoint<C> {
    fn add_assign(&mut self, rhs: AffinePoint<C>) {
        *self = if self.is_zero() {
            rhs
        } else if rhs.is_zero() {
            *self
        } else if self.x == rhs.x && self.y == -rhs.y {
            Self::identity()
        } else if self.x != rhs.x {
            let Self { x: x1, y: y1, .. } = *self;
            let Self { x: x2, y: y2, .. } = rhs;
            let slope = (y2 - y1) / (x2 - x1);
            let x3 = slope.square() - x1 - x2;
            let y3 = slope * (x1 - x3) - y1;
            Self::new_unchecked(x3, y3)
        } else {
            let Self { x: x1, y: y1, .. } = *self;
            let x1_squared = x1.square();
            let slope = (x1_squared + x1_squared + x1_squared + C::A) / y1.double(); // ∂y/∂x
            let x3 = slope.square() - x1.double();
            let y3 = slope * (x1 - x3) - y1;
            Self::new_unchecked(x3, y3)
        }
    }
}

impl<C: WeierstrassCurveConfig> Mul<U256> for AffinePoint<C> {
    type Output = AffinePoint<C>;

    fn mul(self, rhs: U256) -> Self::Output {
        let mut rhs = rhs;
        let mut res = Self::zero();
        let mut accumulator = self;
        while !rhs.is_zero() {
            if (rhs & U256::one()).is_one() {
                res += accumulator;
            }

            accumulator += accumulator;
            rhs >>= U256::one();
        }
        res
    }
}

/// Projective point on a Weierstrass curve
#[derive(Clone, Copy, Debug)]
pub struct ProjectivePoint<C: WeierstrassCurveConfig> {
    pub x: C::BaseFelt,
    pub y: C::BaseFelt,
    pub z: C::BaseFelt,
}

impl<C: WeierstrassCurveConfig> ProjectivePoint<C> {
    pub fn new_unchecked(x: C::BaseFelt, y: C::BaseFelt, z: C::BaseFelt) -> Self {
        ProjectivePoint { x, y, z }
    }
}

impl<C: WeierstrassCurveConfig> Add for ProjectivePoint<C> {
    type Output = Self;

    fn add(self, _rhs: Self) -> Self::Output {
        todo!()
    }
}

impl<C: WeierstrassCurveConfig> Zero for ProjectivePoint<C> {
    /// Returns the point at infinity which always has z=0
    fn zero() -> Self {
        Self::new_unchecked(C::BaseFelt::one(), C::BaseFelt::one(), C::BaseFelt::zero())
    }

    fn is_zero(&self) -> bool {
        self.z.is_zero()
    }
}
