#ifndef felt_u64_h
#define felt_u64_h

#include "u128.h.metal"

// Fields that use prime 18446744069414584321
namespace P18446744069414584321
{

    // Prime field
    class Fp
    {
    public:
        Fp() = default;
        constexpr Fp(unsigned long v) : inner(v) {}

        constexpr Fp operator+(const Fp rhs) const
        {
            return Fp(add(inner, rhs.inner));
        }

        constexpr Fp operator-(const Fp rhs) const
        {
            return Fp(sub(inner, rhs.inner));
        }

        Fp operator*(const Fp rhs) const
        {
            return Fp(mul(inner, rhs.inner));
        }

        Fp pow(unsigned exp)
        {
            Fp res = ONE;

            while (exp > 0)
            {
                if (exp & 1)
                {
                    res = res * *this;
                }
                exp >>= 1;
                *this = *this * *this;
            }

            return res;
        }
        
        // 1 in Montgomery representation
        constexpr static const constant unsigned long ONE = 4294967295;

    private:
        unsigned long inner;

        // Field modulus `p = 2^64 - 2^32 + 1`
        constexpr static const constant unsigned long N = 18446744069414584321;

        // Square of auxiliary modulus R for Montgomery reduction `R2 ≡ (2^64)^2 mod p`
        constexpr static const constant unsigned long R2 = 18446744065119617025;


        inline unsigned long add(const unsigned long a, const unsigned long b) const
        {
            // We compute a + b = a - (p - b).
            unsigned long tmp = N - b;
            unsigned underflow = a < tmp;
            unsigned long x1 = a - tmp;
            unsigned adj = -underflow;
            return x1 - adj;
        }

        inline unsigned long sub(const unsigned long a, const unsigned long b) const
        {
            unsigned underflow = a < b;
            unsigned long x1 = a - b;
            unsigned adj = -underflow;
            return x1 - adj;
        }

        inline unsigned long mul(const unsigned long lhs, const unsigned long rhs) const
        {
            // u128 x = u128(lhs) * u128(rhs);
            unsigned long xl = lhs * rhs; // x.low;
            unsigned long xh = metal::mulhi(lhs, rhs);
            unsigned long tmp = xl << 32;
            unsigned long a_overflow = xl > (0xFFFFFFFFFFFFFFFF - tmp);
            unsigned long a = xl + tmp;
            unsigned long b = a - (a >> 32) - a_overflow;
            unsigned r_underflow = xh < b;
            unsigned long r = xh - b;
            unsigned adj = -r_underflow;
            return r - adj;
        }
    };

    // Cubig extension field over irreducible polynomial x^3 - 2
    // adapted from arkworks
    // TODO: make into generic cubic extension class
    class Fq3
    {
    public:
        Fq3() = default;
        constexpr Fq3(Fp c0, Fp c1, Fp c2) : c0(c0), c1(c1), c2(c2) {}

        constexpr Fq3 operator+(const Fq3 rhs) const
        {
            return Fq3(c0 + rhs.c0, c1 + rhs.c1, c2 + rhs.c2);
        }

        constexpr Fq3 operator-(const Fq3 rhs) const
        {
            return Fq3(c0 - rhs.c0, c1 - rhs.c1, c2 - rhs.c2);
        }

        Fq3 operator*(const Fq3 rhs) const 
        {
            Fp a = rhs.c0;
            Fp b = rhs.c1;
            Fp c = rhs.c2;

            Fp d = c0;
            Fp e = c1;
            Fp f = c2;

            Fp ad = d * a;
            Fp be = e * b;
            Fp cf = f * c;

            Fp x = (e + f) * (b + c) - be - cf;
            Fp y = (d + e) * (a + b) - ad - be;
            Fp z = (d + f) * (a + c) - ad + be - cf;

            return Fq3(
                /* =c0 */ ad + x * Fp(NONREDIDUE), 
                /* =c1 */ y + cf * Fp(NONREDIDUE),
                /* =c2 */ z
            );
        }

        Fq3 operator*(const Fp rhs) const
        {
            return Fq3(c0 * rhs, c1 * rhs, c2 * rhs);
        }

        // TODO: make util function
        Fq3 pow(unsigned exp)
        {
            Fq3 res = Fq3(Fp(Fp::ONE), Fp(0), Fp(0));

            while (exp > 0) {
                if (exp & 1) {
                    res = res * *this;
                }
                exp >>= 1;
                *this = *this * *this;
            }

            return res;
        }

    private:
        Fp c0, c1, c2;

        // Cubic non-residue used to construct the extension field in montgomery representation. 
        // That is, `NONRESIDUE` is such that the cubic polynomial
        // `f(X) = X^3 - NONRESIDUE` in Fp\[X\] is irreducible in `Fp`.
        constexpr static const constant unsigned long NONREDIDUE = /* =2 */ 8589934590;

    };

}

#endif /* felt_u64_h */
