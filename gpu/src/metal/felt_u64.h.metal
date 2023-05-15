#ifndef felt_u64_h
#define felt_u64_h

// Fields that use prime 18446744069414584321
namespace p18446744069414584321
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

        constexpr explicit operator unsigned long() const
        {
            return inner;
        }

        Fp operator*(const Fp rhs) const
        {
            return Fp(mul(inner, rhs.inner));
        }

        // used for S-box in the Rescue Prime Optimized hash function
        Fp pow7() {
            unsigned long t2 = mul(inner, inner);
            return mul(mul(mul(t2, t2), t2), inner);
        }

        // used for S-box-inverse in the Rescue Prime Optimized hash function
        Fp pow10540996611094048183() {
            // ADDCHAIN IMPLEMENTATION
            // // used addchain
            // // https://github.com/mmcloughlin/addchain
            // unsigned long _100 = sqn<2>(inner);
            // unsigned long _100000 = sqn<3>(_100);
            // unsigned long _100001 = mul(_100000, inner);
            // unsigned long i9 = sqn<3>(_100000);
            // unsigned long i10 = mul(_100001, i9);
            // unsigned long i14 = mul(sqn<3>(i9), i10);
            // unsigned long i17 = mul(sqn<1>(mul(_100, i14)), inner);
            // unsigned long i18 = mul(i14, i17);
            // unsigned long i19 = mul(i17, i18);
            // unsigned long i20 = mul(i18, i19);
            // unsigned long i21 = mul(inner, i20);
            // unsigned long i23 = mul(sqn<1>(i21), inner);
            // unsigned long i60 = mul(mul(sqn<17>(mul(sqn<17>(i23), i20)), i23), i21);
            // return mul(sqn<14>(i60), i19);

            // FROM MIDEN IMPLEMENTATION
            // compute base^10540996611094048183 using 72 multiplications per array element
            // 10540996611094048183 = b1001001001001001001001001001000110110110110110110110110110110111
            unsigned long t1 = mul(inner, inner);
            unsigned long t2 = mul(t1, t1);
            unsigned long t3 = exp_acc<3>(t2, t2);
            unsigned long t4 = exp_acc<6>(t3, t3);
            unsigned long t5 = exp_acc<12>(t4, t4);
            unsigned long t6 = exp_acc<6>(t5, t3);
            unsigned long t7 = exp_acc<31>(t6, t6);
            unsigned long t8 = sqn<2>(mul(mul(t7, t7), t6));
            return mul(mul(mul(t1, t2), inner), t8);
        }

        Fp pow(unsigned exp)
        {
            // TODO: consider removing.
            // if (exp == 1) {
            //     return *this;
            // }

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

        Fp inverse() 
        {
            unsigned long t2 = exp_acc<1>(inner, inner);
            unsigned long t3 = exp_acc<1>(t2, inner);
            unsigned long t6 = exp_acc<3>(t3, t3);
            unsigned long t12 = exp_acc<6>(t6, t6);
            unsigned long t24 = exp_acc<12>(t12, t12);
            unsigned long t30 = exp_acc<6>(t24, t6);
            unsigned long t31 = exp_acc<1>(t30, inner);
            unsigned long t63 = exp_acc<32>(t31, t31);
            unsigned long inv = exp_acc<1>(t63, inner);
            return Fp(inv);
        }

        Fp neg()
        {
            // TODO: can improve
            return Fp(sub(0, inner));
        }
        
        // 1 in Montgomery representation
        constexpr static const constant unsigned long ONE = 4294967295;

    private:
        unsigned long inner;

        // Field modulus `p = 2^64 - 2^32 + 1`
        constexpr static const constant unsigned long N = 18446744069414584321;

        // Square of auxiliary modulus R for Montgomery reduction `R2 ≡ (2^64)^2 mod p`
        constexpr static const constant unsigned long R2 = 18446744065119617025;

        template<unsigned N_ACC>
        inline unsigned long exp_acc(unsigned long element, const unsigned long tail) const {
#pragma unroll
            for (unsigned i = 0; i < N_ACC; i++) {
                element = mul(element, element);
            }
            return mul(element, tail);
        }

        template<unsigned N_ACC>
        inline unsigned long sqn(unsigned long element) const {
#pragma unroll
            for (unsigned i = 0; i < N_ACC; i++) {
                element = mul(element, element);
            }
            return element;
        }

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
            unsigned long xl = lhs * rhs;
            unsigned long xh = metal::mulhi(lhs, rhs);
            unsigned long tmp = xl << 32;
            unsigned a_overflow = xl > (0xFFFFFFFFFFFFFFFF - tmp);
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
        constexpr Fq3(Fp c) : c0(c), c1(Fp(0)), c2(Fp(0)) {}
        constexpr Fq3(Fp c0, Fp c1, Fp c2) : c0(c0), c1(c1), c2(c2) {}

        constexpr Fq3 operator+(const Fp rhs) const
        {
            return Fq3(c0 + rhs, c1, c2);
        }

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
            // Devegili OhEig Scott Dahab --- Multiplication and Squaring on
            // AbstractPairing-Friendly
            // Fields.pdf; Section 4 (Karatsuba)
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
            // TODO: investige. for some reason removing this makes tests fail
            if (exp == 1) {
                return *this;
            }

            Fq3 res = Fq3(Fp(Fp::ONE), Fp(0), Fp(0));

            while (exp > 0) {
                if (exp & 1) {
                    res = res * *this;
                }
                exp >>= 1;
                // TODO: this doesn't need to be done on the last itteration.
                // may already be optimized away by the compiler. Not sure.
                *this = *this * *this;
            }

            return res;
        }

        Fq3 neg()
        {
            // TODO: can improve
            return Fq3(Fp(0) - c0, Fp(0) - c1, Fp(0) - c2);
        }

        // Fq3 inverse() 
        // {
        //     // TODO
        // }

    private:
        Fp c0, c1, c2;

        // Cubic non-residue used to construct the extension field in montgomery representation. 
        // That is, `NONRESIDUE` is such that the cubic polynomial
        // `f(X) = X^3 - NONRESIDUE` in Fp\[X\] is irreducible in `Fp`.
        constexpr static const constant unsigned long NONREDIDUE = /* =2 */ 8589934590;
    };

}

#endif /* felt_u64_h */
