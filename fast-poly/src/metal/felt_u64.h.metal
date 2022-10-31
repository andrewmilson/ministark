#ifndef felt_u64_h
#define felt_u64_h

#include "u128.h.metal"

class FP18446744069414584321
{
public:
    FP18446744069414584321() = default;
    constexpr FP18446744069414584321(unsigned long v) : inner(v) {}

    constexpr FP18446744069414584321 operator+(const FP18446744069414584321 rhs) const
    {
        return FP18446744069414584321(add(inner, rhs.inner));
    }

    constexpr FP18446744069414584321 operator-(const FP18446744069414584321 rhs) const
    {
        return FP18446744069414584321(sub(inner, rhs.inner));
    }

    FP18446744069414584321 operator*(const FP18446744069414584321 rhs) const
    {
        return FP18446744069414584321(mul(inner, rhs.inner));
    }

    FP18446744069414584321 pow(unsigned exp) 
    {
        FP18446744069414584321 res = ONE;

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
    unsigned long inner;

    // Field modulus `p = 2^64 - 2^32 + 1`
    constexpr static const constant unsigned long N = 18446744069414584321;

    // Square of auxiliary modulus R for Montgomery reduction `R2 ≡ (2^64)^2 mod p`
    constexpr static const constant unsigned long R2 = 18446744065119617025;

    // // 1 in Montgomery representation
    constexpr static const constant unsigned long ONE = 4294967295;

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
        unsigned long xl = lhs * rhs;// x.low;
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

#endif /* felt_u64_h */
