#ifndef u256_h
#define u256_h

#include "u128.h.metal"

class u256
{
public:
    u256() = default;
    constexpr u256(int l) : low(l), high(0) {}
    constexpr u256(unsigned long l) : low(u128(l)), high(0) {}
    constexpr u256(u128 l) : low(l), high(0) {}
    constexpr u256(bool b) : low(b), high(0) {}
    constexpr u256(u128 h, u128 l) : low(l), high(h) {}
    constexpr u256(unsigned long hh, unsigned long hl, unsigned long lh, unsigned long ll) : 
        low(u128(lh, ll)), high(u128(hh, hl)) {}

    constexpr u256 operator+(const u256 rhs) const
    {
        return u256(high + rhs.high + ((low + rhs.low) < low), low + rhs.low);
    }

    constexpr u256 operator+=(const u256 rhs)
    {
        *this = *this + rhs;
        return *this;
    }

    constexpr inline u256 operator-(const u256 rhs) const
    {
        return u256(high - rhs.high - ((low - rhs.low) > low), low - rhs.low);
    }

    constexpr u256 operator-=(const u256 rhs)
    {
        *this = *this - rhs;
        return *this;
    }

    constexpr bool operator==(const u256 rhs) const
    {
        return high == rhs.high && low == rhs.low;
    }

    constexpr bool operator!=(const u256 rhs) const
    {
        return !(*this == rhs);
    }

    constexpr bool operator<(const u256 rhs) const
    {
        return ((high == rhs.high) && (low < rhs.low)) || (high < rhs.high);
    }

    constexpr u256 operator&(const u256 rhs) const
    {
        return u256(high & rhs.high, low & rhs.low);
    }

    constexpr bool operator>(const u256 rhs) const
    {
        return ((high == rhs.high) && (low > rhs.low)) || (high > rhs.high);
    }

    constexpr bool operator>=(const u256 rhs) const
    {
        return !(*this < rhs);
    }

    constexpr bool operator<=(const u256 rhs) const
    {
        return !(*this > rhs);
    }

    constexpr inline u256 operator>>(unsigned shift) const
    {
        // TODO: reduce branch conditions
        if (shift >= 256)
        {
            return u256(0);
        }
        else if (shift == 128)
        {
            return u256(0, high);
        }
        else if (shift == 0)
        {
            return *this;
        }
        else if (shift < 128)
        {
            return u256(high >> shift, (high << (128 - shift)) | (low >> shift));
        }
        else if ((256 > shift) && (shift > 128))
        {
            return u256(0, (high >> (shift - 128)));
        }
        else
        {
            return u256(0);
        }
    }

    constexpr u256 operator>>=(unsigned rhs)
    {
        *this = *this >> rhs;
        return *this;
    }

    u256 operator*(const bool rhs) const
    {
        return u256(high * rhs, low * rhs);
    }

    u256 operator*(const u256 rhs) const
    {
        // split values into 4 64-bit parts
        u128 top[4] = {u128(high.high), u128(high.low), u128(low.high), u128(low.low)};
        u128 bottom[4] = {u128(rhs.high.high), u128(rhs.high.low), u128(rhs.low.high), u128(rhs.low.low)};
        // u128 top[4] = {high >> 32, high & 0xffffffff, low >> 32, low & 0xffffffff};
        // u128 bottom[4] = {rhs.high >> 32, rhs.high & 0xffffffff, rhs.low >> 32, rhs.low & 0xffffffff};
        u128 products[4][4];

        // // multiply each component of the values
        // Alternative:
        //   for(int y = 3; y > -1; y--){
        //       for(int x = 3; x > -1; x--){
        //           products[3 - x][y] = top[x] * bottom[y];
        //       }
        //   }
        products[0][3] = top[3] * bottom[3];
        products[1][3] = top[2] * bottom[3];
        products[2][3] = top[1] * bottom[3];
        products[3][3] = top[0] * bottom[3];

        products[0][2] = top[3] * bottom[2];
        products[1][2] = top[2] * bottom[2];
        products[2][2] = top[1] * bottom[2];
        // products[3][2] = top[0] * bottom[2];

        products[0][1] = top[3] * bottom[1];
        products[1][1] = top[2] * bottom[1];
        // products[2][1] = top[1] * bottom[1];
        products[3][1] = top[0] * bottom[1];

        products[0][0] = top[3] * bottom[0];
        // products[1][0] = top[2] * bottom[0];
        // products[2][0] = top[1] * bottom[0];
        // products[3][0] = top[0] * bottom[0];

        // first row
        u128 fourth64 = products[0][3].low;
        u128 third64 = u128(products[0][2].low) + u128(products[0][3].high);
        u128 second64 = u128(products[0][1].low) + u128(products[0][2].high);
        u128 first64 = u128(products[0][0].low) + u128(products[0][1].high);

        // second row
        third64 += u128(products[1][3].low);
        second64 += u128(products[1][2].low) + u128(products[1][3].high);
        first64 += u128(products[1][1].low) + u128(products[1][2].high);

        // third row
        second64 += u128(products[2][3].low);
        first64 += u128(products[2][2].low) + u128(products[2][3].high);

        // fourth row
        first64 += u128(products[3][3].low);

        return u256(u128(first64.low, 0), 0) +
               u256(third64.high, u128(third64.low, 0)) +
               u256(second64, 0) +
               u256(fourth64);

        // // move carry to next digit
        // // third64 += fourth64 >> 64; // TODO: figure out if this is a nop
        // second64 += u128(third64.high);
        // first64 += u128(second64.high);

        // // remove carry from current digit
        // // fourth64 &= 0xffffffff; // TODO: figure out if this is a nop
        // // third64 &= 0xffffffff;
        // // second64 = u128(second64.low);
        // // first64 &= 0xffffffff;

        // // combine components
        // // return u256((first64 << 64) | second64, (third64 << 64) | fourth64);
        // return u256(u128(first64.low, second64.low), u128(third64.low, fourth64.low));

        // return u128((first64.high second64, (third64 << 64) | fourth64);
    }

    u256 operator*=(const u256 rhs)
    {
        *this = *this * rhs;
        return *this;
    }

    // TODO: Could get better performance with  smaller limb size
    // Not sure what word size is for M1 GPU
#ifdef __LITTLE_ENDIAN__
    u128 low;
    u128 high;
#endif
#ifdef __BIG_ENDIAN__
    u128 high;
    u128 low;
#endif
};

#endif /* u256_h */
