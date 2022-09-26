# Optimized `arkworks` field implementations

Specialized field implementations that outperform the generic, Montgomery arithmetic, field implementations provided in [arkworks](https://github.com/arkworks-rs/algebra). Currently only one field has been added. Please contribute an implementation! PRs are welcomed!

# Prime field `p=18446744069414584321`

```rust
use ark_ff_optimized::fp64::Fp;
```

An amazing prime with modulus `p=2^64−2^32+1`. This field has some nice properties (1) Multiplying two 32-bit values does not overflow the field modulus and (2) Checking whether four 16-bit values form a valid field element can be done efficiently. This field is used in [Polygon Miden and Polygon Zero](https://twitter.com/0xPolygonMiden/status/1478786573348995075). Implementation was sourced from [EcGFp5: a Specialized Elliptic Curve](https://eprint.iacr.org/2022/274.pdf) and [Facebook's Winterfell repo](https://github.com/novifinancial/winterfell/blob/c7c92620cc7661e38ad58e1a3bdfbd7bba694c5d/math/src/field/f64/mod.rs).

| `Benchmark`                              | `Generic`   | `Specialized` (this repo)         |
| :--------------------------------------- | :---------- | :-------------------------------- |
| **`Sum of products of size 2`**          | `18.04 ns`  | `7.34 ns` (🚀 **2.46x faster**)   |
| **`Inverse`**                            | `556.74 ns` | `283.87 ns` (🚀 **1.96x faster**) |
| **`Legendre for QR`**                    | `1.12 us`   | `596.15 ns` (🚀 **1.88x faster**) |
| **`Naive sum of products of size 2`**    | `15.41 ns`  | `8.68 ns` (🚀 **1.78x faster**)   |
| **`Deserialize Compressed`**             | `8.82 ns`   | `4.99 ns` (🚀 **1.77x faster**)   |
| **`Deserialize Compressed Unchecked`**   | `8.80 ns`   | `4.97 ns` (🚀 **1.77x faster**)   |
| **`Deserialize Uncompressed`**           | `8.86 ns`   | `5.16 ns` (🚀 **1.72x faster**)   |
| **`Deserialize Uncompressed Unchecked`** | `8.81 ns`   | `5.15 ns` (🚀 **1.71x faster**)   |
| **`Square Root for QR`**                 | `4.43 us`   | `2.77 us` (🚀 **1.60x faster**)   |
| **`Multiplication`**                     | `6.15 ns`   | `4.03 ns` (🚀 **1.53x faster**)   |
| **`From BigInt`**                        | `5.32 ns`   | `4.30 ns` (✅ **1.24x faster**)   |
| **`Serialize Uncompressed`**             | `4.72 ns`   | `3.95 ns` (✅ **1.20x faster**)   |
| **`Into BigInt`**                        | `4.72 ns`   | `3.92 ns` (✅ **1.20x faster**)   |
| **`Serialize Compressed`**               | `4.72 ns`   | `3.96 ns` (✅ **1.19x faster**)   |
| **`Square`**                             | `5.60 ns`   | `4.88 ns` (✅ **1.15x faster**)   |
| **`Subtraction`**                        | `4.09 ns`   | `3.77 ns` (✅ **1.09x faster**)   |
| **`Addition`**                           | `4.11 ns`   | `3.79 ns` (✅ **1.08x faster**)   |
| **`Negation`**                           | `4.21 ns`   | `3.90 ns` (✅ **1.08x faster**)   |
| **`Double`**                             | `4.13 ns`   | `4.32 ns` (❌ **1.04x slower**)   |

_Benchmarked on an M1 Max. Markdown generated with
[criterion-table](https://github.com/nu11ptr/criterion-table). More detailed benchmark info is [here](http://andrewmilson.com/optimized-fields/criterion/report/index.html)_
