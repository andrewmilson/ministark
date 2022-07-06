# PolySonic

Optimized for Apple Silicon the PolySonic finite-field polynomial library uses the fast fourier transform, GPU acceleration, and Montgomery arithmetic to achieve ridiculous 100-1000x speedups over CPU polynomial libraries.

- **100x faster polynomial multiplication**
- **100x faster polynomial interpolation**
- **100x faster polynomial evaluation**

# Usage

```bash
# debug mode
export METAL_DEVICE_WRAPPER_TYPE=1
make
cargo test
```
