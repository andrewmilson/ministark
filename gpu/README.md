# GPU optimized STARK/SNARK primitives

This library contains GPU optimized primitives commonly used by STARK/SNARK provers. The code is written in the Metal programming language so currently only supports Apple Silicon. Part of the development of this library was supported by the fantastic [Polygon Miden](https://github.com/0xPolygonMiden/miden-vm) team.

## Usage

```bash
# recompile shaders
make shaders

# run tests
export METAL_DEVICE_WRAPPER_TYPE=1
cargo test

# run benchmarks
cargo bench
```
