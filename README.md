> **Warning**
> this code is new and will change in future versions. Please try it out and provide feedback. If it addresses a use-case that is important to you please open an issue to discuss it or get in touch [andrew.j.milson@gmail.com](mailto:andrew.j.milson@gmail.com)

<div align="center">

# miniSTARK

**GPU accelerated STARK prover and verifier**

[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/andrewmilson/ministark/blob/main/LICENSE)
[![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#experimental)
[![CI](https://github.com/andrewmilson/ministark/actions/workflows/ci.yml/badge.svg)](https://github.com/andrewmilson/ministark/actions/workflows/ci.yml)


</div>

miniSTARK allows you to prove the integrity of arbitrary computations to anyone using the power of [STARKs](https://starkware.co/stark/). The library is written in [Rust](https://www.rust-lang.org/) but [Metal](https://developer.apple.com/metal/) is used to accelerate some polynomial arithmetic on the GPU. The Rust code uses several components from the fantastic [arkworks](https://github.com/arkworks-rs) library. The design of miniSTARK was influenced by [Winterfell](https://github.com/novifinancial/winterfell).

## Demo - proving brainf**k programs

| ![Generating a proof](prover.gif) | ![Verifying a proof](verifier.gif) |
|:--:|:--:|
| *Generating the proof* | *Verifying the proof* 

In this example the prover generates a proof that proves integrity of a brainf**k program that outputs "Hello World". The verifier uses the proof, brainf\*\*k source code and output to verify execution integrity without executing the program at all. To run this demo locally:

```bash
# (optional) compile shaders (M1 Mac only)
# make sure the latest Xcode is installed
# `(cd gpu-poly && make)`

# generate the proof
# use `-F parallel,asm` if not using an M1 Mac
# make sure latest macOS is installed
cargo +nightly run -r -F parallel,asm,gpu --example brainfuck -- \
    prove ./examples/brainfuck/hello_world.bf \
          --dst ./hello_world.proof

# verify the proof
cargo +nightly run -r -F asm --example brainfuck -- \
  verify ./examples/brainfuck/hello_world.bf \
         --output "Hello World" \
         --proof ./hello_world.proof 
```

This is actually a miniSTARK implementation of the [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) tutorial. This is an unrealistic example since verifying by running the program is actually much quicker than verifying by checking the proof. Generating a proof of "Hello World" or proving you can count from 1 to 10 is all fun and games but miniSTARK has much more serious ambitions. A realistic example is [coming soon](#coming-soon).

## Performance

Initial performance carried out on an M1 Max is promising. Compared to a couple of other Rust STARK provers miniSTARK generates proofs around **~2-50x** faster and consumes around **~2-40x** less RAM during proof generation. Since these comparisons were made with unrealistic toy examples they aren't entirely fair and won't be published. Performance results will be published once more realistic examples exist. Also, there are still a few easy performance optimizations to be made 😉.

## Defining AIR constraints

[AIR constraints](https://medium.com/starkware/arithmetization-i-15c046390862) are what the prover and verifier agree on to determine a valid execution trace. These constraints in miniSTARK are represented as multivariate polynomials where each variable abstractly represents either a column of the execution trace or one of the verifier's challenges. There are a lot of cool things the prover and verifier can do when constraints are represented in this way. Below is a contrived example to illustrate how constraints might be represented in Rust:

```rust
#[derive(Column)]
enum ProcessorTable {
    Cycle,
    InstructionPointer,
    CurrentInstruction,
    // ...
}

// all must evaluate to `0`
let constraints = vec![
    // cycle starts at `0`
    Cycle.first(),
    // each row, the cycle increases by `1`
    Cycle.curr() - Cycle.next() - 1,
    // ...
];
```

miniSTARK can use this symbolic form to evaluate a constraint and obtain its polynomial description. This is different from other STARK libraries that separate a constraint's evaluation and corresponding polynomial description. In these libraries, the evaluation is implemented as raw code which compilers are great at making fast but there are a couple of problems:

1. The evaluation and description parts of the codebase need to be maintained in parallel; and,
2. Making sure constraints are properly translated into code is hard - [Bobbin Threadbare talking about this at the ZK Summit](https://www.youtube.com/watch?v=81UAaiIgIYA&t=1383s).

The representation of constraints in miniSTARK is much closer to a representation you might find in a mathematical model and therefore far less error prone. The performance lost in not allowing the compiler to optimize evaluations is offset by evaluating the constraints in parallel on the GPU.


<h2 id="coming-soon">Coming soon (help wanted)</h2>

- Supporting proofs over secp256k1 field: <https://github.com/andrewmilson/ministark/issues/5>
- Polynomial arithmetic implemented in [CUDA](https://en.wikipedia.org/wiki/CUDA): <https://github.com/andrewmilson/ministark/issues/2>
- Speed and memory optimizations: <https://github.com/andrewmilson/ministark/issues/8>
- Using more `arkworks` features
- Reduce proof size using batched Merkle proofs: <https://github.com/andrewmilson/ministark/issues/10>
- Generate proofs with STARK friendly hash functions: <https://github.com/andrewmilson/ministark/issues/11>, <https://github.com/andrewmilson/ministark/issues/9>
- More tests and benchmarks: <https://github.com/andrewmilson/ministark/issues/3>
- More GPU field implementations: <https://github.com/andrewmilson/ministark/issues/1>
- Making gpu-poly less unsafe: <https://github.com/andrewmilson/ministark/issues/12>
- Generating zero knowledge proofs: <https://github.com/andrewmilson/ministark/issues/6>
- Realistic examples

## Acknowledgements

- [StarkWare](https://starkware.co/) - The company started by the people who invented STARKs. There are so many great things to learn from the people at this company. It's great how public and open their educational material is. Check out: GOAT Eli Ben-Sasson's [STARK math thread](https://twitter.com/EliBenSasson/status/1578380154476208131), Cairo whitepaper ([Video](https://www.youtube.com/watch?v=DTVn0oYLVsE), [PDF](https://eprint.iacr.org/2021/1063.pdf)), [DEEP Method Medium article](https://medium.com/starkware/starkdex-deep-dive-the-stark-core-engine-497942d0f0ab).
- [Alan Szepieniec](https://twitter.com/aszepieniec?lang=en) - The [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/) tutorial is (IMHO) the best end-to-end practical resource to understand STARKs. Read through this and the magic of STARKs will start to make sense. [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) is a fantastic sequel to the STARK Anatomy tutorial and is a practical guide to creating an AIR that proves programs written in the Brainf**k programming language. Check out miniSTARK's [brainf\*\*k example](examples/brainfuck/) which is an implementation of the BrainSTARK AIR.
- [Winterfell](https://github.com/novifinancial/winterfell) - A STARK prover and verifier developed at Facebook by [Bobbin Threadbare](https://twitter.com/bobbinth) and others. This repo was heavily used as a reference for several components: fast Merkle Tree, DEEP composition polynomial, public coin, channels, and STARK component traits. Bobbin has some great HackMD articles: [Miden VM program decoder](https://hackmd.io/_aaDBzbWRz6EwQQRtK1pzw), [Memory in Miden VM](https://hackmd.io/@bobbinth/HJr56BKKt), [u32 operations in Miden VM](https://hackmd.io/NC-yRmmtRQSvToTHb96e8Q#u32-operations-in-Miden-VM).
- [OpenZKP](https://github.com/0xProject/OpenZKP) - A STARK prover and verifier developed by [Remco Bloemen](https://twitter.com/recmo) and others. OpenZKP inspired the AIR constraints and some traits used in miniSTARK. Remco also created a [nice short STARK explanation video](https://www.youtube.com/watch?v=H3AKu03AwYc).