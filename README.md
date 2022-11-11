> **Warning**
> this code is new and will change in future versions. Please try it out and provide feedback. If it addresses a use-case that is important to you please open an issue to discuss it or register your interest https://forms.gle/SYrA9hbGQuivbHX5A

<div align="center">

# miniSTARK

**GPU accelerated STARK prover and verifier**

[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/andrewmilson/mini-stark/blob/main/LICENSE)
[![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#experimental)

</div>

miniSTARK allows you to prove the integrity of arbitrary computations to anyone using the power of [STARKs](https://starkware.co/stark/). In the silly example above miniSTARK is used to prove the integrity of a program that outputs "Hello World!" in the [Brainf**k](https://esolangs.org/wiki/Brainfuck) programming language. Anyone can verify the proof instantly.

The library is written [Rust](https://www.rust-lang.org/) but [Metal](https://developer.apple.com/metal/) is used to GPU accelerate polynomial arithmetic. The Rust code uses components from the [arkworks](https://github.com/arkworks-rs) library. The design of miniSTARK was influenced by [Winterfell](https://github.com/novifinancial/winterfell).

## Demo

### Proving "Hello World!" in brainf**k

In this example the prover generates a proof that proves integrity of a brainf**k program that outputs "Hello World!". The verifier uses the proof and program source code to verify the output without having to run the program. It's a fun but silly example since it's quicker to execute the program than verify the proof. To run the example:

```bash
# generate the proof
# use `-F parallel,asm` if not using an M1 MacBook
cargo run -r -F parallel,asm,gpu --example brainfuck -- \
    prove --src ./examples/brainfuck/hello_world.bf \
          --output ./hello_world.proof


# verify the proof
cargo run -r -F asm --example brainfuck -- \
  verify --src ./examples/brainfuck/hello_world.bf \
         --proof ./hello_world.proof 
```

This is actually a miniSTARK implementation of the [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) tutorial. The tutorial provides a [Rust implementation](https://github.com/Neptune-Crypto/twenty-first/tree/master/stark-brainfuck) built using the [twenty-first](https://github.com/Neptune-Crypto/twenty-first) library. Below is a comparison between the twenty-first and miniSTARK prover. It isn't entirely fair for a couple of reasons (1) miniSTARK doesn't generate zero knowledge proofs and (2) each prover commits to constraints in a different way. 



### Proving Fibonacci sequence

## Performance


## Defining AIR Constraints

Constraints in miniSTARK are represented as multivariate polynomials where each variable abstractly represents a column of the execution trace or one of the verifier's challenges. There is a lot of cool things prover and verifier can do when constraints are represented in this way. Bellow is a contrived example to illustrate how constraints are represented. Look at the constraint implementation of [brainf**k](examples/brainfuck/) for a full example.

```rust
#[derive(Hint)]
enum Hints {
    ExpectedCycles,
}

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
    // cycle doesn't exceed expected
    Cycle.last() - ExpectedCycles,
    // ...
];
```

<!-- ## Examples

### Brainf*** virtual machine

Implementation of the [Brainf***](https://esolangs.org/wiki/Brainfuck) virtual machine from [Alan Szepieniec BrainSTARK tutorial](https://aszepieniec.github.io/stark-brainfuck/brainfuck).

```bash
# source: https://esolangs.org
export HELLO_WORLD_BF="++++++++[>++++[>++>+++>+++>+<<<<-]>+>+\
>->>+[<]<-]>>.>---.+++++++..+++.>>.<-.<.+++.------.--------.>>+.>++."

if [[ $(arch) == 'arm64' ]]; then
  # run on the GPU if Apple silicon
  cargo run --release --features parallel,asm,gpu --example bf --src $HELLO_WORLD_BF
else
  # fall back to cpu if not Apple silicon
  cargo run --release --features parallel,asm --example bf --src $HELLO_WORLD_BF
fi

```

### Multiplicative Fibonacci Sequence 

An analogue to the regular fibonacci sequence that uses multiplication rather than addition. Multiplicative fibonacci requires more grunt (more AIR constraints) to prove. Sequence is `1, 2, 2, 4, 8, ...`. The program proves the 

```bash
cargo run --release --features parallel,asm  --example fib
```

## Things I don't like

- remembering what the longest table is
- all terminal and challence and column you have to remember the numerical index. Could implement trait and enums to mitigate this.

## TODO

- debugging memory table. Remove constraint for memory stay the same clock cycle increase -->

<h2 id="coming-soon">Coming soon (help wanted)</h2>

- Supporting proofs over secp256k1 aka Bitcoin's field (see [ecfft-bn254](https://github.com/wborgeaud/ecfft-bn254))
- Polynomial arithmetic implemented in [CUDA](https://en.wikipedia.org/wiki/CUDA)
- Speed and memory optimizations
- Using more `arkworks` features
- Reduce proof size using batched Merkle proofs
- Generate proofs with STARK friendly hash functions
- Cairo VM prover (similar to [giza](https://github.com/maxgillett/giza))
- More tests and benchmarks
- More GPU field implementations
- Making gpu-poly less unsafe
- Generating zero knowledge proofs
- Periodic constraints

## Acknowledgements

- [StarkWare](https://starkware.co/) - The company started by the people who invented STARKs. There's so many great things to learn from the people at this company. It's great how public and open their educational material is. Check out: GOAT Eli Ben-Sasson's [STARK math thread](https://twitter.com/EliBenSasson/status/1578380154476208131), [STARK 101](https://starkware.co/stark-101/), Cairo whitepaper ([Video](https://www.youtube.com/watch?v=DTVn0oYLVsE), [PDF](https://eprint.iacr.org/2021/1063.pdf)), [DEEP Method Medium article](https://medium.com/starkware/starkdex-deep-dive-the-stark-core-engine-497942d0f0ab)
- [Alan Szepieniec](https://twitter.com/aszepieniec?lang=en) - The [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/) tutorial is (IMHO) the best end-to-end practical resource to understand STARKs. Read through this and the magic of STARKs will start to make sense. [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) is a fantastic sequel to the STARK Anatomy tutorial and is a practical guide to creating an AIR that proves programs written in the Brainf**k programming language. Check out miniSTARK's [brainf\*\*k example](examples/brainfuck/) which is an implementation of the BrainSTARK AIR.
- [Winterfell](https://github.com/novifinancial/winterfell) - A STARK prover and verifier developed at Facebook by [Bobbin Threadbare](https://twitter.com/bobbinth) and others. This repo was heavily used as a reference for several components: fast Merkle Tree, DEEP composition polynomial, channels, and STARK component traits. Bobbin has some great HackMD articles: [Miden VM program decoder](https://hackmd.io/_aaDBzbWRz6EwQQRtK1pzw), [Memory in Miden VM](https://hackmd.io/@bobbinth/HJr56BKKt), [u32 operations in Miden VM](https://hackmd.io/NC-yRmmtRQSvToTHb96e8Q#u32-operations-in-Miden-VM).