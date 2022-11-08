> **Warning**
> this code is new and will change in future versions. Please try it out and provide feedback. If it addresses a use-case that is important to you please open an issue to discuss it or register your interest https://forms.gle/SYrA9hbGQuivbHX5A

<div align="center">

# Mini<i>STARK</i> <sup>✨</sup>

**GPU accelerated STARK prover and verifier**

[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/andrewmilson/mini-stark/blob/main/LICENSE)
[![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#experimental)

</div>

MiniSTARK allows you to prove the execution of arbitrary computations to anyone using the power of [STARKs](https://starkware.co/stark/). In the example above MiniSTARK is being used to prove "Hello World!" in the [Brainf**k](https://esolangs.org/wiki/Brainfuck) programming language. Anyone can verify the proof instantly.

The library is written Rust but [Metal](https://developer.apple.com/metal/) is used to GPU accelerate polynomial arithmetic. The Rust code is built on [arkworks](https://github.com/arkworks-rs) and is influenced by [Winterfell](https://github.com/novifinancial/winterfell). MiniSTARK has plenty of opportunities for optimizations. If zero knowledge, optimizations or Eliptic Curve FFTs get you going then check out what's [coming soon](#coming-soon).

## Demo

### Proving Hello World in brainf***

TODO

### Proving Fibonacci sequence

## Performance


## Defining AIR Constraints

MiniSTARK represents constraints as multivariate polynomials where each variable abstractly represents a column of the execution trace or one of the verifier's challenges. Bellow is a contrived example to illustrate how constraints are represented. Look at the constraint implementation of [brainf***](examples/brainfuck/) for a full example.

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

- Support Bitcoin's fields with ECFFT (see [wborgeaud/ecfft-bn254](https://github.com/wborgeaud/ecfft-bn254))
- CUDA support in `gpu-poly`
- Periodic transition constraints
- Speed optimizations that don't sacrifice readability
- Memory reductions that don't sacrifice readability
- Using more `arkworks` features where appropriate
- Reduce proof size with batched Merkle proofs
- Cairo VM prover implemented using MiniSTARK (similar to [giza](https://github.com/maxgillett/giza))
- More tests and benchmarks
- More GPU field implementations
- Making gpu-poly less unsafe
- Zero knowledge 

## Thank yous

- [StarkWare](https://starkware.co/) - The company started by the people who invented STARKs. Learnt so many things from the people at this company. Love how public and open their educational material is. Just a few of the things I've enjoyed: GOAT [Eli Ben-Sasson's](https://twitter.com/EliBenSasson) screencasts, Cairo whitepaper ([Video](https://www.youtube.com/watch?v=DTVn0oYLVsE), [PDF](https://eprint.iacr.org/2021/1063.pdf)), [STARK 101](https://starkware.co/stark-101/), [DEEP Method Medium article](https://medium.com/starkware/starkdex-deep-dive-the-stark-core-engine-497942d0f0ab)
- [Alan Szepieniec](https://twitter.com/aszepieniec?lang=en) - The [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/) was probably the best practical resource I've read on STARKs to date. STARKs remained a mystery to me until I went through this word by word (it took a while to digest). Read through this and the magic of STARKs will start to make sense. [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) is a fantastic sequel to the STARK Anatomy series and is a practical guide to creating an AIR that proves programs written in the Brainf*** programming language.
- [Winterfell](https://github.com/novifinancial/winterfell) - A STARK prover and verifier developed at Facebook by [Bobbin Threadbare](https://twitter.com/bobbinth) and others. This repo was heavily used as a reference for several components: fast Merkle Tree, DEEP composition polynomial, channels, and STARK component traits. Bobbin's HackMD articles are great as well: [Miden VM program decoder](https://hackmd.io/_aaDBzbWRz6EwQQRtK1pzw), [Memory in Miden VM](https://hackmd.io/@bobbinth/HJr56BKKt), [u32 operations in Miden VM](https://hackmd.io/NC-yRmmtRQSvToTHb96e8Q#u32-operations-in-Miden-VM).