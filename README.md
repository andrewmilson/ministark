> **Warning**
> this code is new and may change or be removed in future versions. Please try it out and provide feedback. If it addresses a use-case that is important to you please open an issue to discuss it or register your interest https://forms.gle/SYrA9hbGQuivbHX5A

<div align="center">

# Mini<i>STARK</i> ✨

**GPU accelerated STARK prover and verifier**

[![license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/andrewmilson/mini-stark/blob/main/LICENSE)
[![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/mkenney/software-guides/blob/master/STABILITY-BADGES.md#experimental)

</div>

My time working at Immutable left me fascinated with theoretical proof systems and their applications. I developed an all-consuming amateur interest for scalable transparent arguments of knowledge (STARKs). They give the ability to prove to a 3rd party that some computation with $n$ steps ran correctly. Proofs can be verified in $O(log^{2}(n))$ vs naive $O(n)$<sup>1</sup> ([Learn more about STARKS](https://starkware.co/stark/)). 

The library is primarily written in Rust but has been optimized by moving heavily parallelizable polynomial arithmetic from the CPU to the GPU. Rust has poor support for programming GPUs so all the GPU code is written in [Metal](https://developer.apple.com/metal/). The Rust code is powered by [arkworks](https://github.com/arkworks-rs). There is plenty opportunities for increasing speed and reducing memory overhead.

Constraints in MiniSTARK are represented as multivariate polynomial where each variable abstractly represents a column of the execution trace or one of the verifier's challenges. Bellow is a contrived example to illustrate how constraints are represented. Challenges are represented in a similar way. Look at [Brainf***](mini-stark/examples/brainfuck/) for a full example.

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

// all must evaluate to 0 over the trace domain
let constraints = vec![
    // cycle starts at zero
    Cycle.first(),
    // cycle increases from one row to the next
    Cycle.curr() - Cycle.next() - F::one(),
    // cycle doesn't exceed expected
    Cycle.last() - ExpectedCycles.get_hint(),
    // ...
];
```

## Demo

### Proving Hello World in brainf***

TODO

### Proving Fibonacci squence

## Performance




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

## Coming soon (help wanted)

- More tests and benchmarks
- Making gpu-poly less `unsafe`
- CUDA support in `gpu-poly`
- Periodic transition constraints
- Speed optimizations that don't sacrifice readability
- Memory reductions that don't sacrifice readability
- Using more `arkworks` features where appropriate
- Cairo VM prover implemented using MiniSTARK (similar to [giza](https://github.com/maxgillett/giza))
- Zero knowledge 

## Thank yous

- [StarkWare](https://starkware.co/) - The company started by the people who invented STARKs. Learnt so many things from the people at this company. Love how public and open their educational material is. Just a few of the things I've enjoyed: GOAT [Eli Ben-Sasson's](https://twitter.com/EliBenSasson) screencasts, Cairo whitepaper ([Video](https://www.youtube.com/watch?v=DTVn0oYLVsE), [PDF](https://eprint.iacr.org/2021/1063.pdf)), [STARK 101](https://starkware.co/stark-101/), [DEEP Method Medium article](https://medium.com/starkware/starkdex-deep-dive-the-stark-core-engine-497942d0f0ab)
- [Alan Szepieniec](https://twitter.com/aszepieniec?lang=en) - The [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/) was probably the best practical resource I've read on STARKs to date. STARKs remained a mystery to me until I went through this word by word (it took a while to digest). Read through this and the magic of STARKs will start to make sense. [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) is a fantastic sequel to the STARK Anatomy series and is a practical guide to creating an AIR that proves programs written in the Brainf*** programming language.
- [Winterfell](https://github.com/novifinancial/winterfell) - A STARK prover and verifier developed at Facebook by [Bobbin Threadbare](https://twitter.com/bobbinth) and others. This repo was heavily used as a reference for several components: fast Merkle Tree, DEEP composition polynomial, channels, and STARK component traits. Bobbin's HackMD articles are great as well: [Miden VM program decoder](https://hackmd.io/_aaDBzbWRz6EwQQRtK1pzw), [Memory in Miden VM](https://hackmd.io/@bobbinth/HJr56BKKt), [u32 operations in Miden VM](https://hackmd.io/NC-yRmmtRQSvToTHb96e8Q#u32-operations-in-Miden-VM).

<sup>1</sup> with 0.00000000000000000000000000000000000000000001% error (the same probability that a randomly picked atom in the universe belongs to your body)