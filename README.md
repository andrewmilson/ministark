# MiniSTARK - GPU accelerated STARK engine

My time working at Immutable left me fascinated with theoretical proof systems and their applications. I developed an all-consuming amateur interest for scalable transparent arguments of knowledge (STARKs). They give the ability to prove to a 3rd party that some computation with $n$ steps ran correctly. Proofs can be verified in $O(log^{2}(n))$ vs naive $O(n)$<sup>1</sup>.

The library is primarily written in Rust but has been optimised by moving heavily parellisable polynomial arithmetic from the CPU to the GPU. Rust has poor support for coding GPUs so [Apple's metal shader language](https://developer.apple.com/metal/) was used instead. The Rust code is powered by [arkworks](https://github.com/arkworks-rs) and [winterfell](https://github.com/novifinancial/winterfell) was heavily used as a reference for several components: fast Merkle Tree, DEEP composition polynomial, FRI, channels, overall design. This library runs around twice as fast as Winterfell at the moment but will evolve to be even faster and have lower memory requirements.

MiniSTARK is unique with how constraints are expressed over the execution trace. Constraints are represented symbolically as a multivariate polynomial where each variable represents a column of the execution trace ([Column](mini-stark/src/constraint.rs) trait) or one of the verifier's challenges ([Challenge](mini-stark/src/constraint.rs) trait). Bellow is a contrived example to illustrate how constraints are represented. Challenges are represented in a similar way. Look at the [Brainf***](mini-stark/examples/brainfuck/) for a full example.

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

// cycle starts at zero
let constraint = Cycle.first().equals(F::zero());

// cycle increases from one row to the next
let constraint = Cycle.curr() - Cycle.next() - F::one();

// cycle doesn't exceed expected
let constraint = Cycle.last().equals(ExpectedCycles.get_hint());
```

## Examples

### Brainf*** virtual machine

Implementation of the [Brainf***](https://esolangs.org/wiki/Brainfuck) virtual machine from [Alan Szepieniec BrainSTARK tutorial](https://aszepieniec.github.io/stark-brainfuck/brainfuck).

```bash
# source: https://esolangs.org
export HELLO_WORLD_BF_SOURCE="++++++++[>++++[>++>+++>+++>+<<<<-]>+>+>->>+[<]<-]>>.>---.+++++++..+++.>>.<-.<.+++.------.--------.>>+.>++."
cargo run --release --features parallel,asm  --example bf --src $HELLO_WORLD_BF_SOURCE
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

- debugging memory table. Remove constraint for memory stay the same clock cycle increase

## Help wanted

- More tests
- More benchmarks
- Making gpu-poly less `unsafe`
- CUDA support in `gpu-poly`
- Speed optimizations that don't sacrifice readability
- Memory reductions that don't sacrifice readability
- Using more `arkworks` features where appropriate
- Enabling zero knowledge

## Reference and thanks

- [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/) by [Alan Szepieniec](https://twitter.com/aszepieniec?lang=en) - Probably the best practical resource I've read on STARKs to date. STARKs remained a mystery to me until I went through this word by word (it took a while to digest). Read through this and the magic of STARKs will start to make sense.
- [BrainSTARK](https://aszepieniec.github.io/stark-brainfuck/brainfuck) by [Alan Szepieniec](https://twitter.com/aszepieniec?lang=en) - A practical STARK tutorial to create a STARK vm for the BF programming language. Is a sequel to Anatomy of a STARK and runs through table-lookup arguments
- TODO: winterfel

<sup>1</sup> TODO: with small error