#![feature(allocator_api)]

use air::BrainfuckAir;
use air::ExecutionInfo;
use ark_ff::One;
use gpu_poly::fields::p18446744069414584321::Fp;
use gpu_poly::fields::p18446744069414584321::Fq3;
use ministark::ProofOptions;
use ministark::Prover;
use std::time::Instant;
use trace::BrainfuckTrace;
use vm::compile;
use vm::simulate;

mod air;
mod constraints;
mod tables;
mod trace;
mod vm;

/// Source: http://esoteric.sange.fi/brainfuck/bf-source/prog/fibonacci.txt
const _FIB_TO_55_SOURCE: &str = "
This determines how many numbers to generate:
    +++++++++++

Program:
    >+>>>>++++++++++++++++++++++++++++++++++++++++++++
    >++++++++++++++++++++++++++++++++<<<<<<[>[>>>>>>+>
    +<<<<<<<-]>>>>>>>[<<<<<<<+>>>>>>>-]<[>++++++++++[-
    <-[>>+>+<<<-]>>>[<<<+>>>-]+<[>[-]<[-]]>[<<[>>>+<<<
    -]>>[-]]<<]>>>[>>+>+<<<-]>>>[<<<+>>>-]+<[>[-]<[-]]
    >[<<+>>[-]]<<<<<<<]>>>>>[+++++++++++++++++++++++++
    +++++++++++++++++++++++.[-]]++++++++++<[->-<]>++++
    ++++++++++++++++++++++++++++++++++++++++++++.[-]<<
    <<<<<<<<<<[>>>+>+<<<<-]>>>>[<<<<+>>>>-]<-[>>.>.<<<
    [-]]<<[>>+>+<<<-]>>>[<<<+>>>-]<<[<+>-]>[<+>-]<<<-]
";

/// Source: https://esolangs.org/wiki/Brainfuck
const HELLO_WORLD_SOURCE: &str = "
+++++ +++++             initialize counter (cell #0) to 10
[                       use loop to set 70/100/30/10
    > +++++ ++              add  7 to cell #1
    > +++++ +++++           add 10 to cell #2
    > +++                   add  3 to cell #3
    > +                     add  1 to cell #4
<<<< -                  decrement counter (cell #0)
]
> ++ .                  print 'H'
> + .                   print 'e'
+++++ ++ .              print 'l'
.                       print 'l'
+++ .                   print 'o'
> ++ .                  print ' '
<< +++++ +++++ +++++ .  print 'W'
> .                     print 'o'
+++ .                   print 'r'
----- - .               print 'l'
----- --- .             print 'd'
> + .                   print '!'
> .                     print '\n'
";

const PROVER_SOURCE_CODE: &str = HELLO_WORLD_SOURCE;

struct BrainfuckProver(ProofOptions);

impl Prover for BrainfuckProver {
    type Fp = Fp;
    type Fq = Fq3;
    type Air = BrainfuckAir;
    type Trace = BrainfuckTrace;

    fn new(options: ProofOptions) -> Self {
        BrainfuckProver(options)
    }

    fn options(&self) -> ProofOptions {
        self.0
    }

    fn get_pub_inputs(&self, trace: &BrainfuckTrace) -> ExecutionInfo {
        ExecutionInfo {
            source_code: PROVER_SOURCE_CODE.to_string(),
            input: trace.input_symbols().to_vec(),
            output: trace.output_symbols().to_vec(),
        }
    }
}

fn main() {
    println!("{:?}", Fp::one());

    let now = Instant::now();
    let program = compile(PROVER_SOURCE_CODE);
    let mut output = Vec::new();
    let trace = simulate(&program, &mut std::io::empty(), &mut output);
    println!("Output: {}", String::from_utf8(output).unwrap());

    let options = ProofOptions::new(32, 16, 8, 8, 64);
    let prover = BrainfuckProver::new(options);
    let proof = prover.generate_proof(trace);
    println!("Runtime: {:?}", now.elapsed());
    let proof = proof.unwrap();
    // let mut proof_bytes = Vec::new();
    //     .serialize_compressed(&mut proof_bytes)
    //     .unwrap();
    // println!("Result: {:?}kb", proof_bytes.len() / 1024);
    proof.verify().unwrap();
}
