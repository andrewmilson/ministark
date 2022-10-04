#![feature(allocator_api)]

use ark_ff::One;
use ark_ff_optimized::fp64::Fp;
use fast_poly::allocator::PageAlignedAllocator;
use mini_stark::constraint::helper::are_eq;
use mini_stark::constraint::helper::is_one;
use mini_stark::Air;
use mini_stark::Column;
use mini_stark::Constraint;
use mini_stark::Matrix;
use mini_stark::ProofOptions;
use mini_stark::Prover;
use mini_stark::Trace;
use mini_stark::TraceInfo;
use std::time::Instant;

enum Fib {
    FirstCol,
    SecondCol,
}

impl Column for Fib {
    fn index(&self) -> usize {
        match self {
            Self::FirstCol => 0,
            Self::SecondCol => 1,
        }
    }
}

struct FibTrace(Matrix<Fp>);

impl Trace for FibTrace {
    type Fp = Fp;

    const NUM_BASE_COLUMNS: usize = 2;

    fn len(&self) -> usize {
        println!("YOOO {}", self.0.num_rows());
        self.0.num_rows()
    }

    fn base_columns(&self) -> &Matrix<Self::Fp> {
        &self.0
    }

    fn build_extension_columns(&self, challenges: &[Self::Fp]) -> Option<Matrix<Self::Fp>> {
        todo!()
    }
}

struct FibAir {
    options: ProofOptions,
    trace_info: TraceInfo,
    result: Fp,
}

impl Air for FibAir {
    type Fp = Fp;
    type PublicInputs = Fp;

    fn new(trace_info: TraceInfo, public_input: Fp, options: ProofOptions) -> Self {
        FibAir {
            options,
            trace_info,
            result: public_input,
        }
    }

    fn boundary_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        vec![is_one(Fib::FirstCol.curr()), is_one(Fib::SecondCol.curr())]
    }

    fn transition_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        use Fib::*;
        vec![
            are_eq(FirstCol.curr() + SecondCol.curr(), FirstCol.next()),
            are_eq(FirstCol.next() + SecondCol.curr(), SecondCol.next()),
        ]
    }

    fn terminal_constraints(&self) -> Vec<Constraint<Self::Fp>> {
        vec![Fib::SecondCol.curr() - self.result]
    }

    fn trace_info(&self) -> &TraceInfo {
        &self.trace_info
    }
}

struct FibProver(ProofOptions);

impl Prover for FibProver {
    type Fp = Fp;
    type Air = FibAir;
    type Trace = FibTrace;

    fn new(options: ProofOptions) -> Self {
        FibProver(options)
    }

    fn options(&self) -> ProofOptions {
        self.0
    }

    fn get_pub_inputs(&self, trace: &FibTrace) -> <<Self as Prover>::Air as Air>::PublicInputs {
        let res = *trace.0[Fib::SecondCol].last().unwrap();
        println!("Pub: {}", res);
        res
    }
}

fn gen_trace(len: usize) -> FibTrace {
    assert!(len.is_power_of_two());
    assert!(len > 8);
    let mut col1 = Vec::new_in(PageAlignedAllocator);
    let mut col2 = Vec::new_in(PageAlignedAllocator);
    col1.push(Fp::one());
    col2.push(Fp::one());
    for _ in 1..len / 2 {
        col1.push(*col1.last().unwrap() + col2.last().unwrap());
        col2.push(*col1.last().unwrap() + col2.last().unwrap());
    }
    FibTrace(Matrix::new(vec![col1, col2]))
}

fn main() {
    let now = Instant::now();
    let options = ProofOptions::new(32, 4);
    let prover = FibProver::new(options);
    let trace = gen_trace(4194304);
    let proof = prover.generate_proof(trace);
    println!("Runtime: {:?}", now.elapsed());
    println!("Result: {:?}", proof.unwrap());
}
