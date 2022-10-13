#![feature(allocator_api)]

use ark_ff::One;
use ark_ff_optimized::fp64::Fp;
use fast_poly::allocator::PageAlignedAllocator;
use mini_stark::constraint::are_eq;
use mini_stark::constraint::is_one;
use mini_stark::Air;
use mini_stark::Column;
use mini_stark::Constraint;
use mini_stark::Matrix;
use mini_stark::ProofOptions;
use mini_stark::Prover;
use mini_stark::Trace;
use mini_stark::TraceInfo;
use std::time::Instant;

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
}

struct FibAir {
    options: ProofOptions,
    trace_info: TraceInfo,
    result: Fp,
    boundary_constraints: Vec<Constraint<Fp>>,
    transition_constraints: Vec<Constraint<Fp>>,
    terminal_constraints: Vec<Constraint<Fp>>,
}

impl FibAir {
    fn generate_boundary_constraints() -> Vec<Constraint<Fp>> {
        let v0 = Constraint::from(Fp::one());
        let v1 = &v0 + &v0;
        let v2 = &v0 * &v1;
        let v3 = &v1 * &v2;
        let v4 = &v2 * &v3;
        let v5 = &v3 * &v4;
        let v6 = &v4 * &v5;
        let v7 = &v5 * &v6;

        vec![
            are_eq(0.curr(), v0),
            are_eq(1.curr(), v1),
            are_eq(2.curr(), v2),
            are_eq(3.curr(), v3),
            are_eq(4.curr(), v4),
            are_eq(5.curr(), v5),
            are_eq(6.curr(), v6),
            are_eq(7.curr(), v7),
        ]
    }

    fn generate_transition_constraints() -> Vec<Constraint<Fp>> {
        vec![
            are_eq(0.next(), 6.curr() * 7.curr()),
            are_eq(1.next(), 7.curr() * 0.next()),
            are_eq(2.next(), 0.next() * 1.next()),
            are_eq(3.next(), 1.next() * 2.next()),
            are_eq(4.next(), 2.next() * 3.next()),
            are_eq(5.next(), 3.next() * 4.next()),
            are_eq(6.next(), 4.next() * 5.next()),
            are_eq(7.next(), 5.next() * 6.next()),
        ]
    }

    fn generate_terminal_constraints(result: Fp) -> Vec<Constraint<Fp>> {
        vec![7.curr() - result]
    }
}

impl Air for FibAir {
    type Fp = Fp;
    type PublicInputs = Fp;

    fn new(trace_info: TraceInfo, public_input: Fp, options: ProofOptions) -> Self {
        FibAir {
            options,
            trace_info,
            result: public_input,
            boundary_constraints: Self::generate_boundary_constraints(),
            transition_constraints: Self::generate_transition_constraints(),
            terminal_constraints: Self::generate_terminal_constraints(public_input),
        }
    }

    fn options(&self) -> &ProofOptions {
        &self.options
    }

    fn pub_inputs(&self) -> &Self::PublicInputs {
        &self.result
    }

    fn boundary_constraints(&self) -> &[Constraint<Self::Fp>] {
        &self.boundary_constraints
    }

    fn transition_constraints(&self) -> &[Constraint<Self::Fp>] {
        &self.transition_constraints
    }

    fn terminal_constraints(&self) -> &[Constraint<Self::Fp>] {
        &self.terminal_constraints
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
        let res = *trace.0[7].last().unwrap();
        println!("Pub: {}", res);
        res
    }
}

fn gen_trace(len: usize) -> FibTrace {
    assert!(len.is_power_of_two());
    assert!(len > 8);

    let mut col0 = Vec::new_in(PageAlignedAllocator);
    let mut col1 = Vec::new_in(PageAlignedAllocator);
    let mut col2 = Vec::new_in(PageAlignedAllocator);
    let mut col3 = Vec::new_in(PageAlignedAllocator);
    let mut col4 = Vec::new_in(PageAlignedAllocator);
    let mut col5 = Vec::new_in(PageAlignedAllocator);
    let mut col6 = Vec::new_in(PageAlignedAllocator);
    let mut col7 = Vec::new_in(PageAlignedAllocator);

    let mut v0 = Fp::one();
    let mut v1 = v0 + v0;
    let mut v2 = v0 * v1;
    let mut v3 = v1 * v2;
    let mut v4 = v2 * v3;
    let mut v5 = v3 * v4;
    let mut v6 = v4 * v5;
    let mut v7 = v5 * v6;

    for _ in 0..len / 8 {
        col0.push(v0);
        col1.push(v1);
        col2.push(v2);
        col3.push(v3);
        col4.push(v4);
        col5.push(v5);
        col6.push(v6);
        col7.push(v7);

        v0 = v6 * v7;
        v1 = v7 * v0;
        v2 = v0 * v1;
        v3 = v1 * v2;
        v4 = v2 * v3;
        v5 = v3 * v4;
        v6 = v4 * v5;
        v7 = v5 * v6;
    }

    FibTrace(Matrix::new(vec![
        col0, col1, col2, col3, col4, col5, col6, col7,
    ]))
}

fn main() {
    let now = Instant::now();
    let options = ProofOptions::new(32, 4);
    let prover = FibProver::new(options);
    let trace = gen_trace(1048576 / 4);
    let proof = prover.generate_proof(trace);
    println!("Runtime: {:?}", now.elapsed());
    println!("Result: {:?}", proof.unwrap());
}
