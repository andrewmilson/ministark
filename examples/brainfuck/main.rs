#![feature(allocator_api)]

use air::BrainfuckAir;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ministark::Proof;
use ministark::ProofOptions;
use ministark::Prover;
use ministark::Trace;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;
use structopt::StructOpt;
use trace::BrainfuckTrace;
use vm::simulate;

mod air;
mod constraints;
mod prover;
mod tables;
mod trace;
mod vm;

fn prove(options: ProofOptions, source_code_path: PathBuf, output_path: PathBuf) {
    let source_code = fs::read_to_string(source_code_path).unwrap();
    let mut output = Vec::new();

    let now = Instant::now();
    let trace = simulate(source_code, &mut std::io::empty(), &mut output);
    println!(
        "Generated execution trace (cols={}, rows={}) in {:?}",
        trace.base_columns().num_cols(),
        trace.base_columns().num_rows(),
        now.elapsed(),
    );
    println!(
        "Program output: \"{}\"",
        String::from_utf8(output.clone()).unwrap()
    );

    let prover = prover::BrainfuckProver::new(options);
    let proof = prover.generate_proof(trace).unwrap();
    println!("Proof generated in: {:?}", now.elapsed());
    println!("Proof security: {}bit", proof.conjectured_security_level());

    let mut proof_bytes = Vec::new();
    proof.serialize_compressed(&mut proof_bytes).unwrap();
    println!("Proof size: {:?}KB", proof_bytes.len() / 1024);
    let mut f = File::create(&output_path).unwrap();
    f.write_all(proof_bytes.as_slice()).unwrap();
    f.flush().unwrap();
    println!("Proof written to {}", output_path.as_path().display());
}

fn verify(options: ProofOptions, source_code_path: PathBuf, proof_path: PathBuf) {
    let source_code = fs::read_to_string(source_code_path).unwrap();
    let proof_bytes = fs::read(proof_path).unwrap();
    let proof: Proof<BrainfuckAir> = Proof::deserialize_compressed(proof_bytes.as_slice()).unwrap();
    assert_eq!(source_code, proof.public_inputs.source_code);
    assert_eq!(options, proof.options);

    let now = Instant::now();
    proof.verify().unwrap();
    println!("Proof verified in: {:?}", now.elapsed());
}

#[derive(StructOpt, Debug)]
#[structopt(name = "BrainSTARK", about = "miniSTARK brainfuck prover and verifier")]
pub enum BrainfuckOptions {
    Prove {
        #[structopt(long, parse(from_os_str))]
        src: PathBuf,
        #[structopt(long, parse(from_os_str))]
        output: PathBuf,
    },
    Verify {
        #[structopt(long, parse(from_os_str))]
        src: PathBuf,
        #[structopt(long, parse(from_os_str))]
        proof: PathBuf,
    },
}

fn main() {
    // proof options
    let num_queries = 29;
    let lde_blowup_factor = 16;
    let grinding_factor = 16;
    let fri_folding_factor = 8;
    let fri_max_remainder_size = 64;
    let options = ProofOptions::new(
        num_queries,
        lde_blowup_factor,
        grinding_factor,
        fri_folding_factor,
        fri_max_remainder_size,
    );

    // read command-line args
    match BrainfuckOptions::from_args() {
        BrainfuckOptions::Prove { src, output } => prove(options, src, output),
        BrainfuckOptions::Verify { src, proof } => verify(options, src, proof),
    }
}
