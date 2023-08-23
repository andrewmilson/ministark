#![feature(allocator_api)]

use air::BrainfuckAirConfig;
use ark_serialize::CanonicalDeserialize;
use ark_serialize::CanonicalSerialize;
use ministark::hash::HashFn;
use ministark::hash::Sha256HashFn;
use ministark::merkle::MatrixMerkleTreeImpl;
use ministark::random::PublicCoin;
use ministark::random::PublicCoinImpl;
use ministark::stark::Stark;
use ministark::utils::SerdeOutput;
use ministark::Proof;
use ministark::ProofOptions;
use ministark::Trace;
use ministark_gpu::fields::p18446744069414584321::ark::Fp;
use ministark_gpu::fields::p18446744069414584321::ark::Fq3;
use sha2::Sha256;
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
mod tables;
mod trace;
mod vm;

#[derive(StructOpt, Debug)]
#[structopt(name = "BrainSTARK", about = "miniSTARK brainfuck prover and verifier")]
enum BrainfuckOptions {
    Prove {
        src: PathBuf,
        #[structopt(long, parse(from_os_str))]
        dst: PathBuf,
        #[structopt(long, default_value = "")]
        input: String,
    },
    Verify {
        src: PathBuf,
        #[structopt(long, parse(from_os_str))]
        proof: PathBuf,
        #[structopt(long, default_value = "")]
        input: String,
        #[structopt(long)]
        output: String,
    },
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct BrainfuckClaim {
    pub source_code: String,
    pub input: Vec<u8>,
    pub output: Vec<u8>,
}

impl Stark for BrainfuckClaim {
    type Fp = Fp;
    type Fq = Fq3;
    type AirConfig = BrainfuckAirConfig;
    type Digest = SerdeOutput<Sha256>;
    type PublicCoin = PublicCoinImpl<Fq3, Sha256HashFn>;
    type MerkleTree = MatrixMerkleTreeImpl<Sha256HashFn>;
    type Witness = BrainfuckTrace;
    type Trace = BrainfuckTrace;

    fn gen_public_coin(&self, air: &ministark::Air<Self::AirConfig>) -> Self::PublicCoin {
        let mut seed = Vec::new();
        air.public_inputs().serialize_compressed(&mut seed).unwrap();
        air.trace_len().serialize_compressed(&mut seed).unwrap();
        air.options().serialize_compressed(&mut seed).unwrap();
        PublicCoinImpl::new(Sha256HashFn::hash_chunks([&*seed]))
    }

    fn get_public_inputs(&self) -> Self {
        self.clone()
    }

    fn generate_trace(&self, witness: BrainfuckTrace) -> BrainfuckTrace {
        witness
    }
}

const SECURITY_LEVEL: u32 = 96;

/// Proof options for 96 bit security level
const OPTIONS: ProofOptions = {
    let num_queries = 19;
    let lde_blowup_factor = 16;
    let grinding_factor = 20;
    let fri_folding_factor = 16;
    let fri_max_remainder_coeffs = 16;
    ProofOptions::new(
        num_queries,
        lde_blowup_factor,
        grinding_factor,
        fri_folding_factor,
        fri_max_remainder_coeffs,
    )
};

fn main() {
    // read command-line args
    match BrainfuckOptions::from_args() {
        BrainfuckOptions::Prove { src, dst, input } => prove(src, input, dst),
        BrainfuckOptions::Verify {
            src,
            proof,
            input,
            output,
        } => verify(src, input, output, proof),
    }
}

fn prove(source_code_path: PathBuf, input: String, output_path: PathBuf) {
    let source_code = fs::read_to_string(source_code_path).unwrap();
    let mut output = Vec::new();

    let now = Instant::now();
    let trace = simulate(&source_code, &mut input.as_bytes(), &mut output);
    println!(
        "Generated execution trace (cols={}, rows={}) in {:.0?}",
        trace.base_columns().num_cols(),
        trace.base_columns().num_rows(),
        now.elapsed(),
    );
    println!(
        "Program output: \"{}\"",
        String::from_utf8(output.clone()).unwrap()
    );

    let claim = BrainfuckClaim {
        source_code,
        input: input.into_bytes(),
        output,
    };

    let now = Instant::now();
    let proof = pollster::block_on(claim.prove(OPTIONS, trace)).unwrap();
    println!("Proof generated in: {:.0?}", now.elapsed());
    let security_level = proof.security_level_bits();
    println!("Proof security (conjectured): {security_level}bit",);

    let mut proof_bytes = Vec::new();
    (claim, proof)
        .serialize_compressed(&mut proof_bytes)
        .unwrap();
    println!("Proof size: {:?}KB", proof_bytes.len() / 1024);
    let mut f = File::create(&output_path).unwrap();
    f.write_all(proof_bytes.as_slice()).unwrap();
    f.flush().unwrap();
    println!("Proof written to {}", output_path.as_path().display());
}

fn verify(source_code_path: PathBuf, input: String, output: String, proof_path: PathBuf) {
    let source_code = fs::read_to_string(source_code_path).unwrap();
    let proof_bytes = fs::read(proof_path).unwrap();
    let (execution_info, proof): (BrainfuckClaim, Proof<BrainfuckClaim>) =
        <_>::deserialize_compressed(proof_bytes.as_slice()).unwrap();
    assert_eq!(input.as_bytes(), execution_info.input);
    assert_eq!(output.as_bytes(), execution_info.output);
    assert_eq!(source_code, execution_info.source_code);

    let now = Instant::now();
    execution_info
        .verify(proof, SECURITY_LEVEL)
        .expect("verification failed");
    println!("Proof verified in: {:?}", now.elapsed());
}
