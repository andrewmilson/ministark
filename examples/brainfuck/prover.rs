use crate::air::BrainfuckAirConfig;
use crate::air::ExecutionInfo;
use crate::trace::BrainfuckTrace;
use ministark::ProofOptions;
use ministark::Prover;
use ministark_gpu::fields::p18446744069414584321::ark::Fp;
use ministark_gpu::fields::p18446744069414584321::ark::Fq3;

pub struct BrainfuckProver(ProofOptions);

impl Prover for BrainfuckProver {
    type Fp = Fp;
    type Fq = Fq3;
    type AirConfig = BrainfuckAirConfig;
    type Trace = BrainfuckTrace;

    fn new(options: ProofOptions) -> Self {
        BrainfuckProver(options)
    }

    fn options(&self) -> ProofOptions {
        self.0
    }

    fn get_pub_inputs(&self, trace: &BrainfuckTrace) -> ExecutionInfo {
        let meta = trace.meta();
        ExecutionInfo {
            source_code: meta.source_code.to_string(),
            input: meta.input.to_vec(),
            output: meta.output.to_vec(),
        }
    }
}
