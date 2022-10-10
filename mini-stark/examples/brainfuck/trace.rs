use crate::tables::Challenge;
use ark_ff_optimized::fp64::Fp;
use fast_poly::allocator::PageAlignedAllocator;
use mini_stark::challenges::Challenges;
// use mini_stark::constraint::Challenge as _;
use mini_stark::Matrix;
use mini_stark::Trace;

pub struct BrainfuckTrace {
    processor_base_trace: Matrix<Fp>,
    memory_base_trace: Matrix<Fp>,
    instruction_base_trace: Matrix<Fp>,
    input_base_trace: Matrix<Fp>,
    output_base_trace: Matrix<Fp>,
    base_trace: Matrix<Fp>,
}

impl BrainfuckTrace {
    pub fn new(
        processor_base_trace: Matrix<Fp>,
        memory_base_trace: Matrix<Fp>,
        instruction_base_trace: Matrix<Fp>,
        input_base_trace: Matrix<Fp>,
        output_base_trace: Matrix<Fp>,
    ) -> Self {
        let mut base_trace = processor_base_trace.clone();
        base_trace.append(memory_base_trace.clone());
        base_trace.append(instruction_base_trace.clone());
        base_trace.append(input_base_trace.clone());
        base_trace.append(output_base_trace.clone());
        BrainfuckTrace {
            processor_base_trace,
            memory_base_trace,
            instruction_base_trace,
            input_base_trace,
            output_base_trace,
            base_trace,
        }
    }
}

impl Trace for BrainfuckTrace {
    type Fp = Fp;

    const NUM_BASE_COLUMNS: usize = 35;

    fn len(&self) -> usize {
        println!("YOOO {}", self.base_trace.num_rows());
        self.base_trace.num_rows()
    }

    fn build_extension_columns(
        &self,
        challenges: &Challenges<Self::Fp>,
    ) -> Option<Matrix<Self::Fp>> {
        let Self {
            processor_base_trace,
            memory_base_trace,
            instruction_base_trace,
            input_base_trace,
            output_base_trace,
            ..
        } = self;

        let processor_matrix = gen_processor_ext_matrix(challenges, processor_base_trace);
        let memory_matrix = gen_memory_ext_matrix(challenges, memory_base_trace);
        let instruction_matrix = gen_instruction_ext_matrix(challenges, instruction_base_trace);
        let input_matrix = gen_input_ext_matrix(challenges, input_base_trace);
        let output_matrix = gen_output_ext_matrix(challenges, output_base_trace);

        let mut extension_matrix = processor_matrix;
        extension_matrix.append(memory_matrix);
        extension_matrix.append(instruction_matrix);
        extension_matrix.append(input_matrix);
        extension_matrix.append(output_matrix);

        Some(extension_matrix)
    }

    fn base_columns(&self) -> &Matrix<Self::Fp> {
        &self.base_trace
    }
}

fn gen_processor_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    todo!()
}

fn gen_memory_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    todo!()
}

fn gen_instruction_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    todo!()
}

fn gen_input_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    todo!()
}

fn gen_output_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    todo!()
}
