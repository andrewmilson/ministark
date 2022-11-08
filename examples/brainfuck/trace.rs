use crate::tables::BrainfuckColumn;
use crate::tables::Challenge;
use crate::tables::InputBaseColumn;
use crate::tables::InputExtensionColumn;
use crate::tables::InstructionBaseColumn;
use crate::tables::InstructionExtensionColumn;
use crate::tables::MemoryBaseColumn;
use crate::tables::MemoryExtensionColumn;
use crate::tables::OutputBaseColumn;
use crate::tables::OutputExtensionColumn;
use crate::tables::ProcessorBaseColumn;
use crate::tables::ProcessorExtensionColumn;
use crate::vm::OpCode;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::UniformRand;
use ark_ff::Zero;
use ark_std::rand;
use gpu_poly::allocator::PageAlignedAllocator;
use gpu_poly::fields::p18446744069414584321::Fp;
use gpu_poly::fields::p18446744069414584321::Fq3;
use gpu_poly::GpuVec;
use ministark::challenges::Challenges;
// use ministark::constraint::Challenge as _;
use ministark::Matrix;
use ministark::Trace;

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
        let base_trace = Matrix::join(vec![
            processor_base_trace.clone(),
            memory_base_trace.clone(),
            instruction_base_trace.clone(),
            input_base_trace.clone(),
            output_base_trace.clone(),
        ]);
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
    type Fq = Fq3;

    const NUM_BASE_COLUMNS: usize = 16;
    const NUM_EXTENSION_COLUMNS: usize = 9;

    fn len(&self) -> usize {
        self.base_trace.num_rows()
    }

    fn build_extension_columns(
        &self,
        challenges: &Challenges<Self::Fq>,
    ) -> Option<Matrix<Self::Fq>> {
        let Self {
            processor_base_trace,
            memory_base_trace,
            instruction_base_trace,
            input_base_trace,
            output_base_trace,
            ..
        } = self;

        // TODO: use different method
        let mut rng = rand::thread_rng();
        let instr_initial = Fq3::rand(&mut rng);
        let mem_initial = Fq3::rand(&mut rng);

        let processor_matrix =
            gen_processor_ext_matrix(instr_initial, mem_initial, challenges, processor_base_trace);
        let memory_matrix = gen_memory_ext_matrix(mem_initial, challenges, memory_base_trace);
        let instruction_matrix =
            gen_instruction_ext_matrix(instr_initial, challenges, instruction_base_trace);
        let input_matrix = gen_input_ext_matrix(challenges, input_base_trace);
        let output_matrix = gen_output_ext_matrix(challenges, output_base_trace);

        Some(Matrix::join(vec![
            processor_matrix,
            memory_matrix,
            instruction_matrix,
            input_matrix,
            output_matrix,
        ]))
    }

    fn base_columns(&self) -> &Matrix<Self::Fp> {
        &self.base_trace
    }
}

fn gen_processor_ext_matrix(
    instruction_permutation_initial: Fq3,
    memory_permutation_initial: Fq3,
    challenges: &Challenges<Fq3>,
    base_matrix: &Matrix<Fp>,
) -> Matrix<Fq3> {
    use Challenge::*;
    use ProcessorBaseColumn::*;
    use ProcessorExtensionColumn::*;

    // prepare
    let mut instr_permutation_running_product = instruction_permutation_initial;
    let mut mem_permutation_running_product = memory_permutation_initial;
    let mut input_running_evaluation = Fq3::zero();
    let mut output_running_evaluation = Fq3::zero();

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let curr_base_row = base_matrix.get_row(row).unwrap();
        let next_base_row = base_matrix.get_row(row + 1);
        let mut extension_row = vec![Fq3::zero(); ProcessorExtensionColumn::NUM_TRACE_COLUMNS];

        // Permutations columns
        extension_row[InstructionPermutation as usize] = instr_permutation_running_product;
        extension_row[MemoryPermutation as usize] = mem_permutation_running_product;
        // if not padding
        if !curr_base_row[CurrInstr as usize].is_zero() {
            instr_permutation_running_product *= challenges[Alpha]
                - challenges[A] * curr_base_row[Ip as usize]
                - challenges[B] * curr_base_row[CurrInstr as usize]
                - challenges[C] * curr_base_row[NextInstr as usize];
            mem_permutation_running_product *= challenges[Beta]
                - challenges[D] * curr_base_row[Cycle as usize]
                - challenges[E] * curr_base_row[Mp as usize]
                - challenges[F] * curr_base_row[MemVal as usize];
        }

        // Evaluation columns
        extension_row[InputEvaluation as usize] = input_running_evaluation;
        extension_row[OutputEvaluation as usize] = output_running_evaluation;
        let curr_instr = curr_base_row[CurrInstr as usize].into_bigint().0[0];
        if curr_instr == OpCode::Read as u64 {
            let input_val = next_base_row.unwrap()[MemVal as usize];
            input_running_evaluation = input_running_evaluation * challenges[Gamma] + &input_val;
        } else if curr_instr == OpCode::Write as u64 {
            let output_val = next_base_row.unwrap()[MemVal as usize];
            output_running_evaluation = output_running_evaluation * challenges[Delta] + &output_val;
        }

        extension_rows.push(extension_row);
    }

    // TODO:
    // self.extended_matrix = Some(extended_matrix);
    // self.instr_permutation_terminal = Some(instr_permutation_running_product);
    // self.memory_permutation_terminal = Some(mem_permutation_running_product);
    // self.input_evaluation_terminal = Some(input_running_evaluation);
    // self.output_evaluation_terminal = Some(output_running_evaluation);
    Matrix::new(into_columns(extension_rows))
}

fn gen_memory_ext_matrix(
    memory_permutation_initial: Fq3,
    challenges: &Challenges<Fq3>,
    base_matrix: &Matrix<Fp>,
) -> Matrix<Fq3> {
    use Challenge::*;
    use MemoryBaseColumn::*;
    use MemoryExtensionColumn::*;

    // prepare
    let mut mem_permutation_running_product = memory_permutation_initial;

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let base_row: Vec<Fp> = base_matrix.iter().map(|column| column[row]).collect();
        let mut extension_row = vec![Fq3::zero(); MemoryExtensionColumn::NUM_TRACE_COLUMNS];
        extension_row[Permutation as usize] = mem_permutation_running_product;
        if base_row[Dummy as usize].is_zero() {
            mem_permutation_running_product *= challenges[Beta]
                - challenges[D] * base_row[Cycle as usize]
                - challenges[E] * base_row[Mp as usize]
                - challenges[F] * base_row[MemVal as usize];
        }
        extension_rows.push(extension_row);
    }

    Matrix::new(into_columns(extension_rows))
}

fn gen_instruction_ext_matrix(
    instruction_permutation_initial: Fq3,
    challenges: &Challenges<Fq3>,
    base_matrix: &Matrix<Fp>,
) -> Matrix<Fq3> {
    use Challenge::*;
    use InstructionBaseColumn::*;
    use InstructionExtensionColumn::*;

    // prepare
    let mut permutation_running_product = instruction_permutation_initial;
    let mut evaluation_running_sum = Fq3::zero();
    let mut previous_address = -Fp::one();

    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let curr_base_row = base_matrix.get_row(row).unwrap();
        let prev_base_row = base_matrix.get_row(row.wrapping_sub(1));
        let mut extension_row = vec![Fq3::zero(); InstructionExtensionColumn::NUM_TRACE_COLUMNS];

        if !curr_base_row[CurrInstr as usize].is_zero()
            && row > 0
            && curr_base_row[Ip as usize] == prev_base_row.unwrap()[Ip as usize]
        {
            // permutation argument
            // update running product
            // make sure new row is not padding
            // and that the instruction address didn't just change
            permutation_running_product *= challenges[Alpha]
                - challenges[A] * curr_base_row[Ip as usize]
                - challenges[B] * curr_base_row[CurrInstr as usize]
                - challenges[C] * curr_base_row[NextInstr as usize];
        }
        extension_row[ProcessorPermutation as usize] = permutation_running_product;

        // evaluation argument
        if curr_base_row[Ip as usize] != previous_address {
            evaluation_running_sum = challenges[Eta] * evaluation_running_sum
                + challenges[A] * curr_base_row[Ip as usize]
                + challenges[B] * curr_base_row[CurrInstr as usize]
                + challenges[C] * curr_base_row[NextInstr as usize];
        }
        extension_row[ProgramEvaluation as usize] = evaluation_running_sum;

        previous_address = curr_base_row[Ip as usize];
        extension_rows.push(extension_row);
    }

    Matrix::new(into_columns(extension_rows))
}

fn gen_input_ext_matrix(challenges: &Challenges<Fq3>, base_matrix: &Matrix<Fp>) -> Matrix<Fq3> {
    use Challenge::*;
    use InputBaseColumn::*;
    use InputExtensionColumn::*;

    // prepare
    let mut running_evaluation = Fq3::zero();

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let base_row = base_matrix.get_row(row).unwrap();
        let mut extension_row = vec![Fq3::zero(); InputExtensionColumn::NUM_TRACE_COLUMNS];
        running_evaluation = running_evaluation * challenges[Gamma] + &base_row[Value as usize];
        extension_row[Evaluation as usize] = running_evaluation;
        // TODO:
        // if !self.len().is_zero() && i == self.len() - 1 {
        //     evaluation_terminal = running_evaluation;
        // }
        extension_rows.push(extension_row);
    }

    Matrix::new(into_columns(extension_rows))
}

fn gen_output_ext_matrix(challenges: &Challenges<Fq3>, base_matrix: &Matrix<Fp>) -> Matrix<Fq3> {
    use Challenge::*;
    use OutputBaseColumn::*;
    use OutputExtensionColumn::*;

    // prepare
    let mut running_evaluation = Fq3::zero();
    // let mut evaluation_terminal = Fq3::zero();

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let base_row = base_matrix.get_row(row).unwrap();
        let mut extension_row = vec![Fq3::zero(); OutputExtensionColumn::NUM_TRACE_COLUMNS];
        running_evaluation = running_evaluation * challenges[Delta] + &base_row[Value as usize];
        extension_row[Evaluation as usize] = running_evaluation;
        // TODO:
        // if !self.len().is_zero() && i == self.len() - 1 {
        //     evaluation_terminal = running_evaluation;
        // }
        extension_rows.push(extension_row);
    }

    Matrix::new(into_columns(extension_rows))
}

pub fn into_columns<F: Field>(rows: Vec<Vec<F>>) -> Vec<GpuVec<F>> {
    if rows.is_empty() {
        Vec::new()
    } else {
        let num_cols = rows[0].len();
        let mut cols = Vec::new();
        for _ in 0..num_cols {
            cols.push(Vec::new_in(PageAlignedAllocator));
        }
        for row in rows {
            for (col, val) in cols.iter_mut().zip(row) {
                col.push(val);
            }
        }
        cols
    }
}
