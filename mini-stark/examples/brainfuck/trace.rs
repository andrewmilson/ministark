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
use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff::Zero;
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
    use ProcessorBaseColumn::*;
    use ProcessorExtensionColumn::*;

    // TODO: get random initial values
    let instr_permutation_initial = Fp::one();
    let mem_permutation_initial = Fp::one();

    // prepare
    let mut instr_permutation_running_product = instr_permutation_initial;
    let mut mem_permutation_running_product = mem_permutation_initial;
    let mut input_running_evaluation = Fp::zero();
    let mut output_running_evaluation = Fp::zero();

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let curr_base_row = base_matrix
            .iter()
            .map(|column| column[row])
            .collect::<Vec<Fp>>();
        let next_base_row = base_matrix
            .iter()
            .map(|column| column.get(row + 1).copied())
            .collect::<Vec<Option<Fp>>>();
        let mut extension_row = vec![Fp::zero(); ProcessorExtensionColumn::NUM_TRACE_COLUMNS];

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
            let input_val = next_base_row[MemVal as usize].unwrap();
            input_running_evaluation = input_running_evaluation * challenges[Gamma] + input_val;
        } else if curr_instr == OpCode::Write as u64 {
            let output_val = next_base_row[MemVal as usize].unwrap();
            output_running_evaluation = output_running_evaluation * challenges[Delta] + output_val;
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

fn gen_memory_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    use MemoryBaseColumn::*;
    use MemoryExtensionColumn::*;

    // TODO: get random initial values
    let instr_permutation_initial = Fp::one();
    let mem_permutation_initial = Fp::one();

    // prepare
    let mut mem_permutation_running_product = mem_permutation_initial;

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let base_row = base_matrix
            .iter()
            .map(|column| column[row])
            .collect::<Vec<Fp>>();
        let mut extension_row = vec![Fp::zero(); MemoryExtensionColumn::NUM_TRACE_COLUMNS];
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

fn gen_instruction_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    use InstructionBaseColumn::*;
    use InstructionExtensionColumn::*;

    // TODO: get random initial values
    let instr_permutation_initial = Fp::one();
    let mem_permutation_initial = Fp::one();

    // prepare
    let mut permutation_running_product = instr_permutation_initial;
    let mut evaluation_running_sum = Fp::zero();
    let mut previous_address = -Fp::one();

    let mut extension_rows = Vec::new();
    // TODO: remove
    let mut num_padded_rows = 0usize;
    for row in 0..base_matrix.num_rows() {
        let curr_base_row = base_matrix
            .iter()
            .map(|column| column[row])
            .collect::<Vec<Fp>>();
        let prev_base_row = base_matrix
            .iter()
            .map(|column| column.get(row - 1).copied())
            .collect::<Vec<Option<Fp>>>();
        let mut extension_row = vec![Fp::zero(); InstructionExtensionColumn::NUM_TRACE_COLUMNS];

        if curr_base_row[CurrInstr as usize].is_zero() {
            num_padded_rows += 1;
        } else if row > 0 && curr_base_row[Ip as usize] == prev_base_row[Ip as usize].unwrap() {
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

    // self.extended_matrix = Some(extended_matrix);
    // self.permutation_terminal = Some(permutation_running_product);
    // self.evaluation_terminal = Some(evaluation_running_sum);
    Matrix::new(into_columns(extension_rows))
}

fn gen_input_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    use InputBaseColumn::*;
    use InputExtensionColumn::*;

    // prepare
    let mut running_evaluation = Fp::zero();
    let mut evaluation_terminal = Fp::zero();

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let base_row = base_matrix
            .iter()
            .map(|column| column[row])
            .collect::<Vec<Fp>>();
        let mut extension_row = vec![Fp::zero(); InputExtensionColumn::NUM_TRACE_COLUMNS];
        running_evaluation = running_evaluation * challenges[Gamma] + base_row[Value as usize];
        extension_row[Evaluation as usize] = running_evaluation;
        // TODO:
        // if !self.len().is_zero() && i == self.len() - 1 {
        //     evaluation_terminal = running_evaluation;
        // }
        extension_rows.push(extension_row);
    }

    Matrix::new(into_columns(extension_rows))
}

fn gen_output_ext_matrix(challenges: &Challenges<Fp>, base_matrix: &Matrix<Fp>) -> Matrix<Fp> {
    use Challenge::*;
    use OutputBaseColumn::*;
    use OutputExtensionColumn::*;

    // prepare
    let mut running_evaluation = Fp::zero();
    let mut evaluation_terminal = Fp::zero();

    // loop over all rows
    let mut extension_rows = Vec::new();
    for row in 0..base_matrix.num_rows() {
        let base_row = base_matrix
            .iter()
            .map(|column| column[row])
            .collect::<Vec<Fp>>();
        let mut extension_row = vec![Fp::zero(); OutputExtensionColumn::NUM_TRACE_COLUMNS];
        running_evaluation = running_evaluation * challenges[Delta] + base_row[Value as usize];
        extension_row[Evaluation as usize] = running_evaluation;
        // TODO:
        // if !self.len().is_zero() && i == self.len() - 1 {
        //     evaluation_terminal = running_evaluation;
        // }
        extension_rows.push(extension_row);
    }

    Matrix::new(into_columns(extension_rows))
}

pub fn into_columns(rows: Vec<Vec<Fp>>) -> Vec<Vec<Fp, PageAlignedAllocator>> {
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
