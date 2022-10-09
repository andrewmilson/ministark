use crate::tables::Memory;
use crate::tables::Processor;
use fast_poly::allocator::PageAlignedAllocator;
use fast_poly::GpuField;
use mini_stark::Column;
use mini_stark::Matrix;
use std::os::unix::process;

/// Opcodes determined by the lexer
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OpCode {
    IncrementPointer,
    DecrementPointer,
    Increment,
    Decrement,
    Write,
    Read,
    LoopBegin,
    LoopEnd,
}

impl OpCode {
    pub fn iterator() -> std::slice::Iter<'static, OpCode> {
        static VALUES: [OpCode; 8] = [
            OpCode::IncrementPointer,
            OpCode::DecrementPointer,
            OpCode::Increment,
            OpCode::Decrement,
            OpCode::Write,
            OpCode::Read,
            OpCode::LoopBegin,
            OpCode::LoopEnd,
        ];
        VALUES.iter()
    }
}

// TODO: might remove
impl std::convert::Into<usize> for OpCode {
    fn into(self) -> usize {
        Into::<u64>::into(self) as usize
    }
}

impl std::convert::Into<u64> for OpCode {
    fn into(self) -> u64 {
        match self {
            OpCode::IncrementPointer => b'>'.into(),
            OpCode::DecrementPointer => b'<'.into(),
            OpCode::Increment => b'+'.into(),
            OpCode::Decrement => b'-'.into(),
            OpCode::Write => b'.'.into(),
            OpCode::Read => b','.into(),
            OpCode::LoopBegin => b'['.into(),
            OpCode::LoopEnd => b']'.into(),
        }
    }
}

/// Lexer turns the source code into a sequence of opcodes
fn lex(source: &str) -> Vec<OpCode> {
    let mut operations = Vec::new();

    for symbol in source.chars() {
        let op = match symbol {
            '>' => Some(OpCode::IncrementPointer),
            '<' => Some(OpCode::DecrementPointer),
            '+' => Some(OpCode::Increment),
            '-' => Some(OpCode::Decrement),
            '.' => Some(OpCode::Write),
            ',' => Some(OpCode::Read),
            '[' => Some(OpCode::LoopBegin),
            ']' => Some(OpCode::LoopEnd),
            _ => None,
        };

        // Non-opcode characters are comments
        match op {
            Some(op) => operations.push(op),
            None => (),
        }
    }

    operations
}

pub fn compile(source: &str) -> Vec<usize> {
    let opcodes = lex(source);
    let mut program = Vec::new();
    let mut stack = Vec::new();
    for opcode in opcodes.into_iter() {
        program.push(opcode.clone().into());
        match opcode {
            OpCode::LoopBegin => {
                // Placeholder for position of loop end
                program.push(0);
                stack.push(program.len() - 1);
            }
            OpCode::LoopEnd => {
                let last = stack.pop().expect("loop has no beginning");
                program.push(last + 1); // loop end
                program[last] = program.len(); // loop beginning
            }
            _ => (),
        }
    }
    program
}

/// Registers of the brainfuck VM
#[derive(Default)]
struct Register {
    /// Cycle
    cycle: usize,
    /// Instruction pointer
    ip: usize,
    /// Current instruction
    curr_instr: usize,
    /// Next instruction
    next_instr: usize,
    /// Memory pointer
    mp: usize,
    /// Memory value
    mem_val: usize,
}

// Outputs base execution trace
pub fn simulate<F: GpuField>(
    program: &[usize],
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
) -> Matrix<F> {
    let mut tape = [0u8; 1024];
    let mut register = Register::default();
    register.curr_instr = program[0];
    register.next_instr = if program.len() == 1 { 0 } else { program[1] };

    // execution trace tables in row major
    let mut processor_rows = Vec::new();
    let mut instruction_rows = Vec::new();
    let mut input_rows = Vec::new();
    let mut output_rows = Vec::new();
    let mut memory_rows = Vec::new();

    for i in 0..program.len() {
        instruction_rows.push(vec![
            F::from(i as u64),
            F::from(program[i] as u64),
            F::from(program.get(i + 1).map_or(0, |&x| x as u64)),
        ])
    }

    // main loop
    while register.ip < program.len() {
        let mem_val = F::from(register.mem_val as u64);

        println!("Cycle: {}", register.cycle);

        processor_rows.push(vec![
            F::from(register.cycle as u64),
            F::from(register.ip as u64),
            F::from(register.curr_instr as u64),
            F::from(register.next_instr as u64),
            F::from(register.mp as u64),
            mem_val,
            mem_val.inverse().unwrap_or_else(F::zero),
        ]);

        instruction_rows.push(vec![
            F::from(register.ip as u64),
            F::from(register.curr_instr as u64),
            F::from(register.next_instr as u64),
        ]);

        // Update pointer registers according to instruction
        if register.curr_instr == Into::<usize>::into(OpCode::LoopBegin) {
            register.ip = if register.mem_val == 0 {
                program[register.ip + 1]
            } else {
                register.ip + 2
            };
        } else if register.curr_instr == Into::<usize>::into(OpCode::LoopEnd) {
            register.ip = if register.mem_val != 0 {
                program[register.ip + 1]
            } else {
                register.ip + 2
            }
        } else if register.curr_instr == Into::<usize>::into(OpCode::DecrementPointer) {
            register.ip += 1;
            register.mp -= 1;
        } else if register.curr_instr == Into::<usize>::into(OpCode::IncrementPointer) {
            register.ip += 1;
            register.mp += 1;
        } else if register.curr_instr == Into::<usize>::into(OpCode::Increment) {
            register.ip += 1;
            tape[register.mp] += 1;
        } else if register.curr_instr == Into::<usize>::into(OpCode::Decrement) {
            register.ip += 1;
            tape[register.mp] -= 1;
        } else if register.curr_instr == Into::<usize>::into(OpCode::Write) {
            register.ip += 1;
            let x = &tape[register.mp..register.mp + 1];
            output.write_all(x).expect("failed to write output");
            output_rows.push(vec![x[0].into()]);
        } else if register.curr_instr == Into::<usize>::into(OpCode::Read) {
            register.ip += 1;
            let mut x = [0u8; 1];
            input.read_exact(&mut x).expect("failed to read input");
            tape[register.mp] = x[0];
            input_rows.push(vec![x[0].into()])
        } else {
            panic!("unrecognized instruction at ip:{}", register.ip);
        }

        register.cycle += 1;
        register.curr_instr = program.get(register.ip).map_or(0, |&x| x);
        register.next_instr = program.get(register.ip + 1).map_or(0, |&x| x);
        register.mem_val = tape[register.mp].into(); // TODO: Change to u8
    }

    // Collect final state into execution tables
    let mem_val = F::from(register.mem_val as u64);
    processor_rows.push(vec![
        F::from(register.cycle as u64),
        F::from(register.ip as u64),
        F::from(register.curr_instr as u64),
        F::from(register.next_instr as u64),
        F::from(register.mp as u64),
        mem_val,
        mem_val.inverse().unwrap(),
    ]);

    instruction_rows.push(vec![
        F::from(register.ip as u64),
        F::from(register.curr_instr as u64),
        F::from(register.next_instr as u64),
    ]);

    // sort instructions by address
    instruction_rows.sort_by_key(|row| row[0]);

    memory_rows = derive_memory_rows(&processor_rows);

    let padding_len = {
        let max_length = [
            processor_rows.len(),
            memory_rows.len(),
            instruction_rows.len(),
            input_rows.len(),
            output_rows.len(),
        ]
        .into_iter()
        .max()
        .unwrap();
        ceil_power_of_two(max_length)
    };

    let mut columns = vec![
        into_columns(processor_rows),
        into_columns(memory_rows),
        into_columns(instruction_rows),
        into_columns(input_rows),
        into_columns(output_rows),
    ]
    .into_iter()
    .flat_map(|columns| columns)
    .collect::<Vec<Vec<F, PageAlignedAllocator>>>();

    for column in &mut columns {
        column.resize(padding_len, F::zero());
    }

    let mut dummy_column = Vec::with_capacity(padding_len);
    dummy_column.resize(padding_len, F::one());

    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));
    columns.push(dummy_column.to_vec_in(PageAlignedAllocator));

    Matrix::new(columns)
}

fn derive_memory_rows<F: GpuField>(processor_rows: &[Vec<F>]) -> Vec<Vec<F>> {
    // copy unpadded rows and sort
    // TODO: sorted by IP and then CYCLE. Check to see if processor table sorts by
    // cycle.
    let mut memory_rows = processor_rows
        .iter()
        .filter_map(|row| {
            if row[Processor::CurrInstr as usize].is_zero() {
                None
            } else {
                Some(vec![
                    row[Processor::Cycle as usize],
                    row[Processor::Mp as usize],
                    row[Processor::MemVal as usize],
                    F::zero(), // dummy=no
                ])
            }
        })
        .collect::<Vec<Vec<F>>>();

    // matrix.sort_by_key(|row| );
    memory_rows.sort_by_key(|row| (row[Memory::Mp as usize], row[Memory::Cycle as usize]));

    // insert dummy rows for smooth clk jumps
    let mut i = 0;
    while i < memory_rows.len() - 1 {
        let curr_row = &memory_rows[i];
        let next_row = &memory_rows[i + 1];

        // check sorted by memory address then cycle
        if curr_row[Memory::Mp as usize] == next_row[Memory::Mp as usize] {
            // assert!(curr_row[Self::CYCLE] == next_row[Self::CYCLE] -
            // F::BasePrimeField::one())
        }

        if curr_row[Memory::Mp as usize] == next_row[Memory::Mp as usize]
            && curr_row[Memory::Cycle as usize] + F::one() != next_row[Memory::Cycle as usize]
        {
            memory_rows.insert(
                i + 1,
                vec![
                    curr_row[Memory::Cycle as usize] + F::one(),
                    curr_row[Memory::Mp as usize],
                    curr_row[Memory::MemVal as usize],
                    F::one(), // dummy=yes
                ],
            )
        }

        i += 1;
    }

    memory_rows
}

fn into_columns<F: GpuField>(rows: Vec<Vec<F>>) -> Vec<Vec<F, PageAlignedAllocator>> {
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

/// Rounds the input value up the the nearest power of two
fn ceil_power_of_two(value: usize) -> usize {
    if value.is_power_of_two() {
        value
    } else {
        value.next_power_of_two()
    }
}
