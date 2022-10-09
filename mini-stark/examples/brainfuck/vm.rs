use fast_poly::GpuField;
use mini_stark::Matrix;
use std::os::unix::process;

/// Opcodes determined by the lexer
#[derive(Debug, Clone, PartialEq, Eq)]
enum OpCode {
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
    fn iterator() -> std::slice::Iter<'static, OpCode> {
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

pub fn simulate<F: GpuField>(
    program: &[usize],
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
) -> (Matrix<F>, Matrix<F>, Matrix<F>) {
    let mut tape = [0u8; 1024];
    let mut register = Register::default();
    register.curr_instr = program[0];
    register.next_instr = if program.len() == 1 { 0 } else { program[1] };

    // execution trace tables in row major
    let processor_rows = Vec::new();
    let instruction_rows = Vec::new();
    let input_rows = Vec::new();
    let output_rows = Vec::new();
    let memory_rows = Vec::new();

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
    processor_rows.push([
        F::from(register.cycle as u64),
        F::from(register.ip as u64),
        F::from(register.curr_instr as u64),
        F::from(register.next_instr as u64),
        F::from(register.mp as u64),
        mem_val,
        mem_val.inverse().unwrap(),
    ]);

    instruction_rows.push([
        F::from(register.ip as u64),
        F::from(register.curr_instr as u64),
        F::from(register.next_instr as u64),
    ]);

    // sort instructions by address
    matrices.instruction.sort_by_key(|row| row[0].into_bigint());

    matrices.memory = MemoryTable::<F>::derive_matrix(&matrices.processor);

    matrices
}
