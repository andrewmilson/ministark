use crate::memory_table::MemoryTable;
use crate::OpCode;
use ark_ff::FftField;
use ark_ff::PrimeField;

pub fn compile(source: &str) -> Vec<usize> {
    let opcodes = super::lex(source);
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

// Outputs running time
pub fn run(
    program: &[usize],
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
) -> usize {
    let mut ip = 0; // instruction pointer
    let mut mp = 0; // memory pointer
    let mut tape = [0u8; 1024];
    let mut running_time = 1;

    while ip < program.len() {
        let instruction = program[ip];
        if instruction == Into::<usize>::into(OpCode::LoopBegin) {
            if tape[mp] == 0 {
                ip = program[ip + 1];
            } else {
                ip += 2;
            }
        } else if instruction == Into::<usize>::into(OpCode::LoopEnd) {
            if tape[mp] == 0 {
                ip += 2;
            } else {
                ip = program[ip + 1];
            }
        } else if instruction == Into::<usize>::into(OpCode::IncrementPointer) {
            ip += 1;
            mp += 1;
        } else if instruction == Into::<usize>::into(OpCode::DecrementPointer) {
            ip += 1;
            mp -= 1;
        } else if instruction == Into::<usize>::into(OpCode::Increment) {
            ip += 1;
            tape[mp] += 1;
        } else if instruction == Into::<usize>::into(OpCode::Decrement) {
            ip += 1;
            tape[mp] -= 1;
        } else if instruction == Into::<usize>::into(OpCode::Read) {
            ip += 1;
            let mut x = [0u8; 1];
            input.read_exact(&mut x).expect("failed to read input");
            tape[mp] = x[0];
        } else if instruction == Into::<usize>::into(OpCode::Write) {
            ip += 1;
            output
                .write_all(&tape[mp..mp + 1])
                .expect("failed to write output");
        } else {
            panic!("unrecognized instruction at ip:{}", ip);
        }
        running_time += 1;
    }

    running_time
}

struct Register {
    cycle: usize,
    ip: usize,
    curr_instr: usize,
    next_instr: usize,
    mp: usize,
    mem_val: usize,
}

impl Register {
    fn new() -> Register {
        Self {
            cycle: 0,
            ip: 0,
            curr_instr: 0,
            next_instr: 0,
            mp: 0,
            mem_val: 0,
        }
    }
}

pub struct SimulationMatrices<F> {
    pub processor: Vec<[F; 7]>,
    pub instruction: Vec<[F; 3]>,
    pub input: Vec<[F; 1]>,
    pub output: Vec<[F; 1]>,
    pub memory: Vec<[F; 4]>,
}

pub fn simulate<F: FftField + PrimeField>(
    program: &[usize],
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
) -> SimulationMatrices<F> {
    let mut register = Register::new();
    register.curr_instr = program[0];
    register.next_instr = if program.len() == 1 { 0 } else { program[1] };
    let mut tape = [0u8; 1024];

    // Prepare tables
    let mut matrices = SimulationMatrices::<F> {
        processor: Vec::new(),
        instruction: Vec::new(),
        input: Vec::new(),
        output: Vec::new(),
        memory: Vec::new(),
    };

    for i in 0..program.len() {
        matrices.instruction.push([
            F::from(i as u64),
            F::from(program[i] as u64),
            F::from(program.get(i + 1).map_or(0, |&x| x as u64)),
        ])
    }

    // main loop
    while register.ip < program.len() {
        let mem_val = F::from(register.mem_val as u64);

        matrices.processor.push([
            F::from(register.cycle as u64),
            F::from(register.ip as u64),
            F::from(register.curr_instr as u64),
            F::from(register.next_instr as u64),
            F::from(register.mp as u64),
            mem_val,
            mem_val.inverse().unwrap_or_else(F::zero),
        ]);

        matrices.instruction.push([
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
            matrices.output.push([x[0].into()]);
        } else if register.curr_instr == Into::<usize>::into(OpCode::Read) {
            register.ip += 1;
            let mut x = [0u8; 1];
            input.read_exact(&mut x).expect("failed to read input");
            tape[register.mp] = x[0];
            matrices.input.push([x[0].into()])
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
    matrices.processor.push([
        F::from(register.cycle as u64),
        F::from(register.ip as u64),
        F::from(register.curr_instr as u64),
        F::from(register.next_instr as u64),
        F::from(register.mp as u64),
        mem_val,
        mem_val.inverse().unwrap(),
    ]);

    matrices.instruction.push([
        F::from(register.ip as u64),
        F::from(register.curr_instr as u64),
        F::from(register.next_instr as u64),
    ]);

    // sort instructions by address
    matrices.instruction.sort_by_key(|row| row[0].into_bigint());

    matrices.memory = MemoryTable::<F>::derive_matrix(&matrices.processor);

    matrices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compiles() {
        let source = ">>+++[-].,.<<";

        let program = compile(source);

        assert_eq!(
            program,
            [
                62, // >
                62, // >
                43, // +
                43, // +
                43, // +
                91, // [
                10, // j:10
                45, // -
                93, // ]
                7,  // l:7
                46, // .
                44, // ,
                46, // .
                60, // <
                60  // <
            ]
        );
    }

    #[test]
    fn hello_world() {
        let source = "
            +++++ +++++             initialize counter (cell #0) to 10
            [                       use loop to set 70/100/30/10
                > +++++ ++              add  7 to cell #1
                > +++++ +++++           add 10 to cell #2
                > +++                   add  3 to cell #3
                > +                     add  1 to cell #4
            <<<< -                  decrement counter (cell #0)
            ]
            > ++ .                  print 'H'
            > + .                   print 'e'
            +++++ ++ .              print 'l'
            .                       print 'l'
            +++ .                   print 'o'
            > ++ .                  print ' '
            << +++++ +++++ +++++ .  print 'W'
            > .                     print 'o'
            +++ .                   print 'r'
            ----- - .               print 'l'
            ----- --- .             print 'd'
            > + .                   print '!'
            > .                     print '\n'
        ";
        let mut output = Vec::new();
        let program = compile(source);

        run(&program, &mut std::io::empty(), &mut output);

        assert_eq!(output, "Hello World!\n".as_bytes());
    }
}
