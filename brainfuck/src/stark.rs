use crate::OpCode;
use algebra::fp_u64::BaseFelt;
use algebra::Felt;
use algebra::PrimeFelt;

fn compile(source: &str) -> Vec<usize> {
    let opcodes = super::lex(source);
    let mut program = Vec::new();
    let mut stack = Vec::new();
    for (i, opcode) in opcodes.into_iter().enumerate() {
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

fn run(program: &[usize], input: &mut impl std::io::Read, output: &mut impl std::io::Write) {
    let mut ip = 0; // instruction pointer
    let mut mp = 0; // memory pointer
    let mut tape = [0u8; 1024];
    let mut running_time = 1;

    while ip < program.len() {
        let instruction = program[ip];
        if instruction == OpCode::LoopBegin.into() {
            if tape[mp] == 0 {
                ip = program[ip + 1];
            } else {
                ip += 2;
            }
        } else if instruction == OpCode::LoopEnd.into() {
            if tape[mp] == 0 {
                ip += 2;
            } else {
                ip = program[ip + 1];
            }
        } else if instruction == OpCode::IncrementPointer.into() {
            ip += 1;
            mp += 1;
        } else if instruction == OpCode::DecrementPointer.into() {
            ip += 1;
            mp -= 1;
        } else if instruction == OpCode::Increment.into() {
            ip += 1;
            tape[mp] += 1;
        } else if instruction == OpCode::Decrement.into() {
            ip += 1;
            tape[mp] -= 1;
        } else if instruction == OpCode::Read.into() {
            ip += 1;
            let mut x = [0u8; 1];
            input.read_exact(&mut x).expect("failed to read input");
            tape[mp] = x[0];
        } else if instruction == OpCode::Write.into() {
            ip += 1;
            output
                .write_all(&tape[mp..mp + 1])
                .expect("failed to write output");
        } else {
            panic!("unrecognized instruction at ip:{}", ip);
        }
        running_time += 1;
    }
}

struct Register {
    cycle: usize,
    ip: usize,
    current_instr: usize,
    next_instr: usize,
    mp: usize,
    memory_value: usize,
}

impl Register {
    fn new() -> Register {
        Self {
            cycle: 0,
            ip: 0,
            current_instr: 0,
            next_instr: 0,
            mp: 0,
            memory_value: 0,
        }
    }
}

struct SimulationMatrices<E> {
    processor: Vec<[E; 7]>,
    instruction: Vec<[E; 3]>,
    input: Vec<[E; 1]>,
    output: Vec<[E; 1]>,
    memory: Vec<[E; 10]>,
}

fn simulate<E: PrimeFelt>(
    program: &[usize],
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
) -> SimulationMatrices<E> {
    let mut register = Register::new();
    register.current_instr = program[0];
    register.next_instr = if program.len() == 1 { 0 } else { program[1] };
    let mut tape = [0u8; 1024];

    // Prepare tables
    let mut matrices = SimulationMatrices::<E> {
        processor: Vec::new(),
        instruction: Vec::new(),
        input: Vec::new(),
        output: Vec::new(),
        memory: Vec::new(),
    };

    for i in 0..program.len() {
        matrices.instruction.push([
            i.into(),
            program[i].into(),
            program.get(i + 1).map_or(0, |&x| x).into(),
        ])
    }

    // main loop
    while register.ip < program.len() {
        let memory_value = register.memory_value.into();

        matrices.processor.push([
            register.cycle.into(),
            register.ip.into(),
            register.current_instr.into(),
            register.next_instr.into(),
            register.mp.into(),
            memory_value,
            memory_value.inverse().unwrap(),
        ]);

        matrices.instruction.push([
            register.ip.into(),
            register.current_instr.into(),
            register.next_instr.into(),
        ]);

        // Update pointer registers according to instruction
        if register.current_instr == OpCode::LoopBegin.into() {
            register.ip = if register.mp == 0 {
                program[register.ip + 1]
            } else {
                register.ip + 2
            };
        } else if register.current_instr == OpCode::LoopEnd.into() {
            register.ip = if register.memory_value != 0 {
                program[register.ip + 1]
            } else {
                register.ip + 2
            }
        } else if register.current_instr == OpCode::DecrementPointer.into() {
            register.ip += 1;
            register.mp -= 1;
        } else if register.current_instr == OpCode::IncrementPointer.into() {
            register.ip += 1;
            register.mp += 1;
        } else if register.current_instr == OpCode::Increment.into() {
            register.ip += 1;
            tape[register.mp] += 1;
        } else if register.current_instr == OpCode::Decrement.into() {
            register.ip += 1;
            tape[register.mp] -= 1;
        } else if register.current_instr == OpCode::Write.into() {
            register.ip += 1;
            let x = &tape[register.mp..register.mp + 1];
            output.write_all(x).expect("failed to write output");
            matrices.output.push([x[0].into()]);
        } else if register.current_instr == OpCode::Read.into() {
            register.ip += 1;
            let mut x = [0u8; 1];
            input.read_exact(&mut x).expect("failed to read input");
            tape[register.mp] = x[0];
            matrices.input.push([x[0].into()])
        } else {
            panic!("unrecognized instruction at ip:{}", register.ip);
        }

        register.cycle += 1;
        register.current_instr = program.get(register.ip).map_or(0, |&x| x);
        register.next_instr = program.get(register.ip + 1).map_or(0, |&x| x);
        register.memory_value = tape[register.mp].into(); // TODO: Change to u8
    }

    // Collect final state into execution tables
    let memory_value = register.memory_value.into();
    matrices.processor.push([
        register.cycle.into(),
        register.ip.into(),
        register.current_instr.into(),
        register.next_instr.into(),
        register.mp.into(),
        memory_value,
        memory_value.inverse().unwrap(),
    ]);

    matrices.instruction.push([
        register.ip.into(),
        register.current_instr.into(),
        register.next_instr.into(),
    ]);

    // sort instructions by address
    matrices.instruction.sort_by_key(|row| row[0].into_bigint());

    matrices.mem
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
