use crate::OpCode;
use algebra::fp_u64::BaseFelt;
use algebra::Felt;

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
