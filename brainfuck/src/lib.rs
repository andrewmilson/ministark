//! Implementation inspired by https://github.com/Overv/bf
#![feature(generic_const_exprs, int_log)]

use std::slice::Iter;

pub mod evaluation_argument;
mod instruction_table;
mod io_table;
mod memory_table;
pub mod permutation_argument;
mod processor_table;
pub mod stark;
mod table;
mod util;

pub use instruction_table::InstructionTable;
pub use io_table::InputTable;
pub use io_table::OutputTable;
pub use memory_table::MemoryTable;
pub use processor_table::ProcessorTable;
pub use table::Table;

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
    fn iterator() -> Iter<'static, OpCode> {
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

    //     fn into_felt<E: Felt>(self) -> E {
    //         match self {
    //             OpCode::IncrementPointer => E::from(b'>'),
    //             OpCode::DecrementPointer => E::from(b'<'),
    //             OpCode::Increment => E::from(b'+'),
    //             OpCode::Decrement => E::from(b'-'),
    //             OpCode::Write => E::from(b'.'),
    //             OpCode::Read => E::from(b','),
    //             OpCode::LoopBegin => E::from(b'['),
    //             OpCode::LoopEnd => E::from(b']'),
    //         }
    //     }
}

impl std::convert::Into<usize> for OpCode {
    fn into(self) -> usize {
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

#[derive(Debug, Clone)]
enum Instr {
    IncrementPointer,
    DecrementPointer,
    Increment,
    Decrement,
    Write,
    Read,
    Loop(Vec<Instr>),
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

fn parse(opcodes: &[OpCode]) -> Vec<Instr> {
    let mut program = Vec::new();
    let mut loop_stack = 0;
    let mut loop_start = 0;

    for (i, op) in opcodes.iter().enumerate() {
        if loop_stack == 0 {
            let instr = match op {
                OpCode::IncrementPointer => Some(Instr::IncrementPointer),
                OpCode::DecrementPointer => Some(Instr::DecrementPointer),
                OpCode::Increment => Some(Instr::Increment),
                OpCode::Decrement => Some(Instr::Decrement),
                OpCode::Write => Some(Instr::Write),
                OpCode::Read => Some(Instr::Read),
                OpCode::LoopBegin => {
                    loop_start = i;
                    loop_stack += 1;
                    None
                }
                OpCode::LoopEnd => panic!("loop ending at #{} has no beginning", i),
            };

            match instr {
                Some(instr) => program.push(instr),
                None => (),
            }
        } else {
            match op {
                OpCode::LoopBegin => {
                    loop_stack += 1;
                }
                OpCode::LoopEnd => {
                    loop_stack -= 1;
                    if loop_stack == 0 {
                        program.push(Instr::Loop(parse(&opcodes[loop_start + 1..i])));
                    }
                }
                _ => (),
            }
        }
    }

    if loop_stack != 0 {
        panic!(
            "loop that starts at #{} has no matching ending!",
            loop_start
        );
    }

    program
}

/// Executes a program that was previously parsed
fn run(
    instrs: &[Instr],
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
    tape: &mut [u8],
    data_pointer: &mut usize,
    counter: &mut usize,
) {
    for instr in instrs {
        *counter += 1;
        match instr {
            Instr::IncrementPointer => *data_pointer += 1,
            Instr::DecrementPointer => *data_pointer -= 1,
            Instr::Increment => tape[*data_pointer] += 1,
            Instr::Decrement => tape[*data_pointer] -= 1,
            Instr::Write => output
                .write_all(&tape[*data_pointer..*data_pointer + 1])
                .expect("failed to write output"),
            Instr::Read => {
                let mut x = [0u8; 1];
                input.read_exact(&mut x).expect("failed to read input");
                tape[*data_pointer] = x[0];
            }
            Instr::Loop(nested_instrs) => {
                while tape[*data_pointer] != 0 {
                    run(nested_instrs, input, output, tape, data_pointer, counter);
                }
            }
        }
    }
}

pub fn execute(
    source: &str,
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
) -> usize {
    let opcodes = lex(source);
    let program = parse(&opcodes);
    let mut tape: Vec<u8> = vec![0; 1024];
    let mut data_pointer = 512;
    let mut runtime = 0;
    run(
        &program,
        input,
        output,
        &mut tape,
        &mut data_pointer,
        &mut runtime,
    );
    runtime
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stark::SimulationMatrices;
    use algebra::fp_u64::BaseFelt;
    use algebra::PrimeFelt;

    const HELLO_WORLD_SOURCE: &str = "
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

    #[test]
    fn hello_world() {
        let mut output = Vec::new();

        execute(HELLO_WORLD_SOURCE, &mut std::io::empty(), &mut output);

        assert_eq!(output, "Hello World!\n".as_bytes());
    }
}
