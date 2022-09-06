use brainfuck::stark::compile;
use stark::BrainFuckStark;

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
    let program = compile(HELLO_WORLD_SOURCE);
    let mut output = Vec::new();
    let SimulationMatrices {
        processor: processor_matrix,
        instruction: instruction_matrix,
        input: input_matrix,
        output: output_matrix,
        memory: memory_matrix,
    } = stark::simulate::<BaseFelt>(&program, &mut std::io::empty(), &mut output);

    // let running_time = processor_matrix.len();
    // let memory_length = memory_matrix.len();

    let bfs = BrainFuckStark::new(params);

    assert_eq!(running_time, 10);
}
