#![feature(int_log)]

use algebra::fp_u64::BaseFelt;
use algebra::Felt;
use algebra::StarkFelt;
use algebra::Univariate;
use brainfuck::stark::compile;
use brainfuck::stark::SimulationMatrices;
use mini_stark::number_theory_transform::number_theory_transform;
use num_traits::One;
use num_traits::Zero;
use stark::protocol::StandardProofStream;
use stark::BrainFuckStark;
use stark::StarkParams;

const FIB_TO_100_SOURCE: &str = "
+++++++++++
>+>>>>++++++++++++++++++++++++++++++++++++++++++++
>++++++++++++++++++++++++++++++++<<<<<<[>[>>>>>>+>
+<<<<<<<-]>>>>>>>[<<<<<<<+>>>>>>>-]<[>++++++++++[-
<-[>>+>+<<<-]>>>[<<<+>>>-]+<[>[-]<[-]]>[<<[>>>+<<<
-]>>[-]]<<]>>>[>>+>+<<<-]>>>[<<<+>>>-]+<[>[-]<[-]]
>[<<+>>[-]]<<<<<<<]>>>>>[+++++++++++++++++++++++++
+++++++++++++++++++++++.[-]]++++++++++<[->-<]>++++
++++++++++++++++++++++++++++++++++++++++++++.[-]<<
<<<<<<<<<<[>>>+>+<<<<-]>>>>[<<<<+>>>>-]<-[>>.>.<<<
[-]]<<[>>+>+<<<-]>>>[<<<+>>>-]<<[<+>-]>[<+>-]<<<-]
";

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
    } = brainfuck::stark::simulate::<BaseFelt>(&program, &mut std::io::empty(), &mut output);

    let running_time = processor_matrix.len();
    println!("{running_time}");
    // let memory_length = memory_matrix.len();

    let mut proof_stream = StandardProofStream::<BaseFelt>::new();
    // let params = StarkParams::new(8, 128);
    let params = StarkParams::new(4, 0);
    let mut bfs = BrainFuckStark::<BaseFelt, BaseFelt>::new(params);
    let res = bfs.prove(
        processor_matrix,
        memory_matrix,
        instruction_matrix,
        input_matrix,
        output_matrix,
        &mut proof_stream,
    );

    println!("Size: {}", res.len());
}

#[test]
fn zerofier() {
    let n = 1usize << 4;
    let offset = BaseFelt::GENERATOR;
    let root = BaseFelt::get_root_of_unity(n.ilog2());

    // x - 1
    let poly = Univariate::new(vec![-BaseFelt::one(), BaseFelt::one()]);
    let mut coefficients = poly.scale(offset).coefficients;
    coefficients.resize(n, BaseFelt::zero());

    let eval1 = number_theory_transform(&coefficients);
    let eval2 = (0..n)
        .map(|i| offset * root.pow(&[i as u64]) - BaseFelt::one())
        .collect::<Vec<BaseFelt>>();

    assert_eq!(eval1, eval2);
}
