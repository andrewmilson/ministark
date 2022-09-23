use ark_ff::FftField;
use ark_ff::Field;

pub fn compute_io_terminal<F>(symbols: &[usize], challenge: F) -> F
where
    F: Field,
    F::BasePrimeField: FftField,
{
    let mut acc = F::zero();
    for &s in symbols {
        acc = challenge * acc + F::from(s as u64);
    }
    acc
}

pub fn compute_program_terminal<F>(program: &[usize], challenges: &[F]) -> F
where
    F: Field,
    F::BasePrimeField: FftField,
{
    let mut challenges_iter = challenges.iter().copied();
    let a = challenges_iter.next().unwrap();
    let b = challenges_iter.next().unwrap();
    let c = challenges_iter.next().unwrap();
    let _d = challenges_iter.next().unwrap();
    let _e = challenges_iter.next().unwrap();
    let _f = challenges_iter.next().unwrap();
    let _alpha = challenges_iter.next().unwrap();
    let _beta = challenges_iter.next().unwrap();
    let _gamma = challenges_iter.next().unwrap();
    let _delta = challenges_iter.next().unwrap();
    let eta = challenges_iter.next().unwrap();

    let mut running_sum = F::zero();
    // Placeholder set to some invalid value
    let mut prev_ip = -F::one();
    let mut padded_program = program
        .iter()
        .map(|p| F::from(*p as u64))
        .collect::<Vec<F>>();
    padded_program.push(F::zero());

    for i in 0..padded_program.len() - 1 {
        let ip = F::from(i as u64);
        let curr_instr = padded_program[i];
        let next_instr = padded_program[i + 1];
        // TODO: remove this if. Is there any point to it?
        if prev_ip != ip {
            running_sum = running_sum * eta + a * ip + b * curr_instr + c * next_instr;
        }
        prev_ip = ip;
    }

    let index = padded_program.len() - 1;
    let address = F::from(index as u64);
    let curr_instr = padded_program[index];
    let next_instr = F::zero();
    running_sum = running_sum * eta + a * address + b * curr_instr + c * next_instr;

    running_sum
}
