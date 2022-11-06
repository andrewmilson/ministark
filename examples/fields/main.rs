use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::One;
use ark_ff_optimized::fp64::Fp;

mod fq3;
use fq3::Fq3;

fn main() {
    let mut my_elem = Fq3::from_base_prime_field_elems(&[Fp::one(), Fp::one(), Fp::one()]).unwrap();
    let mut extension_element = Fq3::one();
    let generator = Fp::TWO_ADIC_ROOT_OF_UNITY;

    println!("one: {}", generator);
    println!("generator: {}", generator);
    extension_element *= generator;
    println!("mul: {}", extension_element.inverse().unwrap());
    println!("mul: {}", extension_element.pow([(1 << 32) - 1]));
    my_elem *= generator;
    println!("test: {:?}", (my_elem * my_elem).inverse());
}
