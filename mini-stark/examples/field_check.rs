use mini_stark::{prime_field_u128::BaseElement, StarkElement};

fn main() {
    // let generator = BaseElement::GENERATOR;
    let root2 = BaseElement::get_root_of_unity(1);
    let root4 = BaseElement::get_root_of_unity(2);

    println!("{}", BaseElement::from(0b11011u64));

    println!(
        "{}",
        BaseElement::from(101u64 * 16) / BaseElement::from(16u64)
    )
}
