use std::time::Instant;

mod rescue_prime;

mod signature_scheme;
use signature_scheme::*;

use mini_stark::prime_field_u128;

fn main() {
    let document = String::from("This is a test");

    // TODO: add some love and care to this example: comments, new secret key and public key
    let signature_scheme = StarkSignatureScheme::new();

    let (secret_key, public_key) = (
        prime_field_u128::BaseElement::new(12765326281186373138),
        prime_field_u128::BaseElement::new(6214705293158922468180514843788232123),
    );

    println!("Generating signature...");
    let now = Instant::now();
    let signature = signature_scheme.sign(secret_key, document.clone());
    println!("Signature generated in {:.2?}", now.elapsed());
    println!("Signature size: {}KB\n\n", signature.len() / 1000);

    println!("Verifying signature...");
    let now = Instant::now();
    let is_valid = signature_scheme.verify(public_key, document, &signature);
    println!("Verified validity of signature in {:.2?}", now.elapsed());
    println!(
        "Signature is valid? {}",
        if is_valid.is_ok() { "yes" } else { "no" }
    );
}
