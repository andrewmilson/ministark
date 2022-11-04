use ark_ff::BigInt;
use ark_ff::FftField;
use ark_ff::Field;
use ark_ff::Fp3;
use ark_ff::Fp3Config;
use ark_ff::MontFp;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_ff_optimized::fp64::Fp;
use ark_ff_optimized::fp64::FpParams;
use ark_poly::EvaluationDomain;
use ark_poly::Radix2EvaluationDomain;
use std::marker::PhantomData;
use std::ops::MulAssign;

#[derive(ark_ff::MontConfig)]
#[modulus = "475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137"]
#[generator = "10"]
pub struct FrConfig;
pub type Fq = ark_ff::Fp320<ark_ff::MontBackend<FrConfig, 5>>;

pub type Fq3 = Fp3<Fq3Config>;

pub struct Fq3Config;

impl Fp3Config for Fq3Config {
    type Fp = Fq;

    const NONRESIDUE: Fq = MontFp!("5");

    const TWO_ADICITY: u32 = 34;

    #[rustfmt::skip]
    const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = &[
        0x69232b75663933bd,
        0xca650efcfc00ee0,
        0x77ca3963fe36f720,
        0xe4cb46632f9bcf7e,
        0xef510453f08f9f30,
        0x9dd5b8fc72f02d83,
        0x7f8d017ed86608ab,
        0xeb2219b3697c97a4,
        0xc8663846ab96996f,
        0x833cd532053eac7d,
        0x1d5b73dfb20bd3cc,
        0x6f5f6da606b59873,
        0x62e990f43dfc42d6,
        0x6878f58,
    ];

    const QUADRATIC_NONRESIDUE_TO_T: Fq3 = Fq3::new(
        MontFp!("154361449678783505076984156275977937654331103361174469632346230549735979552469642799720052"),
        Fq::ZERO,
        Fq::ZERO,
    );

    const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[
        Fq::ONE,
        MontFp!("471738898967521029133040851318449165997304108729558973770077319830005517129946578866686956"),
        MontFp!("4183387201740296620308398334599285547820769823264541783190415909159130177461911693276180"),
    ];

    const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[
        Self::FROBENIUS_COEFF_FP3_C1[0],
        Self::FROBENIUS_COEFF_FP3_C1[2],
        Self::FROBENIUS_COEFF_FP3_C1[1],
    ];
}

impl MulAssign<Fq> for Fq3 {
    fn mul_assign(&mut self, rhs: Fq) {
        todo!()
    }
}

// pub struct Fq3Config;

// impl Fp3Config for Fq3Config {
//     type Fp = Fp;
//     const NONRESIDUE: Fp = /* =2 */ ark_ff::Fp(BigInt([8589934590]),
// PhantomData);     const TWO_ADICITY: u32 = Fp::TWO_ADICITY;

//     #[rustfmt::skip]
//     const TRACE_MINUS_ONE_DIV_TWO: &'static [u64] = FpParams::TR;

//     /// (11^T, 0, 0)
//     #[rustfmt::skip]
//     const QUADRATIC_NONRESIDUE_TO_T: Fq3 = Fq3::new(
//         MontFp!("22168644070733283197994897338612733221095941481265408161807376791727499343083607817089033595478370212662133368413166734396127674284827734481031659015434501966360165723728649019457855887066657739809176476252080335185730833468062"),
//         FQ_ZERO,
//         FQ_ZERO,
//     );

//     // Coefficients for the Frobenius automorphism.
//     // c1[0] = 1,
//     // c1[1] =
// 24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132
//     // c1[2] =
// 17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868,
//     #[rustfmt::skip]
//     const FROBENIUS_COEFF_FP3_C1: &'static [Fq] = &[
//         FQ_ONE,
//         MontFp!("24129022407817241407134263419936114379815707076943508280977368156625538709102831814843582780138963119807143081677569721953561801075623741378629346409604471234573396989178424163772589090105392407118197799904755622897541183052132"),
//         MontFp!("17769468560101711995209951371304522748355002843010440790806134764399814103468274958215310983651375801610927890210888755369611256415970113691066895445191924931148019336171640277697829047741006062493737919155152541323243293107868"),
//     ];

//     // c2 = {c1[0], c1[2], c1[1]}
//     #[rustfmt::skip]
//     const FROBENIUS_COEFF_FP3_C2: &'static [Fq] = &[
//         FQ_ONE,
//         Self::FROBENIUS_COEFF_FP3_C1[2],
//         Self::FROBENIUS_COEFF_FP3_C1[1],
//     ];
// }

fn main() {
    let a = Fq3::new(Fq::one(), Fq::one(), Fq::one());
    let b = Fq::one();

    let domain = Radix2EvaluationDomain::<Fp3>::new(8).unwrap();
    let evals = domain.fft(&[
        Fp3::one(),
        Fp3::one(),
        Fp3::one(),
        Fp3::one(),
        Fp3::one(),
        Fp3::one(),
        Fp3::one(),
        Fp3::one(),
    ]);

    println!("Test: {:?}", evals);
}
