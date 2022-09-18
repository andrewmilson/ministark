use fast_poly::allocator::PageAlignedAllocator;
use legacy_algebra::fp_u128::BaseFelt;
use legacy_algebra::Felt;
use mini_stark::polynomial::MultivariatePolynomial;
use mini_stark::polynomial::Polynomial;
use num_traits::Zero;
use std::iter::once;

pub struct RescuePrime {
    // Determines the state width of the hash function
    // Phrased differently, in
    // the evaluation of the function, the state is fully determined by m > 1 field
    // elements.
    pub m: usize,
    capacity: usize,
    // Number of rounds
    pub N: usize,
    alpha: u128,
    alpha_inverse: u128,
    mds: Vec<Vec<BaseFelt>>,
    mds_inverse: Vec<Vec<BaseFelt>>,
    round_constants: Vec<BaseFelt>,
}

impl RescuePrime {
    pub fn new() -> RescuePrime {
        RescuePrime {
            m: 2,
            capacity: 1,
            N: 27,
            alpha: 3,
            alpha_inverse: 180331931428153586757283157844700080811,
            mds: vec![
                vec![
                    BaseFelt::new(270497897142230380135924736767050121214),
                    BaseFelt::new(4),
                ],
                vec![
                    BaseFelt::new(270497897142230380135924736767050121205),
                    BaseFelt::new(13),
                ],
            ],
            mds_inverse: vec![
                vec![
                    BaseFelt::new(210387253332845851216830350818816760948),
                    BaseFelt::new(60110643809384528919094385948233360270),
                ],
                vec![
                    BaseFelt::new(90165965714076793378641578922350040407),
                    BaseFelt::new(180331931428153586757283157844700080811),
                ],
            ],
            round_constants: vec![
                174420698556543096520990950387834928928,
                109797589356993153279775383318666383471,
                228209559001143551442223248324541026000,
                268065703411175077628483247596226793933,
                250145786294793103303712876509736552288,
                154077925986488943960463842753819802236,
                204351119916823989032262966063401835731,
                57645879694647124999765652767459586992,
                102595110702094480597072290517349480965,
                8547439040206095323896524760274454544,
                50572190394727023982626065566525285390,
                87212354645973284136664042673979287772,
                64194686442324278631544434661927384193,
                23568247650578792137833165499572533289,
                264007385962234849237916966106429729444,
                227358300354534643391164539784212796168,
                179708233992972292788270914486717436725,
                102544935062767739638603684272741145148,
                65916940568893052493361867756647855734,
                144640159807528060664543800548526463356,
                58854991566939066418297427463486407598,
                144030533171309201969715569323510469388,
                264508722432906572066373216583268225708,
                22822825100935314666408731317941213728,
                33847779135505989201180138242500409760,
                146019284593100673590036640208621384175,
                51518045467620803302456472369449375741,
                73980612169525564135758195254813968438,
                31385101081646507577789564023348734881,
                270440021758749482599657914695597186347,
                185230877992845332344172234234093900282,
                210581925261995303483700331833844461519,
                233206235520000865382510460029939548462,
                178264060478215643105832556466392228683,
                69838834175855952450551936238929375468,
                75130152423898813192534713014890860884,
                59548275327570508231574439445023390415,
                43940979610564284967906719248029560342,
                95698099945510403318638730212513975543,
                77477281413246683919638580088082585351,
                206782304337497407273753387483545866988,
                141354674678885463410629926929791411677,
                19199940390616847185791261689448703536,
                177613618019817222931832611307175416361,
                267907751104005095811361156810067173120,
                33296937002574626161968730356414562829,
                63869971087730263431297345514089710163,
                200481282361858638356211874793723910968,
                69328322389827264175963301685224506573,
                239701591437699235962505536113880102063,
                17960711445525398132996203513667829940,
                219475635972825920849300179026969104558,
                230038611061931950901316413728344422823,
                149446814906994196814403811767389273580,
                25535582028106779796087284957910475912,
                93289417880348777872263904150910422367,
                4779480286211196984451238384230810357,
                208762241641328369347598009494500117007,
                34228805619823025763071411313049761059,
                158261639460060679368122984607245246072,
                65048656051037025727800046057154042857,
                134082885477766198947293095565706395050,
                23967684755547703714152865513907888630,
                8509910504689758897218307536423349149,
                232305018091414643115319608123377855094,
                170072389454430682177687789261779760420,
                62135161769871915508973643543011377095,
                15206455074148527786017895403501783555,
                201789266626211748844060539344508876901,
                179184798347291033565902633932801007181,
                9615415305648972863990712807943643216,
                95833504353120759807903032286346974132,
                181975981662825791627439958531194157276,
                267590267548392311337348990085222348350,
                49899900194200760923895805362651210299,
                89154519171560176870922732825690870368,
                265649728290587561988835145059696796797,
                140583850659111280842212115981043548773,
                266613908274746297875734026718148328473,
                236645120614796645424209995934912005038,
                265994065390091692951198742962775551587,
                59082836245981276360468435361137847418,
                26520064393601763202002257967586372271,
                108781692876845940775123575518154991932,
                138658034947980464912436420092172339656,
                45127926643030464660360100330441456786,
                210648707238405606524318597107528368459,
                42375307814689058540930810881506327698,
                237653383836912953043082350232373669114,
                236638771475482562810484106048928039069,
                168366677297979943348866069441526047857,
                195301262267610361172900534545341678525,
                2123819604855435621395010720102555908,
                96986567016099155020743003059932893278,
                248057324456138589201107100302767574618,
                198550227406618432920989444844179399959,
                177812676254201468976352471992022853250,
                211374136170376198628213577084029234846,
                105785712445518775732830634260671010540,
                122179368175793934687780753063673096166,
                126848216361173160497844444214866193172,
                22264167580742653700039698161547403113,
                234275908658634858929918842923795514466,
                189409811294589697028796856023159619258,
                75017033107075630953974011872571911999,
                144945344860351075586575129489570116296,
                261991152616933455169437121254310265934,
                18450316039330448878816627264054416127,
            ]
            .into_iter()
            .map(BaseFelt::new)
            .collect(),
        }
    }

    fn rate(&self) -> usize {
        self.m - self.capacity
    }

    pub fn hash(&self, input_element: BaseFelt) -> BaseFelt {
        // absorb
        let mut state = once(input_element)
            .chain(once(BaseFelt::zero()).take(self.m - 1))
            .collect::<Vec<BaseFelt>>();

        let alpha_high = (self.alpha >> 64) as u64;
        let alpha_low = self.alpha as u64;

        let alpha_inv_high = (self.alpha_inverse >> 64) as u64;
        let alpha_inv_low = self.alpha_inverse as u64;

        // permutation
        for r in 0..self.N {
            // forward half round
            // S-box
            for element in state.iter_mut() {
                *element = (*element).pow(&[alpha_low, alpha_high]);
            }
            // matrix
            let mut temp = vec![BaseFelt::zero(); self.m];
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] += self.mds[i][j] * state[j];
                }
            }
            // constants
            state = (0..self.m)
                .map(|i| temp[i] + self.round_constants[2 * r * self.m + i])
                .collect::<Vec<BaseFelt>>();

            // backwards half round
            // S-box
            for element in state.iter_mut() {
                *element = (*element).pow(&[alpha_inv_low, alpha_inv_high]);
            }
            // matrix
            let mut temp = vec![BaseFelt::zero(); self.m];
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] += self.mds[i][j] * state[j];
                }
            }
            // constants
            state = (0..self.m)
                .map(|i| temp[i] + self.round_constants[2 * r * self.m + self.m + i])
                .collect::<Vec<BaseFelt>>();
        }

        state[0]
    }

    fn round_constants_polynomials(
        &self,
        omicron: BaseFelt,
    ) -> (
        Vec<MultivariatePolynomial<BaseFelt>>,
        Vec<MultivariatePolynomial<BaseFelt>>,
    ) {
        let domain = (0..self.N)
            .map(|round| omicron.pow(&[round as u64]))
            .collect::<Vec<BaseFelt>>();

        let mut first_step_constants = vec![];
        for i in 0..self.m {
            let values = (0..self.N)
                .map(|round| self.round_constants[2 * round * self.m + i])
                .collect::<Vec<BaseFelt>>();
            let univariate = Polynomial::interpolate(&domain, &values);
            let multivariate = MultivariatePolynomial::lift(univariate, 0);
            first_step_constants.push(multivariate);
        }

        let mut second_step_constants = vec![];
        for i in 0..self.m {
            let values = (0..self.N)
                .map(|round| self.round_constants[2 * round * self.m + self.m + i])
                .collect::<Vec<BaseFelt>>();
            let univariate = Polynomial::interpolate(&domain, &values);
            let multivariate = MultivariatePolynomial::lift(univariate, 0);
            second_step_constants.push(multivariate);
        }

        (first_step_constants, second_step_constants)
    }

    pub fn transition_constraints(
        &self,
        omicron: BaseFelt,
    ) -> Vec<MultivariatePolynomial<BaseFelt>> {
        // get polynomials that interpolate through the round constants
        let (first_step_constants, second_step_constants) =
            self.round_constants_polynomials(omicron);

        // arithmetize one round of Rescue-Prime
        let variables = MultivariatePolynomial::variables(1 + 2 * self.m);
        let previous_state = &variables[1..(1 + self.m)];
        let next_state = &variables[(1 + self.m)..(1 + 2 * self.m)];
        let mut air = vec![];

        for i in 0..self.m {
            // compute the left hand side symbolically
            // LHS = sum(MultivariatePolynomial::constant(self.MDS[i][k]) *
            // (previous_state[k] ^ self.alpha) for k in 0..self.m)
            let mut lhs = MultivariatePolynomial::constant(BaseFelt::zero());
            for k in 0..self.m {
                lhs = lhs
                    + MultivariatePolynomial::constant(self.mds[i][k])
                        * (previous_state[k].clone() ^ self.alpha);
            }
            lhs = lhs + first_step_constants[i].clone();

            // compute the right hand side symbolically
            let mut rhs = MultivariatePolynomial::constant(BaseFelt::zero());
            for k in 0..self.m {
                rhs = rhs
                    + MultivariatePolynomial::constant(self.mds_inverse[i][k])
                        * (next_state[k].clone() - second_step_constants[k].clone());
            }
            rhs = rhs ^ self.alpha;

            // equate left and right hand sides
            air.push(lhs - rhs);
        }

        air
    }

    // boundary constraints (cycle, register, value)
    pub fn boundary_constraints(&self, output_element: BaseFelt) -> Vec<(usize, usize, BaseFelt)> {
        vec![
            // at the start, capacity is zero
            (0, 1, BaseFelt::zero()),
            // at the end, rate part is the given output element
            (self.N, 0, output_element),
        ]
    }

    pub fn trace(&self, input_element: BaseFelt) -> Vec<Vec<BaseFelt, PageAlignedAllocator>> {
        // absorb
        let mut state = once(input_element)
            .chain(once(BaseFelt::zero()).take(self.m - 1))
            .collect::<Vec<BaseFelt>>();

        // explicit copy and record state into trace

        let mut trace = vec![];
        let mut page_aligned_state = Vec::with_capacity_in(state.len(), PageAlignedAllocator);
        for &value in state.iter() {
            page_aligned_state.push(value)
        }
        trace.push(page_aligned_state);

        let alpha_high = (self.alpha >> 64) as u64;
        let alpha_low = self.alpha as u64;

        let alpha_inv_high = (self.alpha_inverse >> 64) as u64;
        let alpha_inv_low = self.alpha_inverse as u64;

        // permutation
        for round in 0..self.N {
            // forward half round
            // S-box
            for element in state.iter_mut() {
                *element = (*element).pow(&[alpha_low, alpha_high]);
            }
            // matrix
            let mut temp = vec![BaseFelt::zero(); self.m];
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] += self.mds[i][j] * state[j];
                }
            }
            // constants
            state = (0..self.m)
                .map(|i| temp[i] + self.round_constants[2 * round * self.m + i])
                .collect::<Vec<BaseFelt>>();

            // backwards half round
            // S-box
            for element in state.iter_mut() {
                *element = (*element).pow(&[alpha_inv_low, alpha_inv_high]);
            }
            // matrix
            let mut temp = vec![BaseFelt::zero(); self.m];
            for i in 0..self.m {
                for j in 0..self.m {
                    temp[i] += self.mds[i][j] * state[j];
                }
            }
            // constants
            state = (0..self.m)
                .map(|i| temp[i] + self.round_constants[2 * round * self.m + self.m + i])
                .collect::<Vec<BaseFelt>>();

            let mut page_aligned_state = Vec::with_capacity_in(state.len(), PageAlignedAllocator);
            // page_aligned_state.clone_from_slice(&state[..]);
            for &value in state.iter() {
                page_aligned_state.push(value)
            }
            trace.push(page_aligned_state);
        }

        trace
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::StarkElement;
    use num_traits::One;
    use rand::Rng;

    #[test]
    fn test_rescue_prime() {
        let rp = RescuePrime::new();

        // test vectors
        assert_eq!(
            rp.hash(BaseFelt::one()),
            BaseFelt::new(244180265933090377212304188905974087294),
            "rescue prime test vector 1 failed"
        );
        assert_eq!(
            rp.hash(BaseFelt::new(57322816861100832358702415967512842988)),
            BaseFelt::new(89633745865384635541695204788332415101),
            "rescue prime test vector 2 failed"
        );

        // test trace boundaries
        let a = BaseFelt::new(57322816861100832358702415967512842988);
        let b = BaseFelt::new(89633745865384635541695204788332415101);
        let trace = rp.trace(a);
        assert_eq!(
            trace[0][0], a,
            "rescue prime trace does not satisfy boundary conditions"
        );
        assert_eq!(
            trace[trace.len() - 1][0],
            b,
            "rescue prime trace does not satisfy boundary conditions"
        );
    }

    #[test]
    fn test_trace() {
        let rp = RescuePrime::new();

        let input_element = BaseFelt::new(57322816861100832358702415967512842988);
        let b = BaseFelt::new(89633745865384635541695204788332415101);
        let output_element = rp.hash(input_element);
        assert_eq!(b, output_element, "output elements do not match");

        // get trace
        let mut trace = rp.trace(input_element);

        // test boundary constraints
        for (cycle, element, value) in rp.boundary_constraints(output_element) {
            assert_eq!(trace[cycle][element], value, "rescue prime boundary condition error: trace element {} at cycle {} has value {} but should have value {}", element, cycle, trace[cycle][element], value);
        }

        // test transition constraints
        let omicron = BaseFelt::get_root_of_unity(BaseFelt::TWO_ADICITY);
        let transition_constraints = rp.transition_constraints(omicron);
        for o in 0..(trace.len() - 1) {
            for air_poly in transition_constraints.iter() {
                let previous_state = [trace[o][0], trace[o][1]];
                let next_state = [trace[o + 1][0], trace[o + 1][1]];
                let point = once(omicron.pow(o as u128))
                    .chain(previous_state)
                    .chain(next_state)
                    .collect::<Vec<BaseFelt>>();
                assert_eq!(
                    air_poly.evaluate(&point),
                    BaseFelt::zero(),
                    "air polynomial does not evaluate to zero"
                );
            }
        }

        println!("valid Rescue-Prime trace passes tests, testing invalid traces ...");

        let mut rng = rand::thread_rng();

        // insert errors into trace, to make sure errors get noticed
        for k in 0..10 {
            println!("trial {k}...");
            // sample error location and value randomly
            let mut register_index = rng.gen::<usize>() % rp.m;
            let mut cycle_index = rng.gen::<usize>() % (rp.N + 1);
            let mut value_ = BaseFelt::from(rng.gen::<u64>());

            if value_ == BaseFelt::zero() {
                continue;
            }

            // reproduce deterministic error
            if k == 0 {
                register_index = 1;
                cycle_index = 22;
                value_ = BaseFelt::new(17274817952119230544216945715808633996);
            }

            // perturb
            trace[cycle_index][register_index] += value_;

            let mut error_got_noticed = false;

            // test boundary constraints
            for (cycle, element, value) in rp.boundary_constraints(output_element) {
                if trace[cycle][element] != value {
                    error_got_noticed = true;
                    break;
                }
            }

            // test transition constraints
            for o in 0..(trace.len() - 1) {
                for air_poly in rp.transition_constraints(omicron) {
                    let previous_state = [trace[o][0], trace[o][1]];
                    let next_state = [trace[o + 1][0], trace[o + 1][1]];
                    let point = once(omicron.pow(o as u128))
                        .chain(previous_state)
                        .chain(next_state)
                        .collect::<Vec<BaseFelt>>();
                    if !air_poly.evaluate(&point).is_zero() {
                        error_got_noticed = true;
                        break;
                    }
                }
            }

            // if error was not noticed, panic
            if !error_got_noticed {
                println!("error was not noticed.");
                println!("register index: {register_index}");
                println!("cycle index: {cycle_index}");
                println!("value_: {}", value_);
                assert!(false, "error was not noticed");
            }

            // fix the error
            trace[cycle_index][register_index] -= value_;
        }
    }
}
