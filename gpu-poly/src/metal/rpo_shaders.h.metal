#ifndef rpo_shaders_h
#define rpo_shaders_h

#include <metal_stdlib>
#include "felt_u64.h.metal"
using namespace metal;

namespace p18446744069414584321 {

constant const unsigned STATE_WIDTH = 12;
constant const unsigned CAPACITY = 4;
constant const unsigned DIGEST_SIZE = 4;
constant const unsigned NUM_ROUNDS = 7;

// RPO's 12x12 row major MDS matrix in Montgomery domain
constant const Fp MDS[STATE_WIDTH * STATE_WIDTH] = {
    Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), 
    Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), 
    Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), 
    Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), 
    Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), 
    Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), 
    Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), 
    Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), 
    Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), Fp(111669149670), 
    Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), Fp(34359738360), 
    Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065), Fp(98784247785), 
    Fp(98784247785), Fp(34359738360), Fp(111669149670), Fp(55834574835), Fp(42949672950), Fp(38654705655), Fp(30064771065), Fp(25769803770), Fp(94489280490), Fp(90194313195), Fp(34359738360), Fp(30064771065)
};

constant const Fp ROUND_CONSTANTS_0[STATE_WIDTH * NUM_ROUNDS] = {
    Fp(6936159699454947676), Fp(6871277616928621393), Fp(4226339945476756083), Fp(2261225084505152444), Fp(16808067423291017741), Fp(12862191241011323277), Fp(345720808813194915), Fp(10126368034161173654), Fp(840649715788759894), Fp(18155600607269645987), Fp(16577339120870559289), Fp(13749826054300849029),
    Fp(16047969944113931191), Fp(10474334246235299199), Fp(15773847146013662260), Fp(14401231158322525155), Fp(6009395255763488383), Fp(2108579439821148946), Fp(13820200715803196660), Fp(15968614366574245570), Fp(7529997729792773654), Fp(9429194013557833999), Fp(11639903126146281421), Fp(15759666882357935738),
    Fp(14807658266593669785), Fp(17258259860767641342), Fp(9534132615398591413), Fp(358719342502509866), Fp(7123090532818864651), Fp(734193187930710962), Fp(14873184913735487023), Fp(17965359964069906568), Fp(12664837478844326631), Fp(15575491070113731145), Fp(7221479899469196675), Fp(7328957460733188967),
    Fp(15088355010936495340), Fp(16762963605345901631), Fp(15278161326153175940), Fp(6257793333052173411), Fp(8418953127708045776), Fp(6523475766574412380), Fp(15192936988185261803), Fp(1578086224854546096), Fp(10840553425559156784), Fp(7453417405109536362), Fp(5173069484734008228), Fp(3284492202065476384),
    Fp(1724586709636399686), Fp(17997633752581871175), Fp(1284825320737914582), Fp(960534381847281815), Fp(6708901808183456837), Fp(8975591106768797316), Fp(52515315389099119), Fp(10009391031874081397), Fp(3091228317422201238), Fp(1063858230459024983), Fp(3396548655473917480), Fp(15046057790353688034),
    Fp(4867464583127666756), Fp(13816959924674544309), Fp(13931201815459591565), Fp(11494116713280125381), Fp(16823081743980874023), Fp(6760771226809185048), Fp(5346741505458044699), Fp(15124596060558844029), Fp(5332565678905773189), Fp(17640389307200936126), Fp(14049814539797608740), Fp(8882709539093378074),
    Fp(10507930462458090835), Fp(10669463960502417047), Fp(16753662827442720769), Fp(12967456627495301601), Fp(2989815121821278695), Fp(5894674479204135685), Fp(14187454698288462352), Fp(14795723369628125345), Fp(17260571099239679821), Fp(16009836214833755168), Fp(2009092225887788829), Fp(10838446069154019765),
};

constant const Fp ROUND_CONSTANTS_1[STATE_WIDTH * NUM_ROUNDS] = {
    Fp(8939123259393952351), Fp(14708045228210488368), Fp(18125168669810517809), Fp(9309821433754818185), Fp(4714467145607136006), Fp(1302482025306688824), Fp(34829973686821040), Fp(5637233680011148778), Fp(227119480134509573), Fp(2530972937109017559), Fp(7210163798538732239), Fp(955913576003606833), 
    Fp(4449617297638325218), Fp(10843671682695268638), Fp(13198957499160452915), Fp(11541825028620451829), Fp(10963484480734735121), Fp(4752902142121643229), Fp(3015289210993491059), Fp(16344286514680205966), Fp(1811079964700766606), Fp(12735664961476037524), Fp(5775391330037813314), Fp(18223625362487900986), 
    Fp(7222477607687412281), Fp(4215615082079701144), Fp(6177508277476483691), Fp(3491362079220677263), Fp(10961785333913978630), Fp(1935408839283360916), Fp(13974192629927279950), Fp(18013556876298568088), Fp(7565676920589638093), Fp(9265825103386412558), Fp(8061587790235022972), Fp(6806849270604947860), 
    Fp(8066442548506952806), Fp(12791828131640457742), Fp(9268748809821748950), Fp(17496234860625277598), Fp(13583894547367420658), Fp(13920282495726802458), Fp(3933141341199584259), Fp(6658057712176150702), Fp(16812362035931029194), Fp(15160401867587809089), Fp(16411108749946146942), Fp(3390826434320009844), 
    Fp(18405475140095477472), Fp(13864039573264702148), Fp(496144052468360460), Fp(9791523668470936672), Fp(528582340156917005), Fp(15864481364569144493), Fp(682830611952089590), Fp(347158833826327515), Fp(13752775429919623417), Fp(10254722988306758482), Fp(8794150602427420596), Fp(2480344122229837853), 
    Fp(15462337562022968595), Fp(6729968753311049611), Fp(9250220857258211097), Fp(12031447985684644003), Fp(14538803180331344696), Fp(4055445230671851890), Fp(14764039661528567501), Fp(2047787218814287270), Fp(8977863094202715520), Fp(6560450968915612407), Fp(9976241128570886075), Fp(17877509887772213755), 
    Fp(3549624494907837709), Fp(4253629935471652443), Fp(2859199883984623807), Fp(1087607721547343649), Fp(7907517619951970198), Fp(11306402795121903516), Fp(10168009948206732524), Fp(9177440083248248246), Fp(13169036816957726187), Fp(12924186209140199217), Fp(9673006056831483321), Fp(747828276541750689)
};

[[ host_name("rpo_absorb_and_permute") ]] kernel void 
RpoAbsorbAndPermute(constant Fp *col0 [[ buffer(0) ]],
        constant Fp *col1 [[ buffer(1) ]],
        constant Fp *col2 [[ buffer(2) ]],
        constant Fp *col3 [[ buffer(3) ]],
        constant Fp *col4 [[ buffer(4) ]],
        constant Fp *col5 [[ buffer(5) ]],
        constant Fp *col6 [[ buffer(6) ]],
        constant Fp *col7 [[ buffer(7) ]],
        device Fp *states [[ buffer(8) ]],
        device Fp *digests [[ buffer(9) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) { 
    // absorb the input into the state
    // TODO: use 2d grid
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    // only store the first CAPACITY many items of the state in global memory
    unsigned global_state_offset = global_id * CAPACITY;
    
    // current state
    shared[local_state_offset + 0] = states[global_state_offset + 0];
    shared[local_state_offset + 1] = states[global_state_offset + 1];
    shared[local_state_offset + 2] = states[global_state_offset + 2];
    shared[local_state_offset + 3] = states[global_state_offset + 3];
    shared[local_state_offset + CAPACITY + 0] = col0[global_id];
    shared[local_state_offset + CAPACITY + 1] = col1[global_id];
    shared[local_state_offset + CAPACITY + 2] = col2[global_id];
    shared[local_state_offset + CAPACITY + 3] = col3[global_id];
    shared[local_state_offset + CAPACITY + 4] = col4[global_id];
    shared[local_state_offset + CAPACITY + 5] = col5[global_id];
    shared[local_state_offset + CAPACITY + 6] = col6[global_id];
    shared[local_state_offset + CAPACITY + 7] = col7[global_id];

#pragma unroll
    for (unsigned i = 0; i < NUM_ROUNDS; i++) {
#pragma unroll
        for (unsigned m = 0; m < STATE_WIDTH; m++) {
            Fp new_val = Fp(0);
#pragma unroll
            for (unsigned n = 0; n < STATE_WIDTH; n++) {
                Fp old = shared[local_state_offset + n];
                new_val = new_val + old * MDS[m * STATE_WIDTH + n];
            }
            new_val = new_val + ROUND_CONSTANTS_0[i * STATE_WIDTH + m];
            shared[local_state_offset + STATE_WIDTH + m] = new_val.pow7();
        }

#pragma unroll
        for (unsigned m = 0; m < STATE_WIDTH; m++) {
            Fp new_val = Fp(0);
#pragma unroll
            for (unsigned n = 0; n < STATE_WIDTH; n++) {
                Fp old = shared[local_state_offset + STATE_WIDTH + n];
                new_val = new_val + old * MDS[m * STATE_WIDTH + n];
            }
            new_val = new_val + ROUND_CONSTANTS_1[i * STATE_WIDTH + m];
            shared[local_state_offset + m] = new_val.pow10540996611094048183();
        }
    }

    states[global_state_offset + 0] = shared[local_state_offset + 0];
    states[global_state_offset + 1] = shared[local_state_offset + 1];
    states[global_state_offset + 2] = shared[local_state_offset + 2];
    states[global_state_offset + 3] = shared[local_state_offset + 3];

    digests[global_state_offset + 0] = shared[local_state_offset + CAPACITY + 0];
    digests[global_state_offset + 1] = shared[local_state_offset + CAPACITY + 1];
    digests[global_state_offset + 2] = shared[local_state_offset + CAPACITY + 2];
    digests[global_state_offset + 3] = shared[local_state_offset + CAPACITY + 3];




// // absorb the input into the state
//     // TODO: use 2d grid
//     unsigned local_state_offset = local_id * STATE_WIDTH * 2;
//     // only store the first CAPACITY many items of the state in global memory
//     unsigned global_state_offset = global_id * CAPACITY;
    
//     // current state
//     shared[local_state_offset + 0] = states[global_state_offset + 0];
//     shared[local_state_offset + 1] = states[global_state_offset + 1];
//     shared[local_state_offset + 2] = states[global_state_offset + 2];
//     shared[local_state_offset + 3] = states[global_state_offset + 3];
//     shared[local_state_offset + CAPACITY + 0] = col0[global_id];
//     shared[local_state_offset + CAPACITY + 1] = col1[global_id];
//     shared[local_state_offset + CAPACITY + 2] = col2[global_id];
//     shared[local_state_offset + CAPACITY + 3] = col3[global_id];
//     shared[local_state_offset + CAPACITY + 4] = col4[global_id];
//     shared[local_state_offset + CAPACITY + 5] = col5[global_id];
//     shared[local_state_offset + CAPACITY + 6] = col6[global_id];
//     shared[local_state_offset + CAPACITY + 7] = col7[global_id];

//     // next state
//     shared[local_state_offset + STATE_WIDTH + 0] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + 1] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + 2] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + 3] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 0] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 1] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 2] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 3] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 4] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 5] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 6] = Fp(0);
//     shared[local_state_offset + STATE_WIDTH + CAPACITY + 7] = Fp(0);

// #pragma unroll
//     for (unsigned i = 0; i < NUM_ROUNDS; i++) {
// #pragma unroll
//         for (unsigned m = 0; m < STATE_WIDTH; m++) {
//             Fp new_val = Fp(0);
// #pragma unroll
//             for (unsigned n = 0; n < STATE_WIDTH; n++) {
//                 Fp old = shared[local_state_offset + n];
//                 new_val = new_val + old * MDS[m * STATE_WIDTH + n];
//             }
//             shared[local_state_offset + STATE_WIDTH + m] = new_val;
//         }

// #pragma unroll
//         for (unsigned n = 0; n < STATE_WIDTH; n++) {
//             Fp tmp = shared[local_state_offset + STATE_WIDTH + n];
//             tmp = tmp + ROUND_CONSTANTS_0[i * STATE_WIDTH + n];
//             shared[local_state_offset + n] = tmp.pow7();
//             shared[local_state_offset + STATE_WIDTH + n] = Fp(0);
//         }

// #pragma unroll
//         for (unsigned m = 0; m < STATE_WIDTH; m++) {
//             Fp new_val = Fp(0);
// #pragma unroll
//             for (unsigned n = 0; n < STATE_WIDTH; n++) {
//                 Fp old = shared[local_state_offset + n];
//                 new_val = new_val + old * MDS[m * STATE_WIDTH + n];
//             }
//             shared[local_state_offset + STATE_WIDTH + m] = new_val;
//         }

// #pragma unroll
//         for (unsigned n = 0; n < STATE_WIDTH; n++) {
//             Fp tmp = shared[local_state_offset + STATE_WIDTH + n];
//             tmp = tmp + ROUND_CONSTANTS_1[i * STATE_WIDTH + n];
//             shared[local_state_offset + n] = tmp.pow10540996611094048183();
//             shared[local_state_offset + STATE_WIDTH + n] = Fp(0);
//         }
//     }

//     states[global_state_offset + 0] = shared[local_state_offset + 0];
//     states[global_state_offset + 1] = shared[local_state_offset + 1];
//     states[global_state_offset + 2] = shared[local_state_offset + 2];
//     states[global_state_offset + 3] = shared[local_state_offset + 3];

//     digests[global_state_offset + 0] = shared[local_state_offset + CAPACITY + 0];
//     digests[global_state_offset + 1] = shared[local_state_offset + CAPACITY + 1];
//     digests[global_state_offset + 2] = shared[local_state_offset + CAPACITY + 2];
//     digests[global_state_offset + 3] = shared[local_state_offset + CAPACITY + 3];




// FASTER VERSION BUT NOT ALWAYS CORRECT
// // absorb the input into the state
//     // TODO: use 2d grid
//     unsigned local_state_offset = local_id * STATE_WIDTH * 2;
//     // only store the first CAPACITY many items of the state in global memory
//     unsigned global_state_offset = global_id * CAPACITY;
    
//     // current state
//     shared[local_state_offset + 0] = states[global_state_offset + 0];
//     shared[local_state_offset + 1] = states[global_state_offset + 1];
//     shared[local_state_offset + 2] = states[global_state_offset + 2];
//     shared[local_state_offset + 3] = states[global_state_offset + 3];
//     shared[local_state_offset + CAPACITY + 0] = col0[global_id];
//     shared[local_state_offset + CAPACITY + 1] = col1[global_id];
//     shared[local_state_offset + CAPACITY + 2] = col2[global_id];
//     shared[local_state_offset + CAPACITY + 3] = col3[global_id];
//     shared[local_state_offset + CAPACITY + 4] = col4[global_id];
//     shared[local_state_offset + CAPACITY + 5] = col5[global_id];
//     shared[local_state_offset + CAPACITY + 6] = col6[global_id];
//     shared[local_state_offset + CAPACITY + 7] = col7[global_id];


//     // shared[local_state_offset + STATE_WIDTH + 0] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 1] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 2] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 3] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 4] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 5] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 6] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 7] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 8] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 9] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 10] = Fp(0);
//     // shared[local_state_offset + STATE_WIDTH + 11] = Fp(0);

//     Fp new_val = Fp(0);
// #pragma unroll
//     for (unsigned i = 0; i < NUM_ROUNDS; i++) {
// #pragma unroll
//         for (unsigned m = 0; m < STATE_WIDTH; m++) {
//             new_val = Fp(0);
// #pragma unroll
//             for (unsigned n = 0; n < STATE_WIDTH; n++) {
//                 Fp old = shared[local_state_offset + n];
//                 new_val = new_val + old * MDS[m * STATE_WIDTH + n];
//             }
//             new_val = new_val + ROUND_CONSTANTS_0[i * STATE_WIDTH + m];
//             shared[local_state_offset + STATE_WIDTH + m] = new_val.pow7();
//         }

// #pragma unroll
//         for (unsigned m = 0; m < STATE_WIDTH; m++) {
//             new_val = Fp(0);
// #pragma unroll
//             for (unsigned n = 0; n < STATE_WIDTH; n++) {
//                 Fp old = shared[local_state_offset + STATE_WIDTH + n];
//                 new_val = new_val + old * MDS[m * STATE_WIDTH + n];
//             }
//             new_val = new_val + ROUND_CONSTANTS_1[i * STATE_WIDTH + m];
//             shared[local_state_offset + m] = new_val.pow10540996611094048183();
//         }
//     }

//     states[global_state_offset + 0] = shared[local_state_offset + 0];
//     states[global_state_offset + 1] = shared[local_state_offset + 1];
//     states[global_state_offset + 2] = shared[local_state_offset + 2];
//     states[global_state_offset + 3] = shared[local_state_offset + 3];

//     digests[global_state_offset + 0] = shared[local_state_offset + CAPACITY + 0];
//     digests[global_state_offset + 1] = shared[local_state_offset + CAPACITY + 1];
//     digests[global_state_offset + 2] = shared[local_state_offset + CAPACITY + 2];
//     digests[global_state_offset + 3] = shared[local_state_offset + CAPACITY + 3];





// a b c
// d e f
// *
// g
// h
// l
// =
// a * g + b * h + c * l
// d * g + h * e + f * l


//     // absorb the input into the state
//     // TODO: use 2d grid
//     unsigned local_state_offset = (local_id / THREADS_PER_HASHER) * STATE_WIDTH * 2;
//     // only store the first CAPACITY many items of the state in global memory
//     unsigned global_state_offset = (global_id / THREADS_PER_HASHER) * CAPACITY;
    
//     if (local_id % THREADS_PER_HASHER == 0) {
//         // current state
//         shared[local_state_offset + 0] = states[global_state_offset + 0];
//         shared[local_state_offset + 1] = states[global_state_offset + 1];
//         shared[local_state_offset + 2] = states[global_state_offset + 2];
//         shared[local_state_offset + 3] = states[global_state_offset + 3];
//         shared[local_state_offset + CAPACITY + 0] = col0[global_id];
//         shared[local_state_offset + CAPACITY + 1] = col1[global_id];
//         shared[local_state_offset + CAPACITY + 2] = col2[global_id];
//         shared[local_state_offset + CAPACITY + 3] = col3[global_id];
//         shared[local_state_offset + CAPACITY + 4] = col4[global_id];
//         shared[local_state_offset + CAPACITY + 5] = col5[global_id];
//         shared[local_state_offset + CAPACITY + 6] = col6[global_id];
//         shared[local_state_offset + CAPACITY + 7] = col7[global_id];

//         // next state
//         shared[local_state_offset + STATE_WIDTH + 0] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + 1] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + 2] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + 3] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 0] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 1] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 2] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 3] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 4] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 5] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 6] = Fp(0);
//         shared[local_state_offset + STATE_WIDTH + CAPACITY + 7] = Fp(0);
//     }

//     simdgroup_barrier(mem_flags::mem_threadgroup);

//     unsigned hasher_id = local_id % THREADS_PER_HASHER;
//     if (hasher_id < 12) {
//         Fp tmp = Fp(0);
// #pragma unroll
//         for (unsigned i = 0; i < NUM_ROUNDS; i++) {
// #pragma unroll
//             for (unsigned n = 0; n < STATE_WIDTH; n++) {
//                 tmp = shared[local_state_offset + n];
//                 tmp = tmp * MDS[hasher_id * STATE_WIDTH + n];
//                 tmp = tmp + shared[local_state_offset + STATE_WIDTH + hasher_id];
//                 shared[local_state_offset + STATE_WIDTH + hasher_id] = tmp;
//             }

//             tmp = shared[local_state_offset + STATE_WIDTH + hasher_id];
//             tmp = tmp + ROUND_CONSTANTS_0[i * STATE_WIDTH + hasher_id];
//             shared[local_state_offset + hasher_id] = tmp.pow7();
//             shared[local_state_offset + STATE_WIDTH + hasher_id] = Fp(0);

//             simdgroup_barrier(mem_flags::mem_threadgroup);

// #pragma unroll
//             for (unsigned n = 0; n < STATE_WIDTH; n++) {
//                 tmp = shared[local_state_offset + n];
//                 tmp = tmp * MDS[hasher_id * STATE_WIDTH + n];
//                 tmp = tmp + shared[local_state_offset + STATE_WIDTH + hasher_id];
//                 shared[local_state_offset + STATE_WIDTH + hasher_id] = tmp;
//             }

//             tmp = shared[local_state_offset + STATE_WIDTH + hasher_id];
//             tmp = tmp + ROUND_CONSTANTS_1[i * STATE_WIDTH + hasher_id];
//             shared[local_state_offset + hasher_id] = tmp.pow10540996611094048183();
//             shared[local_state_offset + STATE_WIDTH + hasher_id] = Fp(0);

//             simdgroup_barrier(mem_flags::mem_threadgroup);
//         }
//     }

//     if (local_id % THREADS_PER_HASHER == 0) {
//         unsigned digest_offset = (global_id / THREADS_PER_HASHER) * DIGEST_SIZE;
//         digests[digest_offset + 0] = shared[local_state_offset + 0];
//         digests[digest_offset + 1] = shared[local_state_offset + 1];
//         digests[digest_offset + 2] = shared[local_state_offset + 2];
//         digests[digest_offset + 3] = shared[local_state_offset + 3];
//     }

// #pragma unroll
//     for (unsigned i = 0; i < NUM_ROUNDS; i++) {
// // #pragma unroll
// //         for (unsigned m = 0; m < STATE_WIDTH; m++) {
// // #pragma unroll
// //             for (unsigned n = 0; n < STATE_WIDTH; n++) {
// //                 Fp old = shared[local_state_offset + n];
// //                 Fp mds = MDS[m * STATE_WIDTH + n];
// //                 Fp curr = shared[local_state_offset + STATE_WIDTH + m];
// //                 shared[local_state_offset + STATE_WIDTH + m] = curr + old * mds;
// //             }
// //         }

// //         for (unsigned n = 0; n < STATE_WIDTH; n++) {
// //             Fp tmp = shared[local_state_offset + STATE_WIDTH + n];
// //             tmp = tmp + ROUND_CONSTANTS_0[i * STATE_WIDTH + n];
// //             shared[local_state_offset + n] = tmp.pow7();
// //             shared[local_state_offset + STATE_WIDTH + n] = Fp(0);
// //         }


// // #pragma unroll
// //         for (unsigned m = 0; m < STATE_WIDTH; m++) {
// // #pragma unroll
// //             for (unsigned n = 0; n < STATE_WIDTH; n++) {
// //                 Fp old = shared[local_state_offset + n];
// //                 Fp mds = MDS[m * STATE_WIDTH + n];
// //                 Fp curr = shared[local_state_offset + STATE_WIDTH + m];
// //                 shared[local_state_offset + STATE_WIDTH + m] = curr + old * mds;
// //             }
// //         }

//         for (unsigned n = 0; n < STATE_WIDTH / THREADS_PER_HASHER; n++) {
//             unsigned thread_n = n * THREADS_PER_HASHER + (i % THREADS_PER_HASHER);

//             Fp tmp = shared[local_state_offset + STATE_WIDTH + thread_n];
//             tmp = tmp + ROUND_CONSTANTS_1[i * STATE_WIDTH + thread_n];
//             shared[local_state_offset + thread_n] = tmp.pow10540996611094048183();
//             shared[local_state_offset + STATE_WIDTH + thread_n] = Fp(0);
//         }

//         // for (unsigned n = 0; n < STATE_WIDTH / ; n++) {
//         //     Fp tmp = shared[local_state_offset + STATE_WIDTH + n];
//         //     tmp = tmp + ROUND_CONSTANTS_1[i * STATE_WIDTH + n];
//         //     shared[local_state_offset + n] = tmp.pow10540996611094048183();
//         //     shared[local_state_offset + STATE_WIDTH + n] = Fp(0);
//         // }
//     }

    // if i % THREADS_PER_HASHER == 0 {
    //     unsigned digest_offset = global_id * DIGEST_SIZE;
    //     digests[digest_offset + 0] = shared[local_state_offset + 0];
    //     digests[digest_offset + 1] = shared[local_state_offset + 1];
    //     digests[digest_offset + 2] = shared[local_state_offset + 2];
    //     digests[digest_offset + 3] = shared[local_state_offset + 3];
    // }


//     Fp current_state[STATE_WIDTH] = {
//         shared[local_state_offset + 0],
//         shared[local_state_offset + 1],
//         shared[local_state_offset + 2],
//         shared[local_state_offset + 3],
//         shared[local_state_offset + 4],
//         shared[local_state_offset + 5],
//         shared[local_state_offset + 6],
//         shared[local_state_offset + 7],
//         shared[local_state_offset + 8],
//         shared[local_state_offset + 9],
//         shared[local_state_offset + 10],
//         shared[local_state_offset + 11],
//     };

//     Fp next_state[STATE_WIDTH] = {
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//         Fp(0),
//     };

// // #pragma unroll
//     for (unsigned i = 0; i < 1; i++) {
//         // MUL_MDS(next_state, current_state);

//         for (unsigned num_am = 0; i < 12; num_am++) {
//             for (unsigned num_bn = 0; i < 1; num_bn++) {

//             }
//         }

//         for i in 0..AM {
//             for j in 0..BN {
//                 // println!("res[{i}] = ");
//                 for k in 0..AN_BM {
//                     // println!(" + MDS[{i} * STATE_WIDTH + {k}] * state[{k}]");
//                     res[i][j] += a[i][k] * b[k][j];
//                 }
//                 // println!();
//             }
//         }
//     }

    //     for (unsigned m = 0; m < STATE_WIDTH; m++) {
    //         for (unsigned n = 0; n < STATE_WIDTH; n++) {
    //             Fp res = state[n] * MDS[m * STATE_WIDTH + n];
    //             shared[local_state_offset + n] = state[n] + res;
    //         }
    //         // state[m]
    //     }

    //     // for (unsigned n = 0; n < STATE_WIDTH; n++) {
    //     //     shared[local_state_offset + n] = state[n];
    //     // }

    //     // threadgroup_barrier(mem_flags::mem_threadgroup);

    //     // MUL_MDS(state, (&shared[local_state_offset]));
    //     // for (unsigned j = 0; j < STATE_WIDTH; j++) {
    //     //     state[j] = state[j] + ROUND_CONSTANTS_0[i * STATE_WIDTH + j];
    //     // }
    //     // for (unsigned j = 0; j < STATE_WIDTH; j++) {
    //     //     state[j] = state[j].pow7();
    //     // }
    //     // UPDATE_STATE((&shared[local_state_offset]), state)


    //     // MUL_MDS(state, (&shared[local_state_offset]));
    //     // for (unsigned j = 0; j < STATE_WIDTH; j++) {
    //     //     state[j] = state[j] + ROUND_CONSTANTS_1[i * STATE_WIDTH + j];
    //     // }
    //     // for (unsigned j = 0; j < STATE_WIDTH; j++) {
    //     //     state[j] = state[j].pow10540996611094048183();
    //     // }
    //     // UPDATE_STATE((&shared[local_state_offset]), state)
    // }

        // MUL_MDS(state, (&shared[local_state_offset]));
        // ADD_ROUND_CONSTANTS(state, ROUND_CONSTANTS_0, i);
        // state[5] = Fp(30064771065);
        // Fp tmp = state[5].pow7().pow7();
        // state[5] = tmp;
        // SBOX(state);
        // UPDATE_STATE((&shared[local_state_offset]), state)

        // MUL_MDS(state, (&shared[local_state_offset]));
        // ADD_ROUND_CONSTANTS(state, ROUND_CONSTANTS_1, i);
        // SBOX_INV(state);
        // if (i != NUM_ROUNDS) {
        //     UPDATE_STATE((&shared[local_state_offset]), state)
        // }
    // }

    // MUL_MDS(state, (&shared[local_state_offset]));
    // MUL_MDS(state, ((threadgroup Fp*)(shared + local_state_offset)));
    // ADD_ROUND_CONSTANTS(state, ROUND_CONSTANTS_0, 0);

    // ADD_ROUND_CONSTANTS(state, ROUND_CONSTANTS_0, 0);

    // state[0] = ROUND_CONSTANTS_0[0 * STATE_WIDTH + 0];
    // state[1] = ROUND_CONSTANTS_0[0 * STATE_WIDTH + 0];
    // state[2] = ROUND_CONSTANTS_0[0 * STATE_WIDTH + 0];
    // state[3] = ROUND_CONSTANTS_0[0 * STATE_WIDTH + 0];
    // state[4] = ROUND_CONSTANTS_0[0 * STATE_WIDTH + 0];
    // state[5] = ROUND_CONSTANTS_0[0 * STATE_WIDTH + 0];



// #pragma unroll
//     for (unsigned i = 0; i < NUM_ROUNDS; i++) {
        // MUL_MDS(state, (&shared[local_state_offset]));
        // ADD_ROUND_CONSTANTS(state, ROUND_CONSTANTS_0, i);
        // SBOX(state);
        // UPDATE_STATE((&shared[local_state_offset]), state)

        // MUL_MDS(state, (&shared[local_state_offset]));
        // ADD_ROUND_CONSTANTS(state, ROUND_CONSTANTS_1, i);
        // SBOX_INV(state);
        // if (i != ) {
            // UPDATE_STATE((&shared[local_state_offset]), state)
        // }
    // }

    // // store digest
    // unsigned digest_offset = global_id * DIGEST_SIZE;
    // digests[digest_offset + 0] = shared[local_state_offset + CAPACITY + 0];
    // digests[digest_offset + 1] = shared[local_state_offset + CAPACITY + 1];
    // digests[digest_offset + 2] = shared[local_state_offset + CAPACITY + 2];
    // digests[digest_offset + 3] = shared[local_state_offset + CAPACITY + 3];
}

}

#endif /* rpo_shaders_h */
