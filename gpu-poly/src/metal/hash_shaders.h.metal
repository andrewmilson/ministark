#ifndef hash_shaders_h
#define hash_shaders_h

#include <metal_stdlib>
#include "felt_u64.h.metal"
using namespace metal;

namespace p18446744069414584321 {

// Stores the first CAPACITY many items of the state.
// This is the only state that needs to be persisted between
// calls to Rpo128AbsorbColumnsAndPermute
struct Rpo128PartialState {
    Fp s0;
    Fp s1;
    Fp s2;
    Fp s3;
};

// RPO 128 digest output
struct Rpo128Digest {
    Fp e0;
    Fp e1;
    Fp e2;
    Fp e3;
};

// Pair of RPO 128 digests
struct Rpo128DigestPair {
    Rpo128Digest d0;
    Rpo128Digest d1;
};

// RPO parameters
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

// RPO constants used in the first half of each round
constant const Fp ROUND_CONSTANTS_0[STATE_WIDTH * NUM_ROUNDS] = {
    Fp(6936159699454947676), Fp(6871277616928621393), Fp(4226339945476756083), Fp(2261225084505152444), Fp(16808067423291017741), Fp(12862191241011323277), Fp(345720808813194915), Fp(10126368034161173654), Fp(840649715788759894), Fp(18155600607269645987), Fp(16577339120870559289), Fp(13749826054300849029),
    Fp(16047969944113931191), Fp(10474334246235299199), Fp(15773847146013662260), Fp(14401231158322525155), Fp(6009395255763488383), Fp(2108579439821148946), Fp(13820200715803196660), Fp(15968614366574245570), Fp(7529997729792773654), Fp(9429194013557833999), Fp(11639903126146281421), Fp(15759666882357935738),
    Fp(14807658266593669785), Fp(17258259860767641342), Fp(9534132615398591413), Fp(358719342502509866), Fp(7123090532818864651), Fp(734193187930710962), Fp(14873184913735487023), Fp(17965359964069906568), Fp(12664837478844326631), Fp(15575491070113731145), Fp(7221479899469196675), Fp(7328957460733188967),
    Fp(15088355010936495340), Fp(16762963605345901631), Fp(15278161326153175940), Fp(6257793333052173411), Fp(8418953127708045776), Fp(6523475766574412380), Fp(15192936988185261803), Fp(1578086224854546096), Fp(10840553425559156784), Fp(7453417405109536362), Fp(5173069484734008228), Fp(3284492202065476384),
    Fp(1724586709636399686), Fp(17997633752581871175), Fp(1284825320737914582), Fp(960534381847281815), Fp(6708901808183456837), Fp(8975591106768797316), Fp(52515315389099119), Fp(10009391031874081397), Fp(3091228317422201238), Fp(1063858230459024983), Fp(3396548655473917480), Fp(15046057790353688034),
    Fp(4867464583127666756), Fp(13816959924674544309), Fp(13931201815459591565), Fp(11494116713280125381), Fp(16823081743980874023), Fp(6760771226809185048), Fp(5346741505458044699), Fp(15124596060558844029), Fp(5332565678905773189), Fp(17640389307200936126), Fp(14049814539797608740), Fp(8882709539093378074),
    Fp(10507930462458090835), Fp(10669463960502417047), Fp(16753662827442720769), Fp(12967456627495301601), Fp(2989815121821278695), Fp(5894674479204135685), Fp(14187454698288462352), Fp(14795723369628125345), Fp(17260571099239679821), Fp(16009836214833755168), Fp(2009092225887788829), Fp(10838446069154019765),
};

// RPO constants used in the last half of each round
constant const Fp ROUND_CONSTANTS_1[STATE_WIDTH * NUM_ROUNDS] = {
    Fp(8939123259393952351), Fp(14708045228210488368), Fp(18125168669810517809), Fp(9309821433754818185), Fp(4714467145607136006), Fp(1302482025306688824), Fp(34829973686821040), Fp(5637233680011148778), Fp(227119480134509573), Fp(2530972937109017559), Fp(7210163798538732239), Fp(955913576003606833), 
    Fp(4449617297638325218), Fp(10843671682695268638), Fp(13198957499160452915), Fp(11541825028620451829), Fp(10963484480734735121), Fp(4752902142121643229), Fp(3015289210993491059), Fp(16344286514680205966), Fp(1811079964700766606), Fp(12735664961476037524), Fp(5775391330037813314), Fp(18223625362487900986), 
    Fp(7222477607687412281), Fp(4215615082079701144), Fp(6177508277476483691), Fp(3491362079220677263), Fp(10961785333913978630), Fp(1935408839283360916), Fp(13974192629927279950), Fp(18013556876298568088), Fp(7565676920589638093), Fp(9265825103386412558), Fp(8061587790235022972), Fp(6806849270604947860), 
    Fp(8066442548506952806), Fp(12791828131640457742), Fp(9268748809821748950), Fp(17496234860625277598), Fp(13583894547367420658), Fp(13920282495726802458), Fp(3933141341199584259), Fp(6658057712176150702), Fp(16812362035931029194), Fp(15160401867587809089), Fp(16411108749946146942), Fp(3390826434320009844), 
    Fp(18405475140095477472), Fp(13864039573264702148), Fp(496144052468360460), Fp(9791523668470936672), Fp(528582340156917005), Fp(15864481364569144493), Fp(682830611952089590), Fp(347158833826327515), Fp(13752775429919623417), Fp(10254722988306758482), Fp(8794150602427420596), Fp(2480344122229837853), 
    Fp(15462337562022968595), Fp(6729968753311049611), Fp(9250220857258211097), Fp(12031447985684644003), Fp(14538803180331344696), Fp(4055445230671851890), Fp(14764039661528567501), Fp(2047787218814287270), Fp(8977863094202715520), Fp(6560450968915612407), Fp(9976241128570886075), Fp(17877509887772213755), 
    Fp(3549624494907837709), Fp(4253629935471652443), Fp(2859199883984623807), Fp(1087607721547343649), Fp(7907517619951970198), Fp(11306402795121903516), Fp(10168009948206732524), Fp(9177440083248248246), Fp(13169036816957726187), Fp(12924186209140199217), Fp(9673006056831483321), Fp(747828276541750689)
};

inline void rpo_permute(threadgroup Fp* shared, unsigned local_state_offset) {
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
}

// // Rescue Prime Optimized hash function for 128 bit security: https://eprint.iacr.org/2022/1577.pdf
// [[ host_name("rpo_128_absorb_rows_and_permute_p18446744069414584321_fp") ]] kernel void 
// Rpo128AbsorbRowsAndPermute(constant Fp *matrix [[ buffer(0) ]],
//         device Rpo128PartialState *states [[ buffer(1) ]],
//         device Rpo128Digest *digests [[ buffer(2) ]],
//         constant unsigned &num_rows [[ buffer(3) ]],
//         constant unsigned &row_offset [[ buffer(4) ]],
//         constant unsigned &input_width [[ buffer(5) ]],
//         threadgroup Fp *shared [[ threadgroup(0) ]],
//         unsigned global_id [[ thread_position_in_grid ]],
//         unsigned local_id [[ thread_index_in_threadgroup ]]) {   
    


//     // unsigned local_state_offset = local_id * STATE_WIDTH * 2;

//     // // fetch state
//     // unsigned local_state_offset = local_id * STATE_WIDTH * 2;
//     // // faster than loading individual elements
//     // *((threadgroup Rpo128PartialState*) (shared + local_state_offset)) = states[global_id];
//     // // absorb the input into the state
//     // shared[local_state_offset + CAPACITY + 0] = col0[global_id];
//     // shared[local_state_offset + CAPACITY + 1] = col1[global_id];
//     // shared[local_state_offset + CAPACITY + 2] = col2[global_id];
//     // shared[local_state_offset + CAPACITY + 3] = col3[global_id];
//     // shared[local_state_offset + CAPACITY + 4] = col4[global_id];
//     // shared[local_state_offset + CAPACITY + 5] = col5[global_id];
//     // shared[local_state_offset + CAPACITY + 6] = col6[global_id];
//     // shared[local_state_offset + CAPACITY + 7] = col7[global_id];

//     rpo_permute(shared, local_state_offset);

//     // TODO: add flag to only write to one of these buffers
//     // redundant writes here are neglegable on performance <1%
//     digests[global_id] = *((threadgroup Rpo128Digest*) (shared + local_state_offset + CAPACITY));
//     states[global_id] = *((threadgroup Rpo128PartialState*) (shared + local_state_offset));
// }

// TODO: make functions generic over the hash function
// Rescue Prime Optimized hash function for 128 bit security: https://eprint.iacr.org/2022/1577.pdf
// Absorbs 8 columns of equal length. Hashes are generated row-wise.
[[ host_name("rpo_128_absorb_columns_and_permute_p18446744069414584321_fp") ]] kernel void 
Rpo128AbsorbColumnsAndPermute(constant Fp *col0 [[ buffer(0) ]],
        constant Fp *col1 [[ buffer(1) ]],
        constant Fp *col2 [[ buffer(2) ]],
        constant Fp *col3 [[ buffer(3) ]],
        constant Fp *col4 [[ buffer(4) ]],
        constant Fp *col5 [[ buffer(5) ]],
        constant Fp *col6 [[ buffer(6) ]],
        constant Fp *col7 [[ buffer(7) ]],
        device Rpo128PartialState *states [[ buffer(8) ]],
        device Rpo128Digest *digests [[ buffer(9) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {     
    // fetch state
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    // faster than loading individual elements
    *((threadgroup Rpo128PartialState*) (shared + local_state_offset)) = states[global_id];
    // absorb the input into the state
    shared[local_state_offset + CAPACITY + 0] = col0[global_id];
    shared[local_state_offset + CAPACITY + 1] = col1[global_id];
    shared[local_state_offset + CAPACITY + 2] = col2[global_id];
    shared[local_state_offset + CAPACITY + 3] = col3[global_id];
    shared[local_state_offset + CAPACITY + 4] = col4[global_id];
    shared[local_state_offset + CAPACITY + 5] = col5[global_id];
    shared[local_state_offset + CAPACITY + 6] = col6[global_id];
    shared[local_state_offset + CAPACITY + 7] = col7[global_id];

    rpo_permute(shared, local_state_offset);

    // TODO: add flag to only write to one of these buffers
    // redundant writes here are neglegable on performance <1%
    digests[global_id] = *((threadgroup Rpo128Digest*) (shared + local_state_offset + CAPACITY));
    states[global_id] = *((threadgroup Rpo128PartialState*) (shared + local_state_offset));
}

// Generates the first row of merkle tree after the leaf nodes
// using Rescue Prime Optimized hash function.
[[ host_name("rpo_128_gen_merkle_nodes_first_row_p18446744069414584321_fp") ]] kernel void 
Rpo128GenMerkleNodesFirstRow(constant Rpo128DigestPair *leaves [[ buffer(0) ]], 
        device Rpo128Digest *nodes [[ buffer(1) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {
    // fetch state
    // *((threadgroup Rpo128PartialState*) (shared + local_state_offset)) = { .s0 = Fp(0); .s1 = Fp(0); .s2 = Fp(0); .s3 = Fp(0) };
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    shared[local_state_offset + 0] = Fp(0);
    shared[local_state_offset + 1] = Fp(0);
    shared[local_state_offset + 2] = Fp(0);
    shared[local_state_offset + 3] = Fp(0);
    // absorb children as input
    *((threadgroup Rpo128DigestPair*) (shared + local_state_offset + CAPACITY)) = leaves[global_id];

    rpo_permute(shared, local_state_offset);

    // write digest
    nodes[N / 2 + global_id] = *((threadgroup Rpo128Digest*) (shared + local_state_offset + CAPACITY));
}

// Generates a row of merkle tree nodes using the Rescue Prime Optimized hash function.
[[ host_name("rpo_128_gen_merkle_nodes_row_p18446744069414584321_fp") ]] kernel void 
Rpo128GenMerkleNodesRow(device Rpo128Digest *nodes [[ buffer(0) ]],
        constant unsigned &round [[ buffer(1) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {
    // fetch state
    // *((threadgroup Rpo128PartialState*) (shared + local_state_offset)) = { .s0 = Fp(0); .s1 = Fp(0); .s2 = Fp(0); .s3 = Fp(0) };
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    shared[local_state_offset + 0] = Fp(0);
    shared[local_state_offset + 1] = Fp(0);
    shared[local_state_offset + 2] = Fp(0);
    shared[local_state_offset + 3] = Fp(0);
    // absorb children as input
    // TODO: c++ cast for readability
    *((threadgroup Rpo128DigestPair*) (shared + local_state_offset + CAPACITY)) = ((device Rpo128DigestPair*) nodes)[(N >> round) + global_id];

    rpo_permute(shared, local_state_offset);

    // write digest
    nodes[(N >> round) + global_id] = *((threadgroup Rpo128Digest*) (shared + local_state_offset + CAPACITY));
}

}

#endif /* hash_shaders_h */
