#ifndef hash_shaders_h
#define hash_shaders_h

#include <metal_stdlib>
#include "felt_u64.h.metal"
using namespace metal;

namespace p18446744069414584321 {

// Stores the first CAPACITY many items of the state.
// This is the only state that needs to be persisted between
// calls to Rpo256AbsorbColumnsAndPermute
struct Rpo256PartialState {
    Fp s0;
    Fp s1;
    Fp s2;
    Fp s3;
};

// RPO 128 digest output
struct Rpo256Digest {
    Fp e0;
    Fp e1;
    Fp e2;
    Fp e3;
};

// Pair of RPO 128 digests
struct Rpo256DigestPair {
    Rpo256Digest d0;
    Rpo256Digest d1;
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

// TODO: use pair from standard library (can't figure out how to import)
template <class t1, class t2>
struct pair {
  t1 a;
  t2 b;
};

inline ulong2 ifft2_real(long2 in) {
    return ulong2((ulong) (in.x + in.y), (ulong) (in.x - in.y));
}

inline ulong4 ifft4_real(long4 in) {
    ulong2 z0 = ifft2_real(long2(in.x + in.w, in.y));
    ulong2 z1 = ifft2_real(long2(in.x - in.w, -in.z));
    return ulong4(z0.x, z1.x, z0.y, z1.y);
}

inline long2 fft2_real(ulong2 in) {
    return long2((long) (in.x + in.y), (long) in.x - (long) in.y);
}

inline long4 fft4_real(ulong4 in) {
    long2 z0 = fft2_real(ulong2(in.x, in.z));
    long2 z1 = fft2_real(ulong2(in.y, in.w));
    return long4(z0.x + z1.x, z0.y, -z1.y, z0.x - z1.x);
}

constant const long3 MDS_FREQ_BLOCK_ONE = long3(16, 8, 16);
constant const pair<long3, long3> MDS_FREQ_BLOCK_TWO = { .a = long3(-1, -1, 4), .b = long3(2, 1, 8) };
constant const long3 MDS_FREQ_BLOCK_THREE = long3(-8, 1, 1);

inline long3 block1(long3 in) {
    return long3(
        in.x * MDS_FREQ_BLOCK_ONE.x + in.y * MDS_FREQ_BLOCK_ONE.z + in.z * MDS_FREQ_BLOCK_ONE.y,
        in.x * MDS_FREQ_BLOCK_ONE.y + in.y * MDS_FREQ_BLOCK_ONE.x + in.z * MDS_FREQ_BLOCK_ONE.z,
        in.x * MDS_FREQ_BLOCK_ONE.z + in.y * MDS_FREQ_BLOCK_ONE.y + in.z * MDS_FREQ_BLOCK_ONE.x
    );
}

inline pair<long3, long3> block2(pair<long3, long3> in) {
    long x0s = in.a.x + in.b.x;
    long x1s = in.a.y + in.b.y;
    long x2s = in.a.z + in.b.z;
    long y0s = MDS_FREQ_BLOCK_TWO.a.x + MDS_FREQ_BLOCK_TWO.b.x;
    long y1s = MDS_FREQ_BLOCK_TWO.a.y + MDS_FREQ_BLOCK_TWO.b.y;
    long y2s = MDS_FREQ_BLOCK_TWO.a.z + MDS_FREQ_BLOCK_TWO.b.z;

    // Compute x0​y0 ​− ix1​y2​ − ix2​y1​ using Karatsuba for complex numbers multiplication
    long2 m0 = long2(in.a.x * MDS_FREQ_BLOCK_TWO.a.x, in.b.x * MDS_FREQ_BLOCK_TWO.b.x);
    long2 m1 = long2(in.a.y * MDS_FREQ_BLOCK_TWO.a.z, in.b.y * MDS_FREQ_BLOCK_TWO.b.z);
    long2 m2 = long2(in.a.z * MDS_FREQ_BLOCK_TWO.a.y, in.b.z * MDS_FREQ_BLOCK_TWO.b.y);
    long z0r = (m0.x - m0.y) + (x1s * y2s - m1.x - m1.y) + (x2s * y1s - m2.x - m2.y);
    long z0i = (x0s * y0s - m0.x - m0.y) + (-m1.x + m1.y) + (-m2.x + m2.y);
    long2 z0 = long2(z0r, z0i);

    // Compute x0​y1​ + x1​y0​ − ix2​y2 using Karatsuba for complex numbers multiplication
    m0 = long2(in.a.x * MDS_FREQ_BLOCK_TWO.a.y, in.b.x * MDS_FREQ_BLOCK_TWO.b.y);
    m1 = long2(in.a.y * MDS_FREQ_BLOCK_TWO.a.x, in.b.y * MDS_FREQ_BLOCK_TWO.b.x);
    m2 = long2(in.a.z * MDS_FREQ_BLOCK_TWO.a.z, in.b.z * MDS_FREQ_BLOCK_TWO.b.z);
    long z1r = (m0.x - m0.y) + (m1.x - m1.y) + (x2s * y2s - m2.x - m2.y);
    long z1i = (x0s * y1s - m0.x - m0.y) + (x1s * y0s - m1.x - m1.y) + (-m2.x + m2.y);
    long2 z1 = long2(z1r, z1i);

    // Compute x0​y2​ + x1​y1 ​+ x2​y0​ using Karatsuba for complex numbers multiplication
    m0 = long2(in.a.x * MDS_FREQ_BLOCK_TWO.a.z, in.b.x * MDS_FREQ_BLOCK_TWO.b.z);
    m1 = long2(in.a.y * MDS_FREQ_BLOCK_TWO.a.y, in.b.y * MDS_FREQ_BLOCK_TWO.b.y);
    m2 = long2(in.a.z * MDS_FREQ_BLOCK_TWO.a.x, in.b.z * MDS_FREQ_BLOCK_TWO.b.x);
    long z2r = (m0.x - m0.y) + (m1.x - m1.y) + (m2.x - m2.y);
    long z2i = (x0s * y2s - m0.x - m0.y) + (x1s * y1s - m1.x - m1.y) + (x2s * y0s - m2.x - m2.y);
    long2 z2 = long2(z2r, z2i);

    return { .a = long3(z0.x, z1.x, z2.x), .b = long3(z0.y, z1.y, z2.y) };
}

inline long3 block3(long3 in) {
    return long3(
        in.x * MDS_FREQ_BLOCK_THREE.x - in.y * MDS_FREQ_BLOCK_THREE.z - in.z * MDS_FREQ_BLOCK_THREE.y,
        in.x * MDS_FREQ_BLOCK_THREE.y + in.y * MDS_FREQ_BLOCK_THREE.x - in.z * MDS_FREQ_BLOCK_THREE.z,
        in.x * MDS_FREQ_BLOCK_THREE.z + in.y * MDS_FREQ_BLOCK_THREE.y + in.z * MDS_FREQ_BLOCK_THREE.x
    );
}

// Adapted from Miden: 
// https://github.com/0xPolygonMiden/crypto/blob/main/src/hash/rpo/mds_freq.rs
inline void mds_multiply_freq(unsigned long state[STATE_WIDTH]) {
    long4 u0 = fft4_real(ulong4(state[0], state[3], state[6], state[9]));
    long4 u1 = fft4_real(ulong4(state[1], state[4], state[7], state[10]));
    long4 u2 = fft4_real(ulong4(state[2], state[5], state[8], state[11]));

    long3 v0 = block1(long3(u0.x, u1.x, u2.x));
    pair<long3, long3> v1 = block2({ .a = long3(u0.y, u1.y, u2.y), .b = long3(u0.z, u1.z, u2.z) });
    long3 v2 = block3(long3(u0.w, u1.w, u2.w));

    ulong4 s0 = ifft4_real(long4(v0.x, v1.a.x, v1.b.x, v2.x));
    ulong4 s1 = ifft4_real(long4(v0.y, v1.a.y, v1.b.y, v2.y));
    ulong4 s2 = ifft4_real(long4(v0.z, v1.a.z, v1.b.z, v2.z));

    state[0] = s0.x;
    state[1] = s1.x;
    state[2] = s2.x;
    state[3] = s0.y;
    state[4] = s1.y;
    state[5] = s2.y;
    state[6] = s0.z;
    state[7] = s1.z;
    state[8] = s2.z;
    state[9] = s0.w;
    state[10] = s1.w;
    state[11] = s2.w;
}

inline void apply_mds_freq(threadgroup Fp* shared, unsigned local_state_offset) {
    unsigned long state_l[STATE_WIDTH];
    unsigned long state_h[STATE_WIDTH];

#pragma unroll
    for (unsigned j = 0; j < STATE_WIDTH; j++) {
        Fp element = shared[local_state_offset + j];
        unsigned long s = (unsigned long) element;
        state_l[j] = s & 0xFFFFFFFF;
        state_h[j] = s >> 32;
    }

    mds_multiply_freq(state_l);
    mds_multiply_freq(state_h);

#pragma unroll
    for (unsigned j = 0; j < STATE_WIDTH; j++) {
        u128 s = u128(state_l[j]) + (u128(state_h[j]) << 32);
        unsigned long z = (s.high << 32) - s.high;
        unsigned overflow = s.low > (0xFFFFFFFFFFFFFFFF - z);
        unsigned long res = s.low + z;
        unsigned adj = -overflow;
        shared[local_state_offset + j] = Fp(res + adj);
    }
}

inline void rpo_permute_freq(threadgroup Fp* shared, unsigned local_state_offset) {
#pragma unroll
    for (unsigned i = 0; i < NUM_ROUNDS; i++) {
        apply_mds_freq(shared, local_state_offset);
        for (unsigned j = 0; j < STATE_WIDTH; j++) {
            Fp new_val = shared[local_state_offset + j];
            new_val = new_val + ROUND_CONSTANTS_0[i * STATE_WIDTH + j];
            shared[local_state_offset + j] = new_val.pow7();
        }

        apply_mds_freq(shared, local_state_offset);
        for (unsigned j = 0; j < STATE_WIDTH; j++) {
            Fp new_val = shared[local_state_offset + j];
            new_val = new_val + ROUND_CONSTANTS_1[i * STATE_WIDTH + j];
            shared[local_state_offset + j] = new_val.pow10540996611094048183();
        }
    }
}

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

// TODO: make functions generic over the hash function
// Rescue Prime Optimized hash function for 128 bit security: https://eprint.iacr.org/2022/1577.pdf
// Absorbs 8 columns of equal length. Hashes are generated row-wise.
[[ host_name("rpo_256_absorb_columns_and_permute_p18446744069414584321_fp") ]] kernel void 
Rpo256AbsorbColumnsAndPermute(constant Fp *col0 [[ buffer(0) ]],
        constant Fp *col1 [[ buffer(1) ]],
        constant Fp *col2 [[ buffer(2) ]],
        constant Fp *col3 [[ buffer(3) ]],
        constant Fp *col4 [[ buffer(4) ]],
        constant Fp *col5 [[ buffer(5) ]],
        constant Fp *col6 [[ buffer(6) ]],
        constant Fp *col7 [[ buffer(7) ]],
        device Rpo256PartialState *states [[ buffer(8) ]],
        device Rpo256Digest *digests [[ buffer(9) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {     
    // load hasher state
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    *((threadgroup Rpo256PartialState*) (shared + local_state_offset)) = states[global_id];
    // absorb the input into the state
    shared[local_state_offset + CAPACITY + 0] = col0[global_id];
    shared[local_state_offset + CAPACITY + 1] = col1[global_id];
    shared[local_state_offset + CAPACITY + 2] = col2[global_id];
    shared[local_state_offset + CAPACITY + 3] = col3[global_id];
    shared[local_state_offset + CAPACITY + 4] = col4[global_id];
    shared[local_state_offset + CAPACITY + 5] = col5[global_id];
    shared[local_state_offset + CAPACITY + 6] = col6[global_id];
    shared[local_state_offset + CAPACITY + 7] = col7[global_id];

    rpo_permute_freq(shared, local_state_offset);

    // TODO: add flag to only write to one of these buffers
    // redundant writes here are neglegable on performance <1%
    digests[global_id] = *((threadgroup Rpo256Digest*) (shared + local_state_offset + CAPACITY));
    states[global_id] = *((threadgroup Rpo256PartialState*) (shared + local_state_offset));
}

// TODO: make functions generic over the hash function
// Rescue Prime Optimized hash function for 128 bit security: https://eprint.iacr.org/2022/1577.pdf
// Absorbs 8 column of equal length in row major order. Hashes are generated row-wise.
[[ host_name("rpo_256_absorb_rows_and_permute_p18446744069414584321_fp") ]] kernel void 
Rpo256AbsorbRowsAndPermute(constant Fp *rows [[ buffer(0) ]],
        device Rpo256PartialState *states [[ buffer(1) ]],
        device Rpo256Digest *digests [[ buffer(2) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned tg_size [[ threads_per_threadgroup ]],
        unsigned tg_id [[ threadgroup_position_in_grid ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {     
    // load hasher state
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    *((threadgroup Rpo256PartialState*) (shared + local_state_offset)) = states[global_id];
    // absorb the input into the state (done like this for coalleced reads)
    unsigned tg_offset = tg_size * tg_id * 8;
    for (unsigned i = 0; i < 8; i++) {
        unsigned local_offset = tg_size * i + local_id;
        unsigned hasher_id = local_offset / 8;
        unsigned hasher_offset = local_offset % 8;
        shared[hasher_id * STATE_WIDTH * 2 + CAPACITY + hasher_offset] = rows[tg_offset + local_offset];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    rpo_permute_freq(shared, local_state_offset);

    // TODO: add flag to only write to one of these buffers
    // redundant writes here are neglegable on performance <1%
    digests[global_id] = *((threadgroup Rpo256Digest*) (shared + local_state_offset + CAPACITY));
    states[global_id] = *((threadgroup Rpo256PartialState*) (shared + local_state_offset));
}

// Generates the first row of merkle tree after the leaf nodes
// using Rescue Prime Optimized hash function.
[[ host_name("rpo_128_gen_merkle_nodes_first_row_p18446744069414584321_fp") ]] kernel void 
Rpo256GenMerkleNodesFirstRow(constant Rpo256DigestPair *leaves [[ buffer(0) ]], 
        device Rpo256Digest *nodes [[ buffer(1) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {
    // fetch state
    // *((threadgroup Rpo256PartialState*) (shared + local_state_offset)) = { .s0 = Fp(0); .s1 = Fp(0); .s2 = Fp(0); .s3 = Fp(0) };
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    shared[local_state_offset + 0] = Fp(0);
    shared[local_state_offset + 1] = Fp(0);
    shared[local_state_offset + 2] = Fp(0);
    shared[local_state_offset + 3] = Fp(0);
    // absorb children as input
    *((threadgroup Rpo256DigestPair*) (shared + local_state_offset + CAPACITY)) = leaves[global_id];

    rpo_permute(shared, local_state_offset);

    // write digest
    nodes[N / 2 + global_id] = *((threadgroup Rpo256Digest*) (shared + local_state_offset + CAPACITY));
}

// Generates a row of merkle tree nodes using the Rescue Prime Optimized hash function.
[[ host_name("rpo_128_gen_merkle_nodes_row_p18446744069414584321_fp") ]] kernel void 
Rpo256GenMerkleNodesRow(device Rpo256Digest *nodes [[ buffer(0) ]],
        constant unsigned &round [[ buffer(1) ]],
        threadgroup Fp *shared [[ threadgroup(0) ]],
        unsigned global_id [[ thread_position_in_grid ]],
        unsigned local_id [[ thread_index_in_threadgroup ]]) {
    // fetch state
    // *((threadgroup Rpo256PartialState*) (shared + local_state_offset)) = { .s0 = Fp(0); .s1 = Fp(0); .s2 = Fp(0); .s3 = Fp(0) };
    unsigned local_state_offset = local_id * STATE_WIDTH * 2;
    shared[local_state_offset + 0] = Fp(0);
    shared[local_state_offset + 1] = Fp(0);
    shared[local_state_offset + 2] = Fp(0);
    shared[local_state_offset + 3] = Fp(0);
    // absorb children as input
    // TODO: c++ cast for readability
    *((threadgroup Rpo256DigestPair*) (shared + local_state_offset + CAPACITY)) = ((device Rpo256DigestPair*) nodes)[(N >> round) + global_id];

    rpo_permute(shared, local_state_offset);

    // write digest
    nodes[(N >> round) + global_id] = *((threadgroup Rpo256Digest*) (shared + local_state_offset + CAPACITY));
}

}

#endif /* hash_shaders_h */
