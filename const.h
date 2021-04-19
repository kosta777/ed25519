#ifndef CONST_H
#define CONST_H

#include <immintrin.h>
#if STDC
  #if R25
    #include "StdC-optimize/radix25/radix25.h"
  #else
    #include "StdC-optimize/radix/radix51.h"
  #endif
#elif BASE
  #include "straight/radix51.h"
#elif SSE
  #if R25
    #include "SSE/radix25/radix25.h"
  #else
    #include "SSE/radix/radix51.h"
  #endif
#endif
typedef u_int64_t ulli;

// b = 1<<8
// q = 1<<255 - 0x13
// l = 1<<252 + 0x14DEF9DEA2F79CD6 5812631A5CF5D3ED

const ulli ulli_b[] = {0x0,0x0,0x0,0x0,0x100};
const ulli ulli_Q[] = {0x0, 0x7fffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffed};
const ulli ulli_L[] = {0x0, 0x1000000000000000, 0x0000000000000000, 0x14def9dea2f79cd6, 0x5812631a5cf5d3ed};

/*
d = -121665 * inv(121666)     = 37095705934669439343138083508754565189542113879843219016388785533085940283555
I = expmod(2,(q-1)//4,q)      = 19681161376707505956807079304988542015446066515923890162744021073123829784752
By = 4 * inv(5)               = 15112221349535400772501151409588531511454012693041857206046113283949847762202
Bx = xrecover(By)             = 46316835694926478169428394003475163141307993866256225615783033603165251855960
B2 = [Bx % q,By % q]
*/
const ulli ulli_D[] = {0x0, 0x52036cee2b6ffe73, 0x8cc740797779e898, 0x00700a4d4141d8ab, 0x75eb4dca135978a3};
const ulli ulli_I[] = {0x0, 0x2b8324804fc1df0b, 0x2b4d00993dfbd7a7, 0x2f431806ad2fe478, 0xc4ee1b274a0ea0b0};
const ulli ulli_B[2][5] = {
  {0x0, 0x216936d3cd6e53fe, 0xc0a4e231fdd6dc5c, 0x692cc7609525a7b2, 0xc9562d608f25d51a},
  {0x0, 0x6666666666666666, 0x6666666666666666, 0x6666666666666666, 0x6666666666666658}
};

const ulli ulli_D_INV[] = {0x0, 0x40907ed214d5ce43, 0x2b162114cdb9cf66, 0x0b5dd6984279542e, 0x25e0f276cdc9f843}; // 29203323570479832313932222820853031755424015489030106335172286495629888845891


const ulli ulli_ZERO[] = {0x0, 0x0, 0x0, 0x0, 0x0};
const ulli ulli_ONE[] = {0x0, 0x0, 0x0, 0x0, 0x1};
const ulli ulli_TWO[] = {0x0, 0x0, 0x0, 0x0, 0x2};
const ulli ulli_THREE[] = {0x0, 0x0, 0x0, 0x0, 0x3};
const ulli ulli_2p254[] = {0x0, 0x4000000000000000, 0x0, 0x0, 0x0};


const radix_t CURVE_Q(ulli_Q);
const radix_t CURVE_b(ulli_b);
const radix_t CURVE_L(ulli_L);

const radix_t CURVE_D(ulli_D);
const radix_t CURVE_D_INV(ulli_D_INV);
const radix_t CURVE_I(ulli_I);
const radix_t CURVE_Bx(ulli_B[0]);
const radix_t CURVE_By(ulli_B[1]);
const radix_t CURVE_B2[] = {CURVE_Bx, CURVE_By};

const radix_t ZERO(ulli_ZERO);
const radix_t ONE(ulli_ONE);
const radix_t TWO(ulli_TWO);
const radix_t THREE(ulli_THREE);
const radix_t TWOp254(ulli_2p254);

#if R25
/************** RADIX 25.5 **********************/
const r25_t R25_CURVE_Q(ulli_Q);
const r25_t R25_CURVE_b(ulli_b);
const r25_t R25_CURVE_L(ulli_L);

const r25_t R25_CURVE_D(ulli_D);
const r25_t R25_CURVE_D_INV(ulli_D_INV);
const r25_t R25_CURVE_I(ulli_I);
const r25_t R25_CURVE_Bx(ulli_B[0]);
const r25_t R25_CURVE_By(ulli_B[1]);
const r25_t R25_CURVE_B2[] = {R25_CURVE_Bx, R25_CURVE_By};

const r25_t R25_ZERO(ulli_ZERO);
const r25_t R25_ONE(ulli_ONE);
const r25_t R25_TWO(ulli_TWO);
const r25_t R25_THREE(ulli_THREE);
const r25_t R25_TWOp254(ulli_2p254);

const __m256i _00 = {0, 0, 0, 0};
const __m256i _01 = {1, 1, 1, 1};
const __m256i _one[] = {_00, _00, _00, _00, _00, _00, _00, _00, _00, _01};
const __m256i _0 = {21499315, 21499315, 21499315, 21499315};
const __m256i _1 = {48412415, 48412415, 48412415, 48412415};
const __m256i _2 = {30306712, 30306712, 30306712, 30306712};
const __m256i _3 = {60832955, 60832955, 60832955, 60832955};
const __m256i _4 = {24766616, 24766616, 24766616, 24766616};
const __m256i _5 = {114729,   114729,   114729,   114729};
const __m256i _6 = {6949390,  6949390,  6949390,  6949390};
const __m256i _7 = {51736253, 51736253, 51736253, 51736253};
const __m256i _8 = {13857412, 13857412, 13857412, 13857412};
const __m256i _9 = {56195235, 56195235, 56195235, 56195235};
const __m256i _d[] = {_0, _1, _2, _3, _4, _5, _6, _7, _8, _9};
#endif
#endif
