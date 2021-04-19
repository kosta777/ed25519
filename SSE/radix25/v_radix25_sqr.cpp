#include "radix25_sqr.h"
#include <immintrin.h>

void _v_r25_sqr(const v_limb_t *x, v_limb_t *z) {
  __m256i _19 = _mm256_set1_epi64x(19);
  __m256i _38 = _mm256_set1_epi64x(38);
//  __m256i _2 = _mm256_set1_epi64x(2);
  //                uses add to balance ports
  v_limb_t _x0 = x[0]; v_limb_t _2x0 = _mm256_add_epi64(_x0, _x0); v_limb_t _38x0 = _mm256_mul_epu32(_x0, _38);
  v_limb_t _x1 = x[1]; v_limb_t _2x1 = _mm256_add_epi64(_x1, _x1); v_limb_t _19x1 = _mm256_mul_epu32(_x1, _19);
  v_limb_t _x2 = x[2]; v_limb_t _2x2 = _mm256_add_epi64(_x2, _x2); v_limb_t _19x2 = _mm256_mul_epu32(_x2, _19); v_limb_t _38x2 = _mm256_mul_epu32(_x2, _38);
  v_limb_t _x3 = x[3]; v_limb_t _2x3 = _mm256_add_epi64(_x3, _x3); v_limb_t _19x3 = _mm256_mul_epu32(_x3, _19);
  v_limb_t _x4 = x[4]; v_limb_t _2x4 = _mm256_add_epi64(_x4, _x4); v_limb_t _19x4 = _mm256_mul_epu32(_x4, _19); v_limb_t _38x4 = _mm256_mul_epu32(_x4, _38);
  v_limb_t _x5 = x[5]; v_limb_t _2x5 = _mm256_add_epi64(_x5, _x5); v_limb_t _19x5 = _mm256_mul_epu32(_x5, _19);
  v_limb_t _x6 = x[6]; v_limb_t _2x6 = _mm256_add_epi64(_x6, _x6); v_limb_t _19x6 = _mm256_mul_epu32(_x6, _19); v_limb_t _38x6 = _mm256_mul_epu32(_x6, _38);
  v_limb_t _x7 = x[7]; v_limb_t _2x7 = _mm256_add_epi64(_x7, _x7); v_limb_t _19x7 = _mm256_mul_epu32(_x7, _19);
  v_limb_t _x8 = x[8]; v_limb_t _2x8 = _mm256_add_epi64(_x8, _x8); v_limb_t _38x8 = _mm256_mul_epu32(_x8, _38);
  v_limb_t _x9 = x[9]; v_limb_t _2x9 = _mm256_add_epi64(_x9, _x9);

/* Each multiplication produces 58bits numbers.
*  Each limb is a summation of 5-6 products.
*  No nee to split the sum this time as it should always fit in 63 bits
*/
//  z[9] = _x9*_x9  + _38x8*_2x0  + _19x7*_2x1 + _38x6*_2x2 + _19x5*_2x3 + _38x4*_x4;
  z[9] = _mm256_add_epi64(_mm256_mul_epu32(_x9, _x9),
         _mm256_add_epi64(_mm256_mul_epu32(_38x8, _2x0),
         _mm256_add_epi64(_mm256_mul_epu32(_19x7, _2x1),
         _mm256_add_epi64(_mm256_mul_epu32(_38x6, _2x2),
         _mm256_add_epi64(_mm256_mul_epu32(_19x5, _2x3),
                         (_mm256_mul_epu32(_38x4, _x4)))))));
//  z[8] = _2x9*_x8 + _19x7*_2x0  + _19x6*_2x1 + _19x5*_2x2 + _19x4*_2x3;
  z[8] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x8),
         _mm256_add_epi64(_mm256_mul_epu32(_19x7, _2x0),
         _mm256_add_epi64(_mm256_mul_epu32(_19x6, _2x1),
         _mm256_add_epi64(_mm256_mul_epu32(_19x5, _2x2),
                         (_mm256_mul_epu32(_19x4, _2x3))))));
//  z[7] = _2x9*_x7 + _2x8*_x8    + _38x6*_2x0 + _19x5*_2x1 + _38x4*_2x2 + _19x3*_x3;
  z[7] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x7),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _x8),
         _mm256_add_epi64(_mm256_mul_epu32(_38x6, _2x0),
         _mm256_add_epi64(_mm256_mul_epu32(_19x5, _2x1),
         _mm256_add_epi64(_mm256_mul_epu32(_38x4, _2x2),
                         (_mm256_mul_epu32(_19x3, _x3)))))));
//  z[6] = _2x9*_x6 + _2x8*_x7    + _19x5*_2x0 + _19x4*_2x1 + _19x3*_2x2;
  z[6] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x6),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _x7),
         _mm256_add_epi64(_mm256_mul_epu32(_19x5, _2x0),
         _mm256_add_epi64(_mm256_mul_epu32(_19x4, _2x1),
                         (_mm256_mul_epu32(_19x3, _2x2))))));
//  z[5] = _2x9*_x5 + _2x8*_2x6   + _x7*_x7    + _38x4*_2x0 + _19x3*_2x1 + _38x2*_x2;
  z[5] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x5),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _2x6),
         _mm256_add_epi64(_mm256_mul_epu32(_x7, _x7),
         _mm256_add_epi64(_mm256_mul_epu32(_38x4, _2x0),
         _mm256_add_epi64(_mm256_mul_epu32(_19x3, _2x1),
                         (_mm256_mul_epu32(_38x2, _x2)))))));
//  z[4] = _2x9*_x4 + _2x8*_x5    + _2x7*_x6   + _19x3*_2x0 + _19x2*_2x1;
  z[4] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x4),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _x5),
         _mm256_add_epi64(_mm256_mul_epu32(_2x7, _x6),
         _mm256_add_epi64(_mm256_mul_epu32(_19x3, _2x0),
                         (_mm256_mul_epu32(_19x2, _2x1))))));
//  z[3] = _2x9*_x3 + _2x8*_2x4   + _2x7*_x5   + _2x6*_x6   + _38x2*_2x0 + _19x1*_x1;
  z[3] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x3),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _2x4),
         _mm256_add_epi64(_mm256_mul_epu32(_2x7, _x5),
         _mm256_add_epi64(_mm256_mul_epu32(_2x6, _x6),
         _mm256_add_epi64(_mm256_mul_epu32(_38x2, _2x0),
                         (_mm256_mul_epu32(_19x1, _x1)))))));
//  z[2] = _2x9*_x2 + _2x8*_x3    + _2x7*_x4   + _2x6*_x5   + _19x1*_2x0;
  z[2] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x2),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _x3),
         _mm256_add_epi64(_mm256_mul_epu32(_2x7, _x4),
         _mm256_add_epi64(_mm256_mul_epu32(_2x6, _x5),
                         (_mm256_mul_epu32(_19x1, _2x0))))));
//  z[1] = _2x9*_x1 + _2x8*_2x2   + _2x7*_x3   + _2x6*_2x4  + _x5*_x5    + _38x0*_x0;
  z[1] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x1),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _2x2),
         _mm256_add_epi64(_mm256_mul_epu32(_2x7, _x3),
         _mm256_add_epi64(_mm256_mul_epu32(_2x6, _2x4),
         _mm256_add_epi64(_mm256_mul_epu32(_x5, _x5),
                         (_mm256_mul_epu32(_38x0, _x0)))))));
//  z[0] = _2x9*_x0 + _2x8*_x1    + _2x7*_x2   + _2x6*_x3   + _2x5*_x4;
  z[0] = _mm256_add_epi64(_mm256_mul_epu32(_2x9, _x0),
         _mm256_add_epi64(_mm256_mul_epu32(_2x8, _x1),
         _mm256_add_epi64(_mm256_mul_epu32(_2x7, _x2),
         _mm256_add_epi64(_mm256_mul_epu32(_2x6, _x3),
                         (_mm256_mul_epu32(_2x5, _x4))))));
  
  // _carry_and_normalize(z);
  _v_carry_and_normalize(z);
  /* at this point the limb9 is 43bits long, since limb0 had a 38bits carry, then multiplied by 19 becomes 42bits
  *  42 bits + 26bits = 43 bits
  * an additional carry limb9->limb8 of 43-26=17bits
   */
  v_limb_t _t9 = z[9]; v_limb_t _t8 = z[8];
  v_limb_t _r90 = _mm256_srli_epi64(_t9, 26);
  z[8] = _mm256_add_epi64(_t8, _r90); //z8
  __m256i _mask = _mm256_set1_epi64x(RADIXMASK26);
  z[9] = _mm256_and_si256(_t9, _mask);
  
  // _v_fix_degenerate(z);
}
