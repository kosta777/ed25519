#include "radix25_mul.h"
#include <immintrin.h>

void _v_r25_mul(const v_limb_t* x, const v_limb_t* y, v_limb_t* z){
  __m256i c19 = _mm256_set1_epi64x(19);
  __m256i c38 = _mm256_set1_epi64x(38);
  __m256i c2 = _mm256_set1_epi64x(2);
  v_limb_t _x0 = x[0]; v_limb_t _19x0 = _mm256_mul_epu32(c19, _x0);   v_limb_t _38x0 = _mm256_mul_epu32(c38, _x0);
  v_limb_t _x1 = x[1]; v_limb_t _19x1 = _mm256_mul_epu32(c19, _x1);
  v_limb_t _x2 = x[2]; v_limb_t _2x2  = _mm256_add_epi32(_x2, _x2);   v_limb_t _19x2 = _mm256_mul_epu32(c19, _x2); v_limb_t _38x2 = _mm256_mul_epu32(c38, _x2);
  v_limb_t _x3 = x[3]; v_limb_t _19x3 = _mm256_mul_epu32(c19, _x3);
  v_limb_t _x4 = x[4]; v_limb_t _2x4  = _mm256_add_epi32(_x4, _x4);   v_limb_t _19x4 = _mm256_mul_epu32(c19, _x4); v_limb_t _38x4 = _mm256_mul_epu32(c38, _x4);
  v_limb_t _x5 = x[5]; v_limb_t _19x5 = _mm256_mul_epu32(c19, _x5);
  v_limb_t _x6 = x[6]; v_limb_t _2x6  = _mm256_add_epi32(_x6, _x6);   v_limb_t _19x6 = _mm256_mul_epu32(c19, _x6); v_limb_t _38x6 = _mm256_mul_epu32(c38, _x6);
  v_limb_t _x7 = x[7]; v_limb_t _19x7 = _mm256_mul_epu32(c19, _x7);
  v_limb_t _x8 = x[8]; v_limb_t _2x8  = _mm256_add_epi32(_x8, _x8);   v_limb_t _38x8 = _mm256_mul_epu32(c38, _x8);
  v_limb_t _x9 = x[9];
  v_limb_t _y0 = y[0];
  v_limb_t _y1 = y[1];
  v_limb_t _y2 = y[2];
  v_limb_t _y3 = y[3];
  v_limb_t _y4 = y[4];
  v_limb_t _y5 = y[5];
  v_limb_t _y6 = y[6];
  v_limb_t _y7 = y[7];
  v_limb_t _y8 = y[8];
  v_limb_t _y9 = y[9];

/* Each multiplication produces 58bits numbers.
*  Each limb is a summation of 10 products.
*  Split the sum in two to accumulate in two 64bits parts without overflow
*  Then sum together the two halves and carry to the next limb
*/
  // first limb
  v_limb_t _z9lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y9), _mm256_mul_epu32(_38x8, _y0)),    _mm256_add_epi64(_mm256_mul_epu32(_19x7, _y1), _mm256_mul_epu32(_38x6, _y2))),   _mm256_mul_epu32(_19x5, _y3));
  v_limb_t _z9lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_38x4, _y4), _mm256_mul_epu32(_19x3, _y5)),  _mm256_add_epi64(_mm256_mul_epu32(_38x2, _y6), _mm256_mul_epu32(_19x1, _y7))),   _mm256_mul_epu32(_38x0, _y8));
  // second
  v_limb_t _z8lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y8), _mm256_mul_epu32(_x8, _y9)),      _mm256_add_epi64(_mm256_mul_epu32(_19x7, _y0), _mm256_mul_epu32(_19x6, _y1))),   _mm256_mul_epu32(_19x5, _y2));
  v_limb_t _z8lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_19x4, _y3), _mm256_mul_epu32(_19x3, _y4)),  _mm256_add_epi64(_mm256_mul_epu32(_19x2, _y5), _mm256_mul_epu32(_19x1, _y6))),   _mm256_mul_epu32(_19x0, _y7));
  // third
  v_limb_t _z7lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y7), _mm256_mul_epu32(_2x8, _y8)),     _mm256_add_epi64(_mm256_mul_epu32(_x7, _y9), _mm256_mul_epu32(_38x6, _y0))),   _mm256_mul_epu32(_19x5, _y1));
  v_limb_t _z7lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_38x4, _y2), _mm256_mul_epu32(_19x3, _y3)),  _mm256_add_epi64(_mm256_mul_epu32(_38x2, _y4), _mm256_mul_epu32(_19x1, _y5))),   _mm256_mul_epu32(_38x0, _y6));
  // fourth
  v_limb_t _z6lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y6), _mm256_mul_epu32(_x8, _y7)),      _mm256_add_epi64(_mm256_mul_epu32(_x7, _y8), _mm256_mul_epu32(_x6, _y9))),     _mm256_mul_epu32(_19x5, _y0));
  v_limb_t _z6lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_19x4, _y1), _mm256_mul_epu32(_19x3, _y2)),  _mm256_add_epi64(_mm256_mul_epu32(_19x2, _y3), _mm256_mul_epu32(_19x1, _y4))),   _mm256_mul_epu32(_19x0, _y5));
  // fifth
  v_limb_t _z5lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y5), _mm256_mul_epu32(_2x8, _y6)),     _mm256_add_epi64(_mm256_mul_epu32(_x7, _y7), _mm256_mul_epu32(_2x6, _y8))),    _mm256_mul_epu32(_x5, _y9));
  v_limb_t _z5lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_38x4, _y0), _mm256_mul_epu32(_19x3, _y1)),  _mm256_add_epi64(_mm256_mul_epu32(_38x2, _y2), _mm256_mul_epu32(_19x1, _y3))),   _mm256_mul_epu32(_38x0, _y4));
  // sixth
  v_limb_t _z4lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y4), _mm256_mul_epu32(_x8, _y5)),      _mm256_add_epi64(_mm256_mul_epu32(_x7, _y6), _mm256_mul_epu32(_x6, _y7))),     _mm256_mul_epu32(_x5, _y8));
  v_limb_t _z4lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x4, _y9), _mm256_mul_epu32(_19x3, _y0)),  _mm256_add_epi64(_mm256_mul_epu32(_19x2, _y1), _mm256_mul_epu32(_19x1, _y2))),   _mm256_mul_epu32(_19x0, _y3));
  // seventh
  v_limb_t _z3lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y3), _mm256_mul_epu32(_2x8, _y4)),     _mm256_add_epi64(_mm256_mul_epu32(_x7, _y5), _mm256_mul_epu32(_2x6, _y6))),    _mm256_mul_epu32(_x5, _y7));
  v_limb_t _z3lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_2x4, _y8), _mm256_mul_epu32(_x3, _y9)),    _mm256_add_epi64(_mm256_mul_epu32(_38x2, _y0), _mm256_mul_epu32(_19x1, _y1))),   _mm256_mul_epu32(_38x0, _y2));
  // eighth
  v_limb_t _z2lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y2), _mm256_mul_epu32(_x8, _y3)),      _mm256_add_epi64(_mm256_mul_epu32(_x7, _y4), _mm256_mul_epu32(_x6, _y5))),     _mm256_mul_epu32(_x5, _y6));
  v_limb_t _z2lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x4, _y7), _mm256_mul_epu32(_x3, _y8)),    _mm256_add_epi64(_mm256_mul_epu32(_x2, _y9), _mm256_mul_epu32(_19x1, _y0))),   _mm256_mul_epu32(_19x0, _y1));
  // nineth
  v_limb_t _z1lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y1), _mm256_mul_epu32(_2x8, _y2)),     _mm256_add_epi64(_mm256_mul_epu32(_x7, _y3), _mm256_mul_epu32(_2x6, _y4))),    _mm256_mul_epu32(_x5, _y5));
  v_limb_t _z1lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_2x4, _y6), _mm256_mul_epu32(_x3, _y7)),    _mm256_add_epi64(_mm256_mul_epu32(_2x2, _y8), _mm256_mul_epu32(_x1, _y9))),     _mm256_mul_epu32(_38x0, _y0));
  // tenth
  v_limb_t _z0lo0 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x9, _y0), _mm256_mul_epu32(_x8, _y1)),      _mm256_add_epi64(_mm256_mul_epu32(_x7, _y2), _mm256_mul_epu32(_x6, _y3))),     _mm256_mul_epu32(_x5, _y4));
  v_limb_t _z0lo1 = _mm256_add_epi64(_mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(_x4, _y5), _mm256_mul_epu32(_x3, _y6)),    _mm256_add_epi64(_mm256_mul_epu32(_x2, _y7), _mm256_mul_epu32(_x1, _y8))),     _mm256_mul_epu32(_x0, _y9));

  // add two halves and store in limb. Saves carryover in next limb
  // u_int64_t c_in =0;
  // c_in = _addc64(_z9lo0, _z9lo1, &z[9], c_in);
  // c_in = _addc64(_z8lo0, _z8lo1, &z[8], c_in);
  // c_in = _addc64(_z7lo0, _z7lo1, &z[7], c_in);
  // c_in = _addc64(_z6lo0, _z6lo1, &z[6], c_in);
  // c_in = _addc64(_z5lo0, _z5lo1, &z[5], c_in);
  // c_in = _addc64(_z4lo0, _z4lo1, &z[4], c_in);
  // c_in = _addc64(_z3lo0, _z3lo1, &z[3], c_in);
  // c_in = _addc64(_z2lo0, _z2lo1, &z[2], c_in);
  // c_in = _addc64(_z1lo0, _z1lo1, &z[1], c_in);
  // c_in = _addc64(_z0lo0, _z0lo1, &z[0], c_in);
  // unsigned char c_in = 0;
  // c_in = _addcarry_u64(c_in, _z9lo0, _z9lo1, &z[9]);
  // c_in = _addcarry_u64(c_in, _z8lo0, _z8lo1, &z[8]);
  // c_in = _addcarry_u64(c_in, _z7lo0, _z7lo1, &z[7]);
  // c_in = _addcarry_u64(c_in, _z6lo0, _z6lo1, &z[6]);
  // c_in = _addcarry_u64(c_in, _z5lo0, _z5lo1, &z[5]);
  // c_in = _addcarry_u64(c_in, _z4lo0, _z4lo1, &z[4]);
  // c_in = _addcarry_u64(c_in, _z3lo0, _z3lo1, &z[3]);
  // c_in = _addcarry_u64(c_in, _z2lo0, _z2lo1, &z[2]);
  // c_in = _addcarry_u64(c_in, _z1lo0, _z1lo1, &z[1]);
  // c_in = _addcarry_u64(c_in, _z0lo0, _z0lo1, &z[0]);
  // TODO maybe add wrap around
  
  z[9] = _mm256_add_epi64(_z9lo0, _z9lo1);
  z[8] = _mm256_add_epi64(_z8lo0, _z8lo1);
  z[7] = _mm256_add_epi64(_z7lo0, _z7lo1);
  z[6] = _mm256_add_epi64(_z6lo0, _z6lo1);
  z[5] = _mm256_add_epi64(_z5lo0, _z5lo1);
  z[4] = _mm256_add_epi64(_z4lo0, _z4lo1);
  z[3] = _mm256_add_epi64(_z3lo0, _z3lo1);
  z[2] = _mm256_add_epi64(_z2lo0, _z2lo1);
  z[1] = _mm256_add_epi64(_z1lo0, _z1lo1);
  z[0] = _mm256_add_epi64(_z0lo0, _z0lo1);

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
