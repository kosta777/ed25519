#include "radix25_sqr.h"

void _r25_sqr(const limb_t* x, limb_t* z){
//  //                uses add to balance ports
//  limb_t _x0 = x[0]; limb_t _2x0 = _x0+_x0; limb_t _38x0 = 38*_x0;
//  limb_t _x1 = x[1]; limb_t _2x1 = _x1+_x1; limb_t _19x1 = 19*_x1;
//  limb_t _x2 = x[2]; limb_t _2x2 = _x2+_x2; limb_t _19x2 = 19*_x2; limb_t _38x2 = 38*_x2;
//  limb_t _x3 = x[3]; limb_t _2x3 = _x3+_x3; limb_t _19x3 = 19*_x3;
//  limb_t _x4 = x[4]; limb_t _2x4 = _x4+_x4; limb_t _19x4 = 19*_x4; limb_t _38x4 = 38*_x4;
//  limb_t _x5 = x[5]; limb_t _2x5 = _x5+_x5; limb_t _19x5 = 19*_x5;
//  limb_t _x6 = x[6]; limb_t _2x6 = _x6+_x6; limb_t _19x6 = 19*_x6; limb_t _38x6 = 38*_x6;
//  limb_t _x7 = x[7]; limb_t _2x7 = _x7+_x7; limb_t _19x7 = 19*_x7;
//  limb_t _x8 = x[8]; limb_t _2x8 = _x8+_x8; limb_t _38x8 = 38*_x8;
//  limb_t _x9 = x[9]; limb_t _2x9 = _x9+_x9;
  
  __m256i _jjj = _mm256_set_epi64x(19, 19, 19, 38);
  __m256i _jjjj = _mm256_set_epi64x(19, 19, 19, 19);
  __m256i _j = _mm256_set_epi64x(38, 38, 38, 0);
  
  __m256i _x0123 = _mm256_loadu_si256((__m256i const *)&x[0]);
  __m256i _x4567 = _mm256_loadu_si256((__m256i const *)&x[4]);
  __m128i _x89 = _mm_loadu_si128((__m128i const *)&x[8]);
  __m256i _x0246 = _mm256_set_epi64x(x[6], x[4], x[2], 0);
  
  __m256i _2x0123 = _mm256_add_epi32(_x0123, _x0123);   __m256i _jx0123 = _mm256_mul_epu32(_x0123, _jjj);   __m256i _jx0246 = _mm256_mul_epu32(_x0246, _j);
  __m256i _2x4567 = _mm256_add_epi32(_x4567, _x4567);   __m256i _jx4567 = _mm256_mul_epu32(_x4567, _jjjj);
  __m128i _2x89 = _mm_add_epi32(_x89, _x89);
  
  
  limb_t _x0 = x[0]; limb_t _2x0 = _2x0123[0]; limb_t _38x0 = _jx0123[0];
  limb_t _x1 = x[1]; limb_t _2x1 = _2x0123[1]; limb_t _19x1 = _jx0123[1];
  limb_t _x2 = x[2]; limb_t _2x2 = _2x0123[2]; limb_t _19x2 = _jx0123[2]; limb_t _38x2 = _jx0246[1];
  limb_t _x3 = x[3]; limb_t _2x3 = _2x0123[3]; limb_t _19x3 = _jx0123[3];
  limb_t _x4 = x[4]; limb_t _2x4 = _2x4567[0]; limb_t _19x4 = _jx4567[0]; limb_t _38x4 = _jx0246[2];
  limb_t _x5 = x[5]; limb_t _2x5 = _2x4567[1]; limb_t _19x5 = _jx4567[1];
  limb_t _x6 = x[6]; limb_t _2x6 = _2x4567[2]; limb_t _19x6 = _jx4567[2]; limb_t _38x6 = _jx0246[3];
  limb_t _x7 = x[7]; limb_t _2x7 = _2x4567[3]; limb_t _19x7 = _jx4567[3];
  limb_t _x8 = x[8]; limb_t _2x8 = _2x89[0];   limb_t _38x8 = 38*_x8;
  limb_t _x9 = x[9]; limb_t _2x9 = _2x89[1];
  
  

  
//  z[0] = _2x9*_x0 + _2x8*_x1    + _2x7*_x2   + _2x6*_x3   + _2x5*_x4;
//  z[1] = _2x9*_x1 + _2x8*_2x2   + _2x7*_x3   + _2x6*_2x4  + _x5*_x5    + _38x0*_x0;
//  z[2] = _2x9*_x2 + _2x8*_x3    + _2x7*_x4   + _2x6*_x5   + _19x1*_2x0;
//  z[3] = _2x9*_x3 + _2x8*_2x4   + _2x7*_x5   + _2x6*_x6   + _38x2*_2x0 + _19x1*_x1;
//  z[4] = _2x9*_x4 + _2x8*_x5    + _2x7*_x6   + _19x3*_2x0 + _19x2*_2x1;
//  z[5] = _2x9*_x5 + _2x8*_2x6   + _x7*_x7    + _38x4*_2x0 + _19x3*_2x1 + _38x2*_x2;
//  z[6] = _2x9*_x6 + _2x8*_x7    + _19x5*_2x0 + _19x4*_2x1 + _19x3*_2x2;
//  z[7] = _2x9*_x7 + _2x8*_x8    + _38x6*_2x0 + _19x5*_2x1 + _38x4*_2x2 + _19x3*_x3;
//  z[8] = _2x9*_x8 + _19x7*_2x0  + _19x6*_2x1 + _19x5*_2x2 + _19x4*_2x3;
//  z[9] = _x9*_x9  + _38x8*_2x0  + _19x7*_2x1 + _38x6*_2x2 + _19x5*_2x3 + _38x4*_x4;
  
  __m256i _000 = _mm256_set_epi64x(_2x9, _2x9, _2x9, _2x9);     __m256i _001 = _mm256_set_epi64x(_x3, _x2, _x1, _x0);
  __m256i _010 = _mm256_set_epi64x(_2x9, _2x9, _2x9, _2x9);     __m256i _011 = _mm256_set_epi64x(_x7, _x6, _x5, _x4);
  __m128i _020 = _mm_set_epi64x(_x9, _2x9);                     __m128i _021 = _mm_set_epi64x(_x9, _x8);
  
  __m256i _100 = _mm256_set_epi64x(_2x8, _2x8, _2x8, _2x8);     __m256i _101 = _mm256_set_epi64x(_2x4, _x3, _2x2, _x1);
  __m256i _110 = _mm256_set_epi64x(_2x8, _2x8, _2x8, _2x8);     __m256i _111 = _mm256_set_epi64x(_x8, _x7, _2x6, _x5);
  __m128i _120 = _mm_set_epi64x(_38x8, _19x7);                  __m128i _121 = _mm_set_epi64x(_2x0, _2x0);
  
  __m256i _200 = _mm256_set_epi64x(_2x7, _2x7, _2x7, _2x7);     __m256i _201 = _mm256_set_epi64x(_x5, _x4, _x3, _x2);
  __m256i _210 = _mm256_set_epi64x(_38x6, _19x5, _x7, _2x7);    __m256i _211 = _mm256_set_epi64x(_2x0, _2x0, _x7, _x6);
  __m128i _220 = _mm_set_epi64x(_19x7, _19x6);                  __m128i _221 = _mm_set_epi64x(_2x1, _2x1);
  
  __m256i _300 = _mm256_set_epi64x(_2x6, _2x6, _2x6, _2x6);     __m256i _301 = _mm256_set_epi64x(_x6, _x5, _2x4, _x3);
  __m256i _310 = _mm256_set_epi64x(_19x5, _19x4, _38x4, _19x3); __m256i _311 = _mm256_set_epi64x(_2x1, _2x1, _2x0, _2x0);
  __m128i _320 = _mm_set_epi64x(_38x6, _19x5);                  __m128i _321 = _mm_set_epi64x(_2x2, _2x2);
  
  __m256i _400 = _mm256_set_epi64x(_38x2, _19x1, _x5, _2x5);    __m256i _401 = _mm256_set_epi64x(_2x0, _2x0, _x5, _x4);
  __m256i _410 = _mm256_set_epi64x(_38x4, _19x3, _19x3, _19x2); __m256i _411 = _mm256_set_epi64x(_2x2, _2x2, _2x1, _2x1);
  __m128i _420 = _mm_set_epi64x(_19x5, _19x4);                  __m128i _421 = _mm_set_epi64x(_2x3, _2x3);
  
  __m256i _prod00 = _mm256_mul_epu32(_000, _001);   __m256i _prod01 = _mm256_mul_epu32(_010, _011);   __m128i _prod02 = _mm_mul_epu32(_020, _021);
  __m256i _prod10 = _mm256_mul_epu32(_100, _101);   __m256i _prod11 = _mm256_mul_epu32(_110, _111);   __m128i _prod12 = _mm_mul_epu32(_120, _121);
  __m256i _prod20 = _mm256_mul_epu32(_200, _201);   __m256i _prod21 = _mm256_mul_epu32(_210, _211);   __m128i _prod22 = _mm_mul_epu32(_220, _221);
  __m256i _prod30 = _mm256_mul_epu32(_300, _301);   __m256i _prod31 = _mm256_mul_epu32(_310, _311);   __m128i _prod32 = _mm_mul_epu32(_320, _321);
  __m256i _prod40 = _mm256_mul_epu32(_400, _401);   __m256i _prod41 = _mm256_mul_epu32(_410, _411);   __m128i _prod42 = _mm_mul_epu32(_420, _421);
  
  
  
  __m256i _z0123 = _mm256_add_epi64(_prod00, _prod10);  __m256i _z4567 = _mm256_add_epi64(_prod01, _prod11);  __m128i _z89 = _mm_add_epi64(_prod02, _prod12);
  _z0123 = _mm256_add_epi64(_z0123, _prod20);           _z4567 = _mm256_add_epi64(_z4567, _prod21);           _z89 = _mm_add_epi64(_z89, _prod22);
  _z0123 = _mm256_add_epi64(_z0123, _prod30);           _z4567 = _mm256_add_epi64(_z4567, _prod31);           _z89 = _mm_add_epi64(_z89, _prod32);
  _z0123 = _mm256_add_epi64(_z0123, _prod40);           _z4567 = _mm256_add_epi64(_z4567, _prod41);           _z89 = _mm_add_epi64(_z89, _prod42);

  ///*
  _mm256_storeu_si256((__m256i *)&z[0], _z0123);
  _mm256_storeu_si256((__m256i *)&z[4], _z4567);
  _mm_storeu_si128((__m128i *)&z[8], _z89);

  
  __m256i _500 = _mm256_set_epi64x(_19x3, _38x2, _19x1, _38x0); __m256i _501 = _mm256_set_epi64x(_x3, _x2, _x1, _x0);
  __m256i _prod50 = _mm256_mul_epu32(_500, _501);
  
  /*
  __m256i _z1357 = _mm256_set_epi64x(z[7], z[5], z[3], z[1]);
  __m256i _new = _mm256_add_epi64(_z1357, _prod50);
  z[1] = _new[0];
  z[3] = _new[1];
  z[5] = _new[2];
  z[7] = _new[3];
  //*/
  z[1] += _prod50[0];
  z[3] += _prod50[1];
  z[5] += _prod50[2];
  z[7] += _prod50[3];
  z[9] += _38x4*_x4;
  //*/
  
  /*
  z[0] = _2x9*_x0 + _2x8*_x1    + _2x7*_x2   + _2x6*_x3   + _2x5*_x4;
  z[1] = _2x9*_x1 + _2x8*_2x2   + _2x7*_x3   + _2x6*_2x4  + _x5*_x5    + _38x0*_x0;
  z[2] = _2x9*_x2 + _2x8*_x3    + _2x7*_x4   + _2x6*_x5   + _19x1*_2x0;
  z[3] = _2x9*_x3 + _2x8*_2x4   + _2x7*_x5   + _2x6*_x6   + _38x2*_2x0 + _19x1*_x1;
  z[4] = _2x9*_x4 + _2x8*_x5    + _2x7*_x6   + _19x3*_2x0 + _19x2*_2x1;
  z[5] = _2x9*_x5 + _2x8*_2x6   + _x7*_x7    + _38x4*_2x0 + _19x3*_2x1 + _38x2*_x2;
  z[6] = _2x9*_x6 + _2x8*_x7    + _19x5*_2x0 + _19x4*_2x1 + _19x3*_2x2;
  z[7] = _2x9*_x7 + _2x8*_x8    + _38x6*_2x0 + _19x5*_2x1 + _38x4*_2x2 + _19x3*_x3;
  z[8] = _2x9*_x8 + _19x7*_2x0  + _19x6*_2x1 + _19x5*_2x2 + _19x4*_2x3;
  z[9] = _x9*_x9  + _38x8*_2x0  + _19x7*_2x1 + _38x6*_2x2 + _19x5*_2x3 + _38x4*_x4;
  
  assert(z[0] == _z0123[0]);
  assert(z[1] == _z0123[1] + _38x0*_x0);
  assert(z[2] == _z0123[2]);
  assert(z[3] == _z0123[3] + _19x1*_x1);
  assert(z[4] == _z4567[0]);
  assert(z[5] == _z4567[1] + _38x2*_x2);
  assert(z[6] == _z4567[2]);
  assert(z[7] == _z4567[3] + _19x3*_x3);
  assert(z[8] == _z89[0]);
  assert(z[9] == _z89[1] + _38x4*_x4);
  //*/
  
  // _carry_and_normalize(z);
  _carry_and_normalize(z);
  /* at this point the limb9 is 43bits long, since limb0 had a 38bits carry, then multiplied by 19 becomes 42bits
  *  42 bits + 26bits = 43 bits
  * an additional carry limb9->limb8 of 43-26=17bits
  */
  limb_t _t9 = z[9]; limb_t _t8 = z[8];
  limb_t _r90 = _t9 >> 26;
  z[8] = _t8 + _r90; //z8
  z[9] = _t9 & RADIXMASK26;
  // _fix_degenerate(z);
}
