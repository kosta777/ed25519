#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <immintrin.h>

#include "radix25_addsub.h"

/* z = x + y */
void _r25_add(const limb_t* x, const limb_t* y, limb_t* z){
  /*
  // add limbs
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
  z[3] = x[3] + y[3];
  z[4] = x[4] + y[4];
  z[5] = x[5] + y[5];
  z[6] = x[6] + y[6];
  z[7] = x[7] + y[7];
  z[8] = x[8] + y[8];
  z[9] = x[9] + y[9];

  // carry over and normalize
  _carry_and_normalize(z);
  // _fix_degenerate(z);
  //*/
  
  ///*
  __m256i _x0123 = _mm256_loadu_si256((__m256i const *)&x[0]);
  __m256i _x4567 = _mm256_loadu_si256((__m256i const *)&x[4]);
  __m128i _x89 = _mm_loadu_si128((__m128i const *)&x[8]);
  __m256i _y0123 = _mm256_loadu_si256((__m256i const *)&y[0]);
  __m256i _y4567 = _mm256_loadu_si256((__m256i const *)&y[4]);
  __m128i _y89 = _mm_loadu_si128((__m128i const *)&y[8]);
  
  __m256i _z0123 = _mm256_add_epi32(_x0123, _y0123);
  __m256i _z4567 = _mm256_add_epi32(_x4567, _y4567);
  __m128i _z89 = _mm_add_epi32(_x89, _y89);

  _mm256_storeu_si256((__m256i *)&z[0], _z0123);
  _mm256_storeu_si256((__m256i *)&z[4], _z4567);
  _mm_storeu_si128((__m128i *)&z[8], _z89);
  
  // _carry_and_normalize(z);
  //*/
}

void _r25_sub(const limb_t* x, const limb_t* y, limb_t* z){
  /* First add 2*q to avoid underflow
  //*/
  
  /*
  z[0] = (x[0] + 0x3fffffe) - y[0];
  z[1] = (x[1] + 0x7fffffe) - y[1];
  z[2] = (x[2] + 0x3fffffe) - y[2];
  z[3] = (x[3] + 0x7fffffe) - y[3];
  z[4] = (x[4] + 0x3fffffe) - y[4];
  z[5] = (x[5] + 0x7fffffe) - y[5];
  z[6] = (x[6] + 0x3fffffe) - y[6];
  z[7] = (x[7] + 0x7fffffe) - y[7];
  z[8] = (x[8] + 0x3fffffe) - y[8];
  z[9] = (x[9] + 0x7ffffda) - y[9];

  _carry_and_normalize(z);
  // _fix_degenerate(z);
  //*/
  
  ///*
  __m256i _3737 = _mm256_set_epi64x(0x7fffffe, 0x3fffffe, 0x7fffffe, 0x3fffffe);
  __m128i _37 = _mm_set_epi64x(0x7ffffda, 0x3fffffe);
  
  __m256i _x0123 = _mm256_loadu_si256((__m256i const *)&x[0]);
  __m256i _x4567 = _mm256_loadu_si256((__m256i const *)&x[4]);
  __m128i _x89 = _mm_loadu_si128((__m128i const *)&x[8]);
  __m256i _y0123 = _mm256_loadu_si256((__m256i const *)&y[0]);
  __m256i _y4567 = _mm256_loadu_si256((__m256i const *)&y[4]);
  __m128i _y89 = _mm_loadu_si128((__m128i const *)&y[8]);
  
  __m256i _tmp0123 = _mm256_add_epi32(_x0123, _3737);
  __m256i _tmp4567 = _mm256_add_epi32(_x4567, _3737);
  __m128i _tmp89 = _mm_add_epi32(_x89, _37);
  
  __m256i _z0123 = _mm256_sub_epi32(_tmp0123, _y0123);
  __m256i _z4567 = _mm256_sub_epi32(_tmp4567, _y4567);
  __m128i _z89 = _mm_sub_epi32(_tmp89, _y89);
  
  _mm256_storeu_si256((__m256i *)&z[0], _z0123);
  _mm256_storeu_si256((__m256i *)&z[4], _z4567);
  _mm_storeu_si128((__m128i *)&z[8], _z89);
  
  _carry_and_normalize(z);
  _fix_degenerate(z);
  //*/
}

void _r25_neg(const limb_t* x, limb_t* z){
  /*
  z[0] = 0x1ffffff - x[0];
  z[1] = 0x3ffffff - x[1];
  z[2] = 0x1ffffff - x[2];
  z[3] = 0x3ffffff - x[3];
  z[4] = 0x1ffffff - x[4];
  z[5] = 0x3ffffff - x[5];
  z[6] = 0x1ffffff - x[6];
  z[7] = 0x3ffffff - x[7];
  z[8] = 0x1ffffff - x[8];
  z[9] = 0x3ffffed - x[9];

  _carry_and_normalize(z);
  // _fix_degenerate(z);
  //*/
  
  ///*
  __m256i _1313 = _mm256_set_epi64x(0x3ffffff, 0x1ffffff, 0x3ffffff, 0x1ffffff);
  __m128i _13 = _mm_set_epi64x(0x3ffffed, 0x1ffffff);
  
  __m256i _x0123 = _mm256_loadu_si256((__m256i const *)&x[0]);
  __m256i _x4567 = _mm256_loadu_si256((__m256i const *)&x[4]);
  __m128i _x89 = _mm_loadu_si128((__m128i const *)&x[8]);
  
  __m256i _z0123 = _mm256_sub_epi32(_1313, _x0123);
  __m256i _z4567 = _mm256_sub_epi32(_1313, _x4567);
  __m128i _z89 = _mm_sub_epi32(_13, _x89);
  
  _mm256_storeu_si256((__m256i *)&z[0], _z0123);
  _mm256_storeu_si256((__m256i *)&z[4], _z4567);
  _mm_storeu_si128((__m128i *)&z[8], _z89);
  
  _carry_and_normalize(z);
  //*/
}
