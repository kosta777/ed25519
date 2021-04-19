#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <immintrin.h>

#include "radix25_addsub.h"

void _v_r25_add(const v_limb_t *x, const v_limb_t *y, v_limb_t *z) {
  for (int i = 0 ; i < 10 ; i++) {
    z[i] = _mm256_add_epi32(x[i], y[i]);
  }
  
  // _v_carry_and_normalize(z);
  // _v_fix_degenerate(z);
}

void _v_r25_sub(const v_limb_t *x, const v_limb_t *y, v_limb_t *z) {
  __m256i _3ffe = _mm256_set1_epi64x(0x3fffffe);
  __m256i _7ffe = _mm256_set1_epi64x(0x7fffffe);
  __m256i _7fda = _mm256_set1_epi64x(0x7ffffda);
  
  for (int i = 0 ; i < 8 ; i += 2) {
    z[i] = _mm256_sub_epi32(_mm256_add_epi32(x[i], _3ffe), y[i]);
    z[i+1] = _mm256_sub_epi32(_mm256_add_epi32(x[i+1], _7ffe), y[i+1]);
  }
  z[8] = _mm256_sub_epi32(_mm256_add_epi32(x[8], _3ffe), y[8]);
  z[9] = _mm256_sub_epi32(_mm256_add_epi32(x[9], _7fda), y[9]);
  
  _v_carry_and_normalize(z);
  _v_fix_degenerate(z);
}

void _v_r25_neg(const v_limb_t *x, v_limb_t *z) {
  __m256i _1fff = _mm256_set1_epi64x(0x1ffffff);
  __m256i _3fff = _mm256_set1_epi64x(0x3ffffff);
  __m256i _3fed = _mm256_set1_epi64x(0x3ffffed);
  
  for (int i = 0 ; i < 8 ; i += 2) {
    z[i] = _mm256_sub_epi32(_1fff, x[i]);
    z[i+1] = _mm256_sub_epi32(_3fff, x[i+1]);
  }
  z[8] = _mm256_sub_epi32(_1fff, x[8]);
  z[9] = _mm256_sub_epi32(_3fed, x[9]);

  _v_carry_and_normalize(z);
}
