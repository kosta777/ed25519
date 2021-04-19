#include "radix25_logic.h"

void _v_r25_and(const v_limb_t *x, const v_limb_t *y, v_limb_t *z) {
  for (int i = 0 ; i < 10 ; i++) {
    z[i] = _mm256_and_si256(x[i], y[i]);
  }
}

void _v_r25_xor(const v_limb_t *x, const v_limb_t *y, v_limb_t *z) {
  for (int i = 0 ; i < 10 ; i++) {
    z[i] = _mm256_xor_si256(x[i], y[i]);
  }
}

void _v_r25_not(const v_limb_t *x, v_limb_t *z) {
  // ~x = x xor 0xffffffffffffffff
  // x[i] = 0 -> ~x[i] = x[i] xor 1 = 0 xor 1 = 1
  // x[i] = 1 -> ~x[i] = x[i] xor 1 = 1 xor 1 = 0
  __m256i _one = _mm256_set1_epi64x(0xffffffffffffffff);
  for (int i = 0 ; i < 10 ; i++) {
    z[i] = _mm256_xor_si256(x[i], _one);
  }
}

void _v_r25_mask(const v_limb_t *x, const __m256i m, v_limb_t *z) {
  for (int i = 0 ; i < 10 ; i++) {
    z[i] = _mm256_and_si256(x[i], m);
  }
}
