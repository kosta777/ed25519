#include "scalarmult.h"
#include "../../const.h"
#include "../radix25/radix25.h"
#if R51
#include "../radix/radix51.h"
#endif
#include <stdlib.h>

#include <iostream>
#include <immintrin.h>

// ONLY WORKS WITH RADIX 25
inline void v_edwards(const v_limb_t *v_Px_limbs, const v_limb_t *v_Py_limbs, const v_limb_t *v_Qx_limbs, const v_limb_t *v_Qy_limbs, v_limb_t *_outX_limbs, v_limb_t *_outY_limbs) {
//  In terms of the old edwards funnction:
//  v_Px_limbs == x1;
//  v_Py_limbs == y1;
//  v_Qx_limbs == x2;
//  v_Qy_limbs == y2;

  
  v_limb_t _x1y2[10]; _v_r25_mul(v_Px_limbs, v_Qy_limbs, _x1y2); // x1*y2
  v_limb_t _x2y1[10]; _v_r25_mul(v_Qx_limbs, v_Py_limbs, _x2y1); // x2*y1
  v_limb_t _x3[10]; _v_r25_add(_x1y2, _x2y1, _x3); // x1*y2 + x2*y1
  
  v_limb_t _y1y2[10]; _v_r25_mul(v_Py_limbs, v_Qy_limbs, _y1y2); // y1*y2
  v_limb_t _x1x2[10]; _v_r25_mul(v_Px_limbs, v_Qx_limbs, _x1x2); // x1*x2
  v_limb_t _y3[10]; _v_r25_add(_y1y2, _x1x2, _y3); // y1*y2 + x1*x2
  
  v_limb_t _all[10]; _v_r25_mul(_x1x2, _y1y2, _all); // x1 * x2 * y1 * y2
  v_limb_t _dall[10]; _v_r25_mul(&_d[0], _all, _dall); // CURVE_D * x1 * x2 * y1 * y2

  v_limb_t _1pdall[10]; _v_r25_add(&_one[0], _dall, _1pdall); // 1 + (CURVE_D * x1 * x2 * y1 * y2)
  v_limb_t _1mdall[10]; _v_r25_sub(&_one[0], _dall, _1mdall); // 1 - (CURVE_D * x1 * x2 * y1 * y2)

  _v_r25_div(_x3, _1pdall, _outX_limbs);
  _v_r25_div(_y3, _1mdall, _outY_limbs);
}




/*
* Conditional move
* Sets P<-Q if b != 0, otherwise P remains unchanged
*/
inline void v_CMOV_point(short b[4], v_limb_t *v_Px_limbs, v_limb_t *v_Py_limbs, const v_limb_t *v_Qx_limbs, const v_limb_t *v_Qy_limbs) {
  __m256i _one = _mm256_set1_epi64x(0xffffffffffffffff);
  __m256i _m = _mm256_set_epi64x((int64_t)(0 - b[3]), (int64_t)(0 - b[2]), (int64_t)(0 - b[1]), (int64_t)(0 - b[0]));
  __m256i _notm = _mm256_xor_si256(_m, _one);

  v_limb_t P0m[10]; _v_r25_mask(v_Px_limbs, _notm, &P0m[0]);
  v_limb_t Q0m[10]; _v_r25_mask(v_Qx_limbs, _m, &Q0m[0]);
  v_limb_t P1m[10]; _v_r25_mask(v_Py_limbs, _notm, &P1m[0]);
  v_limb_t Q1m[10]; _v_r25_mask(v_Qy_limbs, _m, &Q1m[0]);
  
  _v_r25_xor(P0m, Q0m, v_Px_limbs);
  _v_r25_xor(P1m, Q1m, v_Py_limbs);
}




void v_scalar_mult(const radix_t P[4][2], const radix_t e[4], radix_t *P00, radix_t *P01, radix_t *P02, radix_t *P03) {
#ifndef R25
  std::cerr << "Vectorised computations only supported with radix25.5. Running regular computations instead." << std::endl;
  scalarmult(&P[0][0], e[0], P00);
  scalarmult(&P[1][0], e[1], P01);
  scalarmult(&P[2][0], e[2], P02);
  scalarmult(&P[3][0], e[3], P03);
#endif

  // In terms of the old scalarmult:
  // P0 (in the old one) is the output value (so P0x here, x = 0..3)
  // P1 (in the old one) is the input value (so P[x] here, x = 0..3)
  
  v_limb_t P0[2][10];
  for (int i = 0 ; i < 9 ; i++) {
    P0[0][i] = _mm256_setzero_si256();
    P0[1][i] = _mm256_setzero_si256();
  }
  P0[0][9] = _mm256_setzero_si256();
  P0[1][9] = _mm256_set1_epi64x(1);
  
  v_limb_t P1[2][10];
  for (int i = 0 ; i < 10 ; i++) {
    P1[0][i] = _mm256_set_epi64x(P[3][0].limbs[i], P[2][0].limbs[i], P[1][0].limbs[i], P[0][0].limbs[i]);
    P1[1][i] = _mm256_set_epi64x(P[3][1].limbs[i], P[2][1].limbs[i], P[1][1].limbs[i], P[0][1].limbs[i]);
  }

  short b[4] = {0, 0, 0, 0};
  short ki1[4] = {0, 0, 0, 0};
  short ki[4] = {0, 0, 0, 0};
  v_limb_t P2[2][10];
  for (int i = 0 ; i < 10 ; i++) {
    P2[0][i] = _mm256_setzero_si256();
    P2[1][i] = _mm256_setzero_si256();
  }

  for (int i = 255 ; i >= 0 ; i--) {
    for (int j = 0 ; j < 4 ; j++) {
      ki[j] = get_bit(e[j], i);
      b[j] = ki[j] ^ ki1[j];
    }

    v_edwards(&P0[0][0], &P0[1][0], &P1[0][0], &P1[1][0], &P2[0][0], &P2[1][0]); // P2 <- P0 + P1
    v_CMOV_point(b, &P0[0][0], &P0[1][0], &P1[0][0], &P1[1][0]); // P0 <- CMOV(b, P0, P1)
    v_edwards(&P0[0][0], &P0[1][0], &P0[0][0], &P0[1][0], &P0[0][0], &P0[1][0]); // P0 <- 2 * P0
    for (int j = 0 ; j < 10 ; j++) {
      P1[0][j] = P2[0][j]; // P1 <- P2
      P1[1][j] = P2[1][j]; // P1 <- P2
    }

    for (int j = 0 ; j < 4 ; j++) {
      ki1[j] = ki[j];
    }
  }

  v_CMOV_point(ki, P0[0], P0[1], P1[0], P1[1]); // P0 <- CMOV(k_0, P0, P1)
  // Writing to output
  for (int i = 0; i < 10; i++)
  {
    P00[0].limbs[i] = P0[0][i][0];
    P00[1].limbs[i] = P0[1][i][0];
    P01[0].limbs[i] = P0[0][i][1];
    P01[1].limbs[i] = P0[1][i][1];
    P02[0].limbs[i] = P0[0][i][2];
    P02[1].limbs[i] = P0[1][i][2];
    P03[0].limbs[i] = P0[0][i][3];
    P03[1].limbs[i] = P0[1][i][3];
  }
}
