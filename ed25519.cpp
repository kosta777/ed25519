#include <openssl/sha.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "ed25519.h"
// #include "scalarmult/scalarmult.h"
// #include "radix/radix51.h"
// #include "radix/radix51_expmod.h"
#include "const.h"
#include "util.h"

unsigned int bit_at(unsigned char *m, int i) {
  int idx = i / 8;
  int bit = i % 8;
  return (m[idx] >> bit) & 1;
}


radix_t x_recover(radix_t y) {
  radix_t xx = y * (y-ONE);
  xx = xx * CURVE_D_INV;
  xx = xx * y * (y+ONE);
  radix_t x = expmod1(xx.limbs);
  
  if ((x*x - xx) != ZERO) {
    x = x * CURVE_I;
  }
  if ((x.limbs[4] & 1) != 0) {
    x = CURVE_Q - x;
  }

  return x;
}


void encode_int(radix_t y, char *dst) {
  short bits[256]; // CURVE_B = 256
  ulli yulli[4];
  radixtoulli(&y, &yulli[0]);
  for (int i = 0 ; i < 256 ; i++) {
    // bits[i] = (y >> i).limbs[4] & 1; // Only LSB
    bits[i] = (yulli[3-i/64] >> i % 64) & 1;
    // bits[i] = get_bit(y, i);
  }
  
  int top = 32; // CURVE_B = 256, 256 / 8 = 32
  for (int i = 0 ; i < top ; i++) {
    int sum = 0;
    for (int j = 0 ; j < 8 ; j++) {
      sum += ((int)bits[i*8 + j] << j);
    }
    snprintf(dst+(i*2), 8, "%02x", sum);
  }
}

inline void encode_point(radix_t *P, char *dst) {
  radix_t x = P[0];
  radix_t y = P[1];
  ulli yulli[4];
  radixtoulli(&y, &yulli[0]);
  short bits[256]; // CURVE_B = 256
  for (int i = 0 ; i < 255 ; i++) {
    // short foo = (y >> i).limbs[4] & 1; // Only LSB
    bits[i] = (yulli[3-i/64] >> i % 64) & 1;
    // bits[i] = get_bit(y, i);
  }
#if R25
  bits[255] = x.limbs[9] & 1; // Only LSB
#else
  bits[255] = x.limbs[4] & 1; // Only LSB
#endif
  int top = 32; // CURVE_B / 8 = 32
  for (int i = 0 ; i < top ; i++) {
    int sum = 0;
    for (int j = 0 ; j < 8 ; j++) {
      sum += ((int)bits[i*8 + j] << j);
    }
    snprintf(dst+(i*2), 4, "%02x", sum);
  }
}

inline void H(const unsigned char *m, int n, unsigned char *dst) {
  SHA512(m, n, (unsigned char *)dst);
}


void publickey(char *sk, char *pk) {
  unsigned char h[SHA512_DIGEST_LENGTH];
  radix_t P0[2];
  
  H((const unsigned char *)sk, 32, h);
  ulli _a[5];
  u_int64_t *hint = (u_int64_t *)h;
  _a[4] = (hint)[0] & 0xfffffffffffffff8; // mask out lower 3 bits
  _a[3] = (hint)[1];
  _a[2] = (hint)[2];
  _a[1] = (hint)[3] | ulli_2p254[1]; // + TWOp254
  _a[1] &= 0x7fffffffffffffff;
  _a[0] = 0;
  radix_t a(_a);
  
  P0[0] = ZERO;
  P0[1] = ZERO;
  
  radix_t *A = scalar_mult((radix_t *)CURVE_B2, a, &P0[0]);
  encode_point(A, pk);
}

// Computes n public keys
void publickeys(int n, char sk[][32], char pk[][64]) {
#if SSE && R25
  v_publickeys(n, sk, pk);
  return;
#endif
  
  unsigned char h[SHA512_DIGEST_LENGTH];
  radix_t P0[2];
  for (int k = 0 ; k < n ; k++) {
    H((const unsigned char *)(sk[k]), 32, h);
    ulli _a[5];
    u_int64_t *hint = (u_int64_t *)h;
    _a[4] = hint[0] & 0xfffffffffffffff8; // mask out lower 3 bits
    _a[3] = hint[1];
    _a[2] = hint[2];
    _a[1] = hint[3] | ulli_2p254[1]; // + TWOp254
    _a[1] &= 0x7fffffffffffffff;
    _a[0] = 0;
    radix_t a(_a);

#if R25
    P0[0] = R25_ZERO;
    P0[1] = R25_ZERO;
#else
    P0[0] = ZERO;
    P0[1] = ZERO;
#endif
    
#if BASE
    radix_t *A = scalar_mult((radix_t *)CURVE_B2, a, &P0[0]);
#elif R25
    radix_t *A = scalar_mult_branchless((radix_t *)R25_CURVE_B2, a, &P0[0]);
#else
    radix_t *A = scalar_mult_branchless((radix_t *)CURVE_B2, a, &P0[0]);
#endif
    encode_point(A, pk[k]);
  }
}

#if SSE && R25
void v_publickeys(int n, char sk[][32], char pk[][64]) {
  unsigned char h[SHA512_DIGEST_LENGTH];
  radix_t P0[2];
  int peel = n % 4; // How many cannot be vectorised
  for (int k = 0 ; k < peel ; k++) {
    H((const unsigned char *)(sk[k]), 32, h);
    ulli _a[5];
    u_int64_t *hint = (u_int64_t *)h;
    _a[4] = hint[0] & 0xfffffffffffffff8; // mask out lower 3 bits
    _a[3] = hint[1];
    _a[2] = hint[2];
    _a[1] = hint[3] | ulli_2p254[1]; // + TWOp254
    _a[1] &= 0x7fffffffffffffff;
    _a[0] = 0;
    radix_t a(_a);
    
    P0[0] = R25_ZERO;
    P0[1] = R25_ZERO;
    
    radix_t *A = scalar_mult_branchless((radix_t *)R25_CURVE_B2, a, &P0[0]);
    encode_point(A, pk[k]);
  }
    
  unsigned char hv[4][SHA512_DIGEST_LENGTH];
  radix_t P0v[4][2];
  radix_t B2v[4][2];
  for (int i = 0 ; i < 4 ; i++) {
    B2v[i][0] = CURVE_B2[0];
    B2v[i][1] = CURVE_B2[1];
  }
    
  for (int k = peel ; k < n ; k += 4) {
    for (int i = 0 ; i < 4 ; i++) {
      H((const unsigned char *)(sk[k+i]), 32, &hv[i][0]);
    }
    ulli _a[4][5];
    u_int64_t *hint[4];
    for (int i = 0 ; i < 4 ; i++) {
      hint[i] = (u_int64_t *)hv[i];
    }
    for (int i = 0 ; i < 4 ; i++) {
      _a[i][4] = hint[i][0] & 0xfffffffffffffff8; // mask out lower 3 bits
      _a[i][3] = hint[i][1];
      _a[i][2] = hint[i][2];
      _a[i][1] = hint[i][3] | ulli_2p254[1]; // + TWOp254
      _a[i][1] &= 0x7fffffffffffffff;
      _a[i][0] = 0;
    }
      
    radix_t a[4];
    for (int i = 0 ; i < 4 ; i++) {
      a[i] = radix_t(_a[i]);
    }
      
    v_scalar_mult(B2v, a, &P0v[0][0], &P0v[1][0], &P0v[2][0], &P0v[3][0]);
    
    for (int i = 0 ; i < 4 ; i++) {
      encode_point(&P0v[i][0], pk[k+i]);
    }
  }
}
#endif




/*
r51_t H_int(const unsigned char *m) {
  char *hash = (char *)malloc(SHA512_DIGEST_LENGTH);
  H(m, 32, hash);
  Radix51 acc = ZERO;
  Radix51 twopi = ONE;
  for (int i = 0 ; i < 512 ; i++) { // 2 * CURVE_B
    if (bit_at(hash, i)) {
      acc = acc + twopi;
    }
    twopi = twopi * TWO;
  }
  // TODO: this fails. acc is not correct (because it is mod Q and not on 512 bits), but hash is
  free(hash);
  return acc;
}


void signature(char *m, char *sk, char *pk, char *sig) {
  char *h = (char *)malloc(SHA512_DIGEST_LENGTH);
  H((const unsigned char *)sk, 32, h);
  Radix51 a = TWOp254;
  for (int i = 3 ; i < 254 ; i++) {
    if (bit_at(h, i)) {
      a = a + (ONE << i);
    }
  }
//  std::cerr << "a = 0x" << str_rep(&a) << std::endl;
  
  char *bits = (char *)calloc(32 + strlen(m), sizeof(char)); // h[32..63], m
  memcpy(bits, &h[32], 32);
  memcpy(bits+32, m, strlen(m));
//  std::cerr << "bits = " << hexlify(bits, 32 + strlen(m)) << std::endl;

  Radix51 r = H_int((const unsigned char *)bits); // TODO: find a way to get the hash as a number (512 bits, not moduloed)
//  std::cerr << "r % q = 0x" << str_rep(&r) << std::endl;
  free(bits);

  Radix51 *R = scalar_mult((r51_t *)CURVE_B2, r);
  char *Rpkm = (char *)calloc(128 + strlen(m), sizeof(char)); // encoded R + pk + m = 64 + 64 + strlen(m)
  encode_point(R, Rpkm);
//  std::cerr << "R         = 0x" << Rpkm << std::endl;
  memcpy(Rpkm+64, pk, 64);
  memcpy(Rpkm+128, m, strlen(m));
  
  Radix51 S = r + (H_int((const unsigned char *)Rpkm) * a);
  free(Rpkm);
  while (S > CURVE_L || S == CURVE_L) {
    S = S - CURVE_L;
  }
  char *S_enc = (char *)calloc(64, sizeof(char));
  encode_int(S, S_enc);
//  std::cerr << "S         = 0x"; for (int i = 0 ; i < 64 ; i++) std::cerr << " "; std::cerr << S_enc << std::endl;
  free(S_enc);
  
  encode_point(R, sig);
  free(R);
  encode_int(S, sig+64);
//  std::cerr << "Signature = 0x" << sig << std::endl;
  
}


bool is_on_curve(Radix51 *P) {
  Radix51 x = P[0];
  Radix51 y = P[1];
  return (-x*x + y*y - ONE - CURVE_D*x*x*y*y) == ZERO;
}

Radix51 decode_int(char *s) {
  Radix51 res = ZERO;
  for (int i = 0 ; i < 256 ; i++) { // CURVE_B
    // res = res + ((1 << i) * bit_at(s, i));
    res = res + ((ONE << i) * (bit_at(s, i) ? ONE : ZERO));
  }
  return res;
}

Radix51 *decode_point(char *s) {
  Radix51 y = ZERO;
  Radix51 twopi = ONE;
  for (int i = 0 ; i < 255 ; i++) { // CURVE_B
    // y = y + ((1 << i) * bit_at(s, i));
    if (bit_at(s, i)) {
      y = y + twopi;
    }
    twopi = twopi * TWO;
  }
  std::cerr << str_rep(&y) << std::endl;
  Radix51 x = x_recover(y);
  if ((x.limbs[4] & 1) != bit_at(s, 255)) {
    x = CURVE_Q - x;
  }
  Radix51 *P = (Radix51 *)malloc(10 * sizeof(ulli));
  P[0] = x;
  P[1] = y;
  if (!is_on_curve(P)) {
//    assert(0 == 1); // I know, this could be more elegant...
  }
  return P;
}

int verify(char *sig, char *m, char *pk){
  int bo8 = 32; // CURVE_B / 8
  int bo4 = 64; // CURVE_B / 4
//  if (strlen(sig) != bo4) {
//    std::cerr << "Wrong signature length (expected 64, got " << strlen(sig) << " instead)" << std::endl;
//    return -1;
//  }
//  if (strlen(pk) != bo8) {
//    std::cerr << "Wrong public key length (expected 32, got " << strlen(pk) << " instead)" << std::endl;
//    return -2;
//  }
  char to_decode_1[bo8];
  for (int i = 0 ; i < bo8 ; i++) {
    to_decode_1[i] = sig[i];
  }
//  std::cerr << "to_decode_1 = " << hexlify(to_decode_1, bo8) << std::endl;
  Radix51 *R = decode_point(to_decode_1);
  Radix51 *A = decode_point(pk);
  
  char to_decode_2[bo4 - bo8];
  for (int i = bo8 ; i < bo4 ; i++) {
    to_decode_2[i - bo8] = sig[i];
  }
//  std::cerr << "to_decode_2 = " << hexlify(to_decode_2, bo4-bo8) << std::endl;
  Radix51 S = decode_int(to_decode_2);
  
  char *Rpkm = (char *)calloc(128 + strlen(m), sizeof(char)); // encoded R + pk + m = 64 + 64 + strlen(m)
  encode_point(R, Rpkm);
  memcpy(Rpkm+64, pk, 64);
  memcpy(Rpkm+128, m, strlen(m));
  Radix51 h = H_int((const unsigned char *)Rpkm); // TODO: again, need 512 bits and no mod

  limb_t empty[] = {0, 0, 0, 0, 0};
  r51_t tempX(empty), tempY(empty);
  edwards(R, scalar_mult(A, h), &tempX, &tempY);
  r51_t *smul_result = scalar_mult((r51_t *)CURVE_B2, S);
  
  int retval = 0;
  if (!(smul_result[0] == tempX && smul_result[1] == tempY)) {
    // std::cerr << "Incorrect parameters" << std::endl;
    retval = -3;
  }
  
  free(R);
  free(A);
  free(Rpkm);
  free(smul_result);
  return retval;
}
*/
