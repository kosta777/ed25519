#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "radix25_common.h"

/* From the Go language:
*  Add64 returns the sum with carry of x, y and carry: sum = x + y + carry.
*  The carry input must be 0 or 1; otherwise the behavior is undefined.
*  The carryOut output is guaranteed to be 0 or 1.
*/
u_int64_t _addc64(const u_int64_t x, const u_int64_t y, u_int64_t* z, u_int64_t c){

  *z = x + y + c;

  return ((x & y) | ((x | y) & ~(*z))) >> 63;
};

void _carry_and_normalize(limb_t* z){
  // ulli tmp;
  // tmp = z[4];
  // tmp = tmp >>51; // isolate carries
  // z[3]+=tmp; // carry over
  // z[4]&=RADIXMASK51; // zero the carries
  limb_t _t9 = z[9]; limb_t _t8 = z[8];
  limb_t _r90 = _t9 >> 26;
  _t8 = _t8 + _r90; //z8
  _t9 = _t9 & RADIXMASK26;
  
  limb_t _t7 = z[7];
  limb_t _r80 = _t8 >> 25;
  _t7 = _t7 + _r80;
  limb_t _z8 = _t8 & RADIXMASK25; //z8

  limb_t _t6 = z[6];
  limb_t _r70 = _t7 >> 26;
  _t6 = _t6 + _r70;
  limb_t _z7 = _t7 & RADIXMASK26; //z7

  limb_t _t5 = z[5];
  limb_t _r60 = _t6 >> 25;
  _t5 = _t5 + _r60;
  limb_t _z6 = _t6 & RADIXMASK25; //z6

  limb_t _t4 = z[4];
  limb_t _r50 = _t5 >> 26;
  _t4 = _t4 + _r50;
  limb_t _z5 = _t5 & RADIXMASK26; //z5
  
  limb_t _t3 = z[3];
  limb_t _r40 = _t4 >> 25;
  _t3 = _t3 + _r40;
  limb_t _z4 = _t4 & RADIXMASK25; //z4

  limb_t _t2 = z[2];
  limb_t _r30 = _t3 >> 26;
  _t2 = _t2 + _r30;
  limb_t _z3 = _t3 & RADIXMASK26; //z3

  limb_t _t1 = z[1];
  limb_t _r20 = _t2 >> 25;
  _t1 = _t1 + _r20;
  limb_t _z2 = _t2 & RADIXMASK25; //z2

  limb_t _t0 = z[0];
  limb_t _r10 = _t1 >> 26;
  _t0 = _t0 + _r10;
  limb_t _z1 = _t1 & RADIXMASK26; //z1

  limb_t _r00 = _t0 >> 25;
  limb_t _z9 = _t9 + _r00*19;     //z9
  limb_t _z0 = _t0 & RADIXMASK25; //z0

  //store
  z[0] = _z0;
  z[1] = _z1;
  z[2] = _z2;
  z[3] = _z3;
  z[4] = _z4;
  z[5] = _z5;
  z[6] = _z6;
  z[7] = _z7;
  z[8] = _z8;
  z[9] = _z9;
}

void _v_carry_and_normalize(v_limb_t* z){
  // ulli tmp;
  // tmp = z[4];
  // tmp = tmp >>51; // isolate carries
  // z[3]+=tmp; // carry over
  // z[4]&=RADIXMASK51; // zero the carries

  __m256i rm25 = _mm256_set1_epi64x(RADIXMASK25);
  __m256i rm26 = _mm256_set1_epi64x(RADIXMASK26);
  __m256i _19 = _mm256_set1_epi64x(19);

  v_limb_t _t9 = z[9]; v_limb_t _t8 = z[8];
  v_limb_t _r90 = _mm256_srli_epi64(_t9, 26);
  _t8 = _mm256_add_epi64(_t8, _r90); //z8
  _t9 = _mm256_and_si256(_t9, rm26);
  
  v_limb_t _t7 = z[7];
  v_limb_t _r80 = _mm256_srli_epi64(_t8, 25);
  _t7 = _mm256_add_epi64(_t7, _r80);
  v_limb_t _z8 = _mm256_and_si256(_t8, rm25); //z8

  v_limb_t _t6 = z[6];
  v_limb_t _r70 = _mm256_srli_epi64(_t7, 26);
  _t6 = _mm256_add_epi64(_t6, _r70);
  v_limb_t _z7 = _mm256_and_si256(_t7, rm26); //z7

  v_limb_t _t5 = z[5];
  v_limb_t _r60 = _mm256_srli_epi64(_t6, 25);
  _t5 = _mm256_add_epi64(_t5, _r60);
  v_limb_t _z6 = _mm256_and_si256(_t6, rm25); //z6

  v_limb_t _t4 = z[4];
  v_limb_t _r50 = _mm256_srli_epi64(_t5, 26);
  _t4 = _mm256_add_epi64(_t4, _r50);
  v_limb_t _z5 = _mm256_and_si256(_t5, rm26); //z5
  
  v_limb_t _t3 = z[3];
  v_limb_t _r40 = _mm256_srli_epi64(_t4, 25);
  _t3 = _mm256_add_epi64(_t3, _r40);
  v_limb_t _z4 = _mm256_and_si256(_t4, rm25); //z4

  v_limb_t _t2 = z[2];
  v_limb_t _r30 = _mm256_srli_epi64(_t3, 26);
  _t2 = _mm256_add_epi64(_t2, _r30);
  v_limb_t _z3 = _mm256_and_si256(_t3, rm26); //z3

  v_limb_t _t1 = z[1];
  v_limb_t _r20 = _mm256_srli_epi64(_t2, 25);
  _t1 = _mm256_add_epi64(_t1, _r20);
  v_limb_t _z2 = _mm256_and_si256(_t2, rm25); //z2

  v_limb_t _t0 = z[0];
  v_limb_t _r10 = _mm256_srli_epi64(_t1, 26);
  _t0 = _mm256_add_epi64(_t0, _r10);
  v_limb_t _z1 = _mm256_and_si256(_t1, rm26); //z1

  v_limb_t _r00 = _mm256_srli_epi64(_t0, 25);
  v_limb_t _z9 = _mm256_add_epi64(_t9, _mm256_mul_epu32(_r00, _19));     //z9
  v_limb_t _z0 = _mm256_and_si256(_t0, rm25); //z0

  //store
  z[0] = _z0;
  z[1] = _z1;
  z[2] = _z2;
  z[3] = _z3;
  z[4] = _z4;
  z[5] = _z5;
  z[6] = _z6;
  z[7] = _z7;
  z[8] = _z8;
  z[9] = _z9;
}

void _fix_degenerate(limb_t* z){
  /* Reversed the logic for faster escaping
  */
  if (z[0]>= 0x1ffffff &&
      z[1]>= 0x3ffffff &&
      z[2]>= 0x1ffffff &&
      z[3]>= 0x3ffffff &&
      z[4]>= 0x1ffffff &&
      z[5]>= 0x3ffffff &&
      z[6]>= 0x1ffffff &&
      z[7]>= 0x3ffffff &&
      z[8]>= 0x1ffffff &&
      z[9]>= 0x3ffffed)
  {
    z[9] = z[9] - 0x3ffffed;
    z[0] = z[1] = z[2] = z[3] = z[4] = z[5] = z[6] = z[7] = z[8] = 0;
    // z[4] += 18;
  }
}

void _v_fix_degenerate(v_limb_t *z) {
  for (int i = 0 ; i < 4 ; i++) {
    if (z[0][i]>= 0x1ffffff &&
        z[1][i]>= 0x3ffffff &&
        z[2][i]>= 0x1ffffff &&
        z[3][i]>= 0x3ffffff &&
        z[4][i]>= 0x1ffffff &&
        z[5][i]>= 0x3ffffff &&
        z[6][i]>= 0x1ffffff &&
        z[7][i]>= 0x3ffffff &&
        z[8][i]>= 0x1ffffff &&
        z[9][i]>= 0x3ffffed)
    {
      z[9][i] = z[9][i] - 0x3ffffed;
      z[0][i] = z[1][i] = z[2][i] = z[3][i] = z[4][i] = z[5][i] = z[6][i] = z[7][i] = z[8][i] = 0;
      // z[4] += 18;
    }
  }
}
