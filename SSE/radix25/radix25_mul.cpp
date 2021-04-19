#include "radix25_mul.h"

void _r25_mul(const limb_t* x, const limb_t* y, limb_t* z){
  /*
   // Somehow this is about the same time
  limb_t _x0 = x[0]; limb_t _19x0 = 19*_x0; limb_t _38x0 = 38*_x0;
  limb_t _x1 = x[1]; limb_t _19x1 = 19*_x1;
  limb_t _x2 = x[2]; limb_t _2x2 = 2*_x2;   limb_t _19x2 = 19*_x2; limb_t _38x2 = 38*_x2;
  limb_t _x3 = x[3]; limb_t _19x3 = 19*_x3;
  limb_t _x4 = x[4]; limb_t _2x4 = 2*_x4;   limb_t _19x4 = 19*_x4; limb_t _38x4 = 38*_x4;
  limb_t _x5 = x[5]; limb_t _19x5 = 19*_x5;
  limb_t _x6 = x[6]; limb_t _2x6 = 2*_x6;   limb_t _19x6 = 19*_x6; limb_t _38x6 = 38*_x6;
  limb_t _x7 = x[7]; limb_t _19x7 = 19*_x7;
  limb_t _x8 = x[8]; limb_t _2x8 = 2*_x8;   limb_t _38x8 = 38*_x8;
  limb_t _x9 = x[9];
  limb_t _y0 = y[0];
  limb_t _y1 = y[1];
  limb_t _y2 = y[2];
  limb_t _y3 = y[3];
  limb_t _y4 = y[4];
  limb_t _y5 = y[5];
  limb_t _y6 = y[6];
  limb_t _y7 = y[7];
  limb_t _y8 = y[8];
  limb_t _y9 = y[9];
  //*/
 
  ///*
    limb_t _x0 = x[0]; // limb_t _19x0 = 19*_x0; limb_t _38x0 = 38*_x0;
    limb_t _x1 = x[1]; // limb_t _19x1 = 19*_x1;
    limb_t _x2 = x[2]; // limb_t _2x2 = 2*_x2;   limb_t _19x2 = 19*_x2; limb_t _38x2 = 38*_x2;
    limb_t _x3 = x[3]; // limb_t _19x3 = 19*_x3;
    limb_t _x4 = x[4]; // limb_t _2x4 = 2*_x4;   limb_t _19x4 = 19*_x4; limb_t _38x4 = 38*_x4;
    limb_t _x5 = x[5]; // limb_t _19x5 = 19*_x5;
    limb_t _x6 = x[6]; // limb_t _2x6 = 2*_x6;   limb_t _19x6 = 19*_x6; limb_t _38x6 = 38*_x6;
    limb_t _x7 = x[7]; // limb_t _19x7 = 19*_x7;
    limb_t _x8 = x[8]; limb_t _2x8 = 2*_x8;   limb_t _38x8 = 38*_x8;
    limb_t _x9 = x[9];
    limb_t _y0 = y[0];
    limb_t _y1 = y[1];
    limb_t _y2 = y[2];
    limb_t _y3 = y[3];
    limb_t _y4 = y[4];
    limb_t _y5 = y[5];
    limb_t _y6 = y[6];
    limb_t _y7 = y[7];
    limb_t _y8 = y[8];
    limb_t _y9 = y[9];

  __m256i _jj2j = _mm256_set_epi64x(19, 2, 19, 19);
  __m256i _2j2j = _mm256_set_epi64x(19, 2, 19, 2);
  __m256i _jjj = _mm256_set_epi64x(19, 19, 19, 38);
  __m256i _j = _mm256_set_epi64x(38, 38, 38, 1);

  __m256i _x0123 = _mm256_loadu_si256((__m256i const *)&x[0]);
  __m256i _x4567 = _mm256_loadu_si256((__m256i const *)&x[4]);
  __m256i _x0246 = _mm256_set_epi64x(x[6], x[4], x[2], x[0]);

  __m256i _prod00 = _mm256_mul_epu32(_x0123, _jj2j);
  __m256i _prod01 = _mm256_mul_epu32(_x4567, _2j2j);
  __m256i _prod1 = _mm256_mul_epu32(_x0246, _jjj);
  __m256i _prod2 = _mm256_mul_epu32(_x0246, _j);

  limb_t _19x0 = _prod00[0]; limb_t _38x0 = _prod1[0];
  limb_t _19x1 = _prod00[1];
  limb_t _2x2 = _prod00[2];   limb_t _19x2 = _prod1[1]; limb_t _38x2 = _prod2[1];
  limb_t _19x3 = _prod00[3];
  limb_t _2x4 = _prod01[0];   limb_t _19x4 = _prod1[2]; limb_t _38x4 = _prod2[2];
  limb_t _19x5 = _prod01[1];
  limb_t _2x6 = _prod01[2];   limb_t _19x6 = _prod1[3]; limb_t _38x6 = _prod2[3];
  limb_t _19x7 = _prod01[3];
  //*/
  
  // lo0
  __m256i _0000 = _mm256_set_epi64x(_x9, _x9, _x9, _x9);         __m256i _0001 = _mm256_set_epi64x(_y6, _y7, _y8, _y9);
  __m256i _0010 = _0000;                                         __m256i _0011 = _mm256_set_epi64x(_y2, _y3, _y4, _y5);
  __m128i _0020 = _mm_set_epi64x(_x9, _x9);                      __m128i _0021 = _mm_set_epi64x(_y0, _y1);

  __m256i _0100 = _mm256_set_epi64x(_x8, _2x8, _x8, _38x8);      __m256i _0101 = _mm256_set_epi64x(_y7, _y8, _y9, _y0);
  __m256i _0110 = _mm256_set_epi64x(_x8, _2x8, _x8, _2x8);       __m256i _0111 = _mm256_set_epi64x(_y3, _y4, _y5, _y6);
  __m128i _0120 = _mm_set_epi64x(_x8, _2x8);                     __m128i _0121 = _mm_set_epi64x(_y1, _y2);
  
  __m256i _0200 = _mm256_set_epi64x(_x7, _x7, _19x7, _19x7);     __m256i _0201 = _mm256_set_epi64x(_y8, _y9, _y0, _y1);
  __m256i _0210 = _mm256_set_epi64x(_x7, _x7, _x7, _x7);         __m256i _0211 = _mm256_set_epi64x(_y4, _y5, _y6, _y7);
  __m128i _0220 = _mm_set_epi64x(_x7, _x7);                      __m128i _0221 = _mm_set_epi64x(_y2, _y3);
  
  __m256i _0300 = _mm256_set_epi64x(_x6, _38x6, _19x6, _38x6);   __m256i _0301 = _mm256_set_epi64x(_y9, _y0, _y1, _y2);
  __m256i _0310 = _mm256_set_epi64x(_x6, _2x6, _x6, _2x6);       __m256i _0311 = _mm256_set_epi64x(_y5, _y6, _y7, _y8);
  __m128i _0320 = _mm_set_epi64x(_x6, _2x6);                     __m128i _0321 = _mm_set_epi64x(_y3, _y4);
  
  __m256i _0400 = _mm256_set_epi64x(_19x5, _19x5, _19x5, _19x5); __m256i _0401 = _mm256_set_epi64x(_y0, _y1, _y2, _y3);
  __m256i _0410 = _mm256_set_epi64x(_x5, _x5, _x5, _x5);         __m256i _0411 = _mm256_set_epi64x(_y6, _y7, _y8, _y9);
  __m128i _0420 = _mm_set_epi64x(_x5, _x5);                      __m128i _0421 = _mm_set_epi64x(_y4, _y5);
  
  
  __m256i _a00 = _mm256_mul_epu32(_0000, _0001);  __m256i _a01 = _mm256_mul_epu32(_0010, _0011);  __m128i _a02 = _mm_mul_epu32(_0020, _0021);
  __m256i _a10 = _mm256_mul_epu32(_0100, _0101);  __m256i _a11 = _mm256_mul_epu32(_0110, _0111);  __m128i _a12 = _mm_mul_epu32(_0120, _0121);
  __m256i _a20 = _mm256_mul_epu32(_0200, _0201);  __m256i _a21 = _mm256_mul_epu32(_0210, _0211);  __m128i _a22 = _mm_mul_epu32(_0220, _0221);
  __m256i _a30 = _mm256_mul_epu32(_0300, _0301);  __m256i _a31 = _mm256_mul_epu32(_0310, _0311);  __m128i _a32 = _mm_mul_epu32(_0320, _0321);
  __m256i _a40 = _mm256_mul_epu32(_0400, _0401);  __m256i _a41 = _mm256_mul_epu32(_0410, _0411);  __m128i _a42 = _mm_mul_epu32(_0420, _0421);
  
  __m256i _z9876lo0 = _mm256_add_epi64(_a00, _a10);   __m256i _z5432lo0 = _mm256_add_epi64(_a01, _a11);   __m128i _z10lo0 = _mm_add_epi64(_a02, _a12);
  _z9876lo0 = _mm256_add_epi64(_z9876lo0, _a20);      _z5432lo0 = _mm256_add_epi64(_z5432lo0, _a21);      _z10lo0 = _mm_add_epi64(_z10lo0, _a22);
  _z9876lo0 = _mm256_add_epi64(_z9876lo0, _a30);      _z5432lo0 = _mm256_add_epi64(_z5432lo0, _a31);      _z10lo0 = _mm_add_epi64(_z10lo0, _a32);
  _z9876lo0 = _mm256_add_epi64(_z9876lo0, _a40);      _z5432lo0 = _mm256_add_epi64(_z5432lo0, _a41);      _z10lo0 = _mm_add_epi64(_z10lo0, _a42);
  
  u_int64_t _z9876lo0arr[4], _z5432lo0arr[4], _z10lo0arr[2];
  _mm256_storeu_si256((__m256i *)&_z9876lo0arr[0], _z9876lo0);
  _mm256_storeu_si256((__m256i *)&_z5432lo0arr[0], _z5432lo0);
  _mm_store_si128((__m128i *)&_z10lo0arr[0], _z10lo0);
  
  // lo1
  __m256i _1000 = _mm256_set_epi64x(_19x4, _38x4, _19x4, _38x4);  __m256i _1001 = _mm256_set_epi64x(_y1, _y2, _y3, _y4);
  __m256i _1010 = _mm256_set_epi64x(_x4, _2x4, _x4, _38x4);       __m256i _1011 = _mm256_set_epi64x(_y7, _y8, _y9, _y0);
  __m128i _1020 = _mm_set_epi64x(_x4, _2x4);                      __m128i _1021 = _mm_set_epi64x(_y5, _y6);

  __m256i _1100 = _mm256_set_epi64x(_19x3, _19x3, _19x3, _19x3);  __m256i _1101 = _mm256_set_epi64x(_y2, _y3, _y4, _y5);
  __m256i _1110 = _mm256_set_epi64x(_x3, _x3, _19x3, _19x3);      __m256i _1111 = _mm256_set_epi64x(_y8, _y9, _y0, _y1);
  __m128i _1120 = _mm_set_epi64x(_x3, _x3);                       __m128i _1121 = _mm_set_epi64x(_y6, _y7);
  
  __m256i _1200 = _mm256_set_epi64x(_19x2, _38x2, _19x2, _38x2);  __m256i _1201 = _mm256_set_epi64x(_y3, _y4, _y5, _y6);
  __m256i _1210 = _mm256_set_epi64x(_x2, _38x2, _19x2, _38x2);    __m256i _1211 = _mm256_set_epi64x(_y9, _y0, _y1, _y2);
  __m128i _1220 = _mm_set_epi64x(_x2, _2x2);                      __m128i _1221 = _mm_set_epi64x(_y7, _y8);
  
  __m256i _1300 = _mm256_set_epi64x(_19x1, _19x1, _19x1, _19x1);  __m256i _1301 = _mm256_set_epi64x(_y4, _y5, _y6, _y7);
  __m256i _1310 = _mm256_set_epi64x(_19x1, _19x1, _19x1, _19x1);  __m256i _1311 = _mm256_set_epi64x(_y0, _y1, _y2, _y3);
  __m128i _1320 = _mm_set_epi64x(_x1, _x1);                       __m128i _1321 = _mm_set_epi64x(_y8, _y9);
  
  __m256i _1400 = _mm256_set_epi64x(_19x0, _38x0, _19x0, _38x0);  __m256i _1401 = _mm256_set_epi64x(_y5, _y6, _y7, _y8);
  __m256i _1410 = _mm256_set_epi64x(_19x0, _38x0, _19x0, _38x0);  __m256i _1411 = _mm256_set_epi64x(_y1, _y2, _y3, _y4);
  __m128i _1420 = _mm_set_epi64x(_x0, _38x0);                     __m128i _1421 = _mm_set_epi64x(_y9, _y0);
  
  
  __m256i _b00 = _mm256_mul_epu32(_1000, _1001);  __m256i _b01 = _mm256_mul_epu32(_1010, _1011);  __m128i _b02 = _mm_mul_epu32(_1020, _1021);
  __m256i _b10 = _mm256_mul_epu32(_1100, _1101);  __m256i _b11 = _mm256_mul_epu32(_1110, _1111);  __m128i _b12 = _mm_mul_epu32(_1120, _1121);
  __m256i _b20 = _mm256_mul_epu32(_1200, _1201);  __m256i _b21 = _mm256_mul_epu32(_1210, _1211);  __m128i _b22 = _mm_mul_epu32(_1220, _1221);
  __m256i _b30 = _mm256_mul_epu32(_1300, _1301);  __m256i _b31 = _mm256_mul_epu32(_1310, _1311);  __m128i _b32 = _mm_mul_epu32(_1320, _1321);
  __m256i _b40 = _mm256_mul_epu32(_1400, _1401);  __m256i _b41 = _mm256_mul_epu32(_1410, _1411);  __m128i _b42 = _mm_mul_epu32(_1420, _1421);
  
  __m256i _z9876lo1 = _mm256_add_epi64(_b00, _b10);   __m256i _z5432lo1 = _mm256_add_epi64(_b01, _b11);   __m128i _z10lo1 = _mm_add_epi64(_b02, _b12);
  _z9876lo1 = _mm256_add_epi64(_z9876lo1, _b20);      _z5432lo1 = _mm256_add_epi64(_z5432lo1, _b21);      _z10lo1 = _mm_add_epi64(_z10lo1, _b22);
  _z9876lo1 = _mm256_add_epi64(_z9876lo1, _b30);      _z5432lo1 = _mm256_add_epi64(_z5432lo1, _b31);      _z10lo1 = _mm_add_epi64(_z10lo1, _b32);
  _z9876lo1 = _mm256_add_epi64(_z9876lo1, _b40);      _z5432lo1 = _mm256_add_epi64(_z5432lo1, _b41);      _z10lo1 = _mm_add_epi64(_z10lo1, _b42);
  
  u_int64_t _z9876lo1arr[4], _z5432lo1arr[4], _z10lo1arr[2];
  _mm256_storeu_si256((__m256i *)&_z9876lo1arr[0], _z9876lo1);
  _mm256_storeu_si256((__m256i *)&_z5432lo1arr[0], _z5432lo1);
  _mm_store_si128((__m128i *)&_z10lo1arr[0], _z10lo1);

  /* Each multiplication produces 58bits numbers.
  *  Each limb is a summation of 10 products.
  *  Split the sum in two to accumulate in two 64bits parts without overflow
  *  Then sum together the two halves and carry to the next limb
  */
  /*
    // first limb
   u_int64_t _z9lo0 = _x9*_y9 + _38x8*_y0   + _19x7*_y1 + _38x6*_y2 + _19x5*_y3;
   u_int64_t _z9lo1 = _38x4*_y4 + _19x3*_y5 + _38x2*_y6 + _19x1*_y7 + _38x0*_y8;
//   assert(_z9lo0 == _z9876lo0arr[0]);
//  assert(_z9lo1 == _z9876lo1arr[0]);
    // second
    u_int64_t _z8lo0 = _x9*_y8 + _x8*_y9     + _19x7*_y0 + _19x6*_y1 + _19x5*_y2;
    u_int64_t _z8lo1 = _19x4*_y3 + _19x3*_y4 + _19x2*_y5 + _19x1*_y6 + _19x0*_y7;
//    assert(_z8lo0 == _z9876lo0arr[1]);
//    assert(_z8lo1 == _z9876lo1arr[1]);
    // third
    u_int64_t _z7lo0 = _x9*_y7 + _2x8*_y8    + _x7*_y9   + _38x6*_y0 + _19x5*_y1;
    u_int64_t _z7lo1 = _38x4*_y2 + _19x3*_y3 + _38x2*_y4 + _19x1*_y5 + _38x0*_y6;
//    assert(_z7lo0 == _z9876lo0arr[2]);
//    assert(_z7lo1 == _z9876lo1arr[2]);
    // fourth
    u_int64_t _z6lo0 = _x9*_y6 + _x8*_y7     + _x7*_y8   + _x6*_y9   + _19x5*_y0;
    u_int64_t _z6lo1 = _19x4*_y1 + _19x3*_y2 + _19x2*_y3 + _19x1*_y4 + _19x0*_y5;
//    assert(_z6lo0 == _z9876lo0arr[3]);
//    assert(_z6lo1 == _z9876lo1arr[3]);
    // fifth
    u_int64_t _z5lo0 = _x9*_y5 + _2x8*_y6    + _x7*_y7   + _2x6*_y8  + _x5*_y9  ;
    u_int64_t _z5lo1 = _38x4*_y0 + _19x3*_y1 + _38x2*_y2 + _19x1*_y3 + _38x0*_y4;
//    assert(_z5lo0 == _z5432lo0arr[0]);
//    assert(_z5lo1 == _z5432lo1arr[0]);
    // sixth
    u_int64_t _z4lo0 = _x9*_y4 + _x8*_y5     + _x7*_y6   + _x6*_y7   + _x5*_y8  ;
    u_int64_t _z4lo1 = _x4*_y9   + _19x3*_y0 + _19x2*_y1 + _19x1*_y2 + _19x0*_y3;
//    assert(_z4lo0 == _z5432lo0arr[1]);
//    assert(_z4lo1 == _z5432lo1arr[1]);
    // seventh
    u_int64_t _z3lo0 = _x9*_y3 + _2x8*_y4    + _x7*_y5   + _2x6*_y6  + _x5*_y7  ;
    u_int64_t _z3lo1 = _2x4*_y8  + _x3*_y9   + _38x2*_y0 + _19x1*_y1 + _38x0*_y2;
//    assert(_z3lo0 == _z5432lo0arr[2]);
//    assert(_z3lo1 == _z5432lo1arr[2]);
    // eighth
    u_int64_t _z2lo0 = _x9*_y2 + _x8*_y3     + _x7*_y4   + _x6*_y5   + _x5*_y6  ;
    u_int64_t _z2lo1 = _x4*_y7   + _x3*_y8   + _x2*_y9   + _19x1*_y0 + _19x0*_y1;
//    assert(_z2lo0 == _z5432lo0arr[3]);
//    assert(_z2lo1 == _z5432lo1arr[3]);
    // nineth
    u_int64_t _z1lo0 = _x9*_y1 + _2x8*_y2    + _x7*_y3   + _2x6*_y4  + _x5*_y5  ;
    u_int64_t _z1lo1 = _2x4*_y6  + _x3*_y7   + _2x2*_y8  + _x1*_y9   + _38x0*_y0;
//    assert(_z1lo0 == _z10lo0[0]);
//    assert(_z1lo1 == _z10lo1arr[0]);
    // tenth
    u_int64_t _z0lo0 = _x9*_y0 + _x8*_y1     + _x7*_y2   + _x6*_y3   + _x5*_y4  ;
    u_int64_t _z0lo1 = _x4*_y5   + _x3*_y6   + _x2*_y7   + _x1*_y8   + _x0*_y9;
//    assert(_z0lo0 == _z10lo0[1]);
//    assert(_z0lo1 == _z10lo1arr[1]);
  //*/

    // add two halves and store in limb. Saves carryover in next limb
    u_int64_t c_in =0;
    c_in = _addcarry_u64(c_in, _z9876lo0arr[0], _z9876lo1arr[0], &z[9]);
    c_in = _addcarry_u64(c_in, _z9876lo0arr[1], _z9876lo1arr[1], &z[8]);
    c_in = _addcarry_u64(c_in, _z9876lo0arr[2], _z9876lo1arr[2], &z[7]);
    c_in = _addcarry_u64(c_in, _z9876lo0arr[3], _z9876lo1arr[3], &z[6]);
    c_in = _addcarry_u64(c_in, _z5432lo0arr[0], _z5432lo1arr[0], &z[5]);
    c_in = _addcarry_u64(c_in, _z5432lo0arr[1], _z5432lo1arr[1], &z[4]);
    c_in = _addcarry_u64(c_in, _z5432lo0arr[2], _z5432lo1arr[2], &z[3]);
    c_in = _addcarry_u64(c_in, _z5432lo0arr[3], _z5432lo1arr[3], &z[2]);
    c_in = _addcarry_u64(c_in, _z10lo0arr[0], _z10lo1arr[0], &z[1]);
    c_in = _addcarry_u64(c_in, _z10lo0arr[1], _z10lo1arr[1], &z[0]);
    
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
