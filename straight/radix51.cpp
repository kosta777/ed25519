#include "radix51.h"
#include "../const.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <string>
#include <sstream>


void _umul64wide (u_int64_t , u_int64_t , u_int64_t*, u_int64_t*);
void _ureduce128wide(limb_t*, limb_t*, limb_t*, limb_t*, limb_t*, limb_t*);
void _carry_and_normalize(limb_t*);
void _init_uint64a(limb_t*, limb_t*);
u_int64_t _addc64(u_int64_t, u_int64_t, u_int64_t*, u_int64_t);


std::string str_rep(const r51_t *r) {
  ulli r64[4];
  radixtoulli(r, &r64[0]);
  std::stringstream ss;
  for (int i = 0 ; i < 4 ; i++) {
    ss << std::hex << r64[i];
  }
  free(r64);
  return ss.str();
}

/* From the Go language:
*  Add64 returns the sum with carry of x, y and carry: sum = x + y + carry.
*  The carry input must be 0 or 1; otherwise the behavior is undefined.
*  The carryOut output is guaranteed to be 0 or 1.
*/
u_int64_t _addc64(const u_int64_t x, const u_int64_t y, u_int64_t* z, u_int64_t c){

  *z = x + y + c;

  return ((x & y) | ((x | y) & ~(*z))) >> 63;
};

void _umul64wide (const u_int64_t a,const u_int64_t b, u_int64_t *chi, u_int64_t *clo)
{
  u_int64_t a_lo = (u_int64_t)(u_int32_t)a;
  u_int64_t a_hi = a >> 32;
  u_int64_t b_lo = (u_int64_t)(u_int32_t)b;
  u_int64_t b_hi = b >> 32;

  u_int64_t p0 = a_lo * b_lo;
  u_int64_t p1 = a_lo * b_hi;
  u_int64_t p2 = a_hi * b_lo;
  u_int64_t p3 = a_hi * b_hi;

  u_int32_t cy = (u_int32_t)(((p0 >> 32) + (u_int32_t)p1 + (u_int32_t)p2) >> 32);

  *clo = p0 + (p1 << 32) + (p2 << 32);
  *chi = p3 + (p1 >> 32) + (p2 >> 32) + cy;
}

void _ureduce128wide(limb_t* _z0, limb_t* _z1, limb_t* _z2, limb_t* _z3, limb_t* _z4, limb_t*z){
  _z4[0] = (_z4[0] << 13) | (_z4[1] >> 51);
  _z4[1] &= RADIXMASK51;
  // z3
  _z3[0] = (_z3[0] << 13) | (_z3[1] >>51);
  _z3[1] &= RADIXMASK51;
  z[3] = _z3[1] + _z4[0];
  // z2
  _z2[0] = (_z2[0] << 13) | (_z2[1] >>51);
  _z2[1] &= RADIXMASK51;
  z[2] = _z2[1] + _z3[0];
  // z1
  _z1[0] = (_z1[0] << 13) | (_z1[1] >>51);
  _z1[1] &= RADIXMASK51;
  z[1] = _z1[1] + _z2[0];
  // z0
  _z0[0] = (_z0[0] << 13) | (_z0[1] >>51);
  _z0[1] &= RADIXMASK51;
  z[0] = _z0[1] + _z1[0];
  // z4
  z[4] = _z4[1] + 19*_z0[0];
}

void _carry_and_normalize(limb_t* z){
  ulli tmp;
  tmp = z[4];
  tmp = tmp >>51; // isolate carries
  z[3]+=tmp; // carry over
  z[4]&=RADIXMASK51; // zero the carries

  tmp = z[3];
  tmp = tmp >>51;
  z[2]+=tmp;
  z[3]&=RADIXMASK51;

  tmp = z[2];
  tmp = tmp >>51;
  z[1]+=tmp;
  z[2]&=RADIXMASK51;

  tmp = z[1];
  tmp = tmp >>51;
  z[0]+=tmp;
  z[1]&=RADIXMASK51;

  tmp = z[0];
  tmp = tmp >> 51;
  z[4]+= tmp * 19;
  z[0]&=RADIXMASK51;
}

/* Translates number in radix64 to radix51 
*  a[0] should always be 0 as a consequence
*/
void _init_uint64a(const ulli* a, limb_t* limbs){
  if(a[0] > 0) {
    std::cerr << "Tried to convert to radix51 a number larger than 256 bits!" << std::endl;
    std::cerr << std::showbase << std::internal << std::setfill('0') << std::hex;
    std::cerr << std::setw(18) << a[0] << " "
              << std::setw(18) << a[1] << " "
              << std::setw(18) << a[2] << " "
              << std::setw(18) << a[3] << " "
              << std::setw(18) << a[4] << std::endl;
    exit(1);
    }

  limbs[4] = (a[4]          + (a[0] >> 13)) & RADIXMASK51; // redundant?
  limbs[3] = ((a[3] << 13)  + (a[4] >> 51)) & RADIXMASK51;
  limbs[2] = ((a[2] << 26)  + (a[3] >> 38)) & RADIXMASK51; // unsure if | or + between the two
  limbs[1] = ((a[1] << 39)  + (a[2] >> 25)) & RADIXMASK51; // same
  limbs[0] = 0x0            + (a[1] >> 12);
}

/* z = x + y */
void _r51_add(const limb_t* x, const limb_t* y, limb_t* z){
  // add limbs
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
  z[3] = x[3] + y[3];
  z[4] = x[4] + y[4];

  // carry over and normalize
  _carry_and_normalize(z);
}

void _fix_degenerate(limb_t* z){
  if (z[0]>= 0x0007ffffffffffff &&
      z[1]>= 0x0007ffffffffffff &&
      z[2]>= 0x0007ffffffffffff &&
      z[3]>= 0x0007ffffffffffff &&
      z[4]>= 0x0007ffffffffffed)
  {
    z[4] = z[4] - CURVE_Q.limbs[4];
    z[0] = z[1] = z[2] = z[3] = 0;
    // z[4] += 18;
  }
}

// z = x*y mod (2^255 -19)
void _r51_mul(const limb_t* x, const limb_t* y, limb_t* z){
  //**********************
  // compute each 128bit limb
  // z0 = x0*y0 + 19x1y1 + 19x2y3 + 19x3y2 + 19x4y1
  u_int64_t _z40[2]; // {high, low}
  u_int64_t _z41[2];
  u_int64_t _z42[2];
  u_int64_t _z43[2];
  u_int64_t _z44[2];
  limb_t _z4[2] = {0,0}; //{high, low}
  u_int64_t c = 0;
  u_int64_t trash;
  _umul64wide(x[4], y[4], &_z40[0], &_z40[1]);
  c = _addc64(_z40[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z40[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(19*x[3], y[0], &_z41[0], &_z41[1]);
  c = _addc64(_z41[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z41[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(19*x[2], y[1], &_z42[0], &_z42[1]);
  c = _addc64(_z42[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z42[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(19*x[1], y[2], &_z43[0], &_z43[1]);
  c = _addc64(_z43[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z43[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(19*x[0], y[3], &_z44[0], &_z44[1]);
  c = _addc64(_z44[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z44[0], _z4[0], &_z4[0], c);
  c=0;
  // _z4[0] = _z40[0] +_z41[0] +_z42[0] +_z43[0] +_z44[0];
  // _z4[1] = _z40[1] +_z41[1] +_z42[1] +_z43[1] +_z44[1];

  // z1 = x0y1 + x1y0 + 19x2y4 + 19x3y3 + 19x4y2
  u_int64_t _z30[2]; // {high, low}
  u_int64_t _z31[2];
  u_int64_t _z32[2];
  u_int64_t _z33[2];
  u_int64_t _z34[2];
  limb_t _z3[2] = {0,0}; //{high, low}
  _umul64wide(x[4], y[3], &_z30[0], &_z30[1]);
  c = _addc64(_z30[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z30[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(x[3], y[4], &_z31[0], &_z31[1]);
  c = _addc64(_z31[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z31[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(19*x[2], y[0], &_z32[0], &_z32[1]);
  c = _addc64(_z32[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z32[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(19*x[1], y[1], &_z33[0], &_z33[1]);
  c = _addc64(_z33[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z33[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(19*x[0], y[2], &_z34[0], &_z34[1]);
  c = _addc64(_z34[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z34[0], _z3[0], &_z3[0], c);
  c=0;
  // _z3[0] = _z30[0] +_z31[0] +_z32[0] +_z33[0] +_z34[0];
  // _z3[1] = _z30[1] +_z31[1] +_z32[1] +_z33[1] +_z34[1];

  // z2 = x0y2 + x1y1 + x2y0 +19x3y4 + 19x4y3
  u_int64_t _z20[2]; // {high, low}
  u_int64_t _z21[2];
  u_int64_t _z22[2];
  u_int64_t _z23[2];
  u_int64_t _z24[2];
  limb_t _z2[2] = {0,0}; //{high, low}
  _umul64wide(x[4], y[2], &_z20[0], &_z20[1]);
  c = _addc64(_z20[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z20[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(x[3], y[3], &_z21[0], &_z21[1]);
  c = _addc64(_z21[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z21[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(x[2], y[4], &_z22[0], &_z22[1]);
  c = _addc64(_z22[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z22[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(19*x[1], y[0], &_z23[0], &_z23[1]);
  c = _addc64(_z23[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z23[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(19*x[0], y[1], &_z24[0], &_z24[1]);
  c = _addc64(_z24[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z24[0], _z2[0], &_z2[0], c);
  c=0;
  // _z2[0] = _z20[0] +_z21[0] +_z22[0] +_z23[0] +_z24[0];
  // _z2[1] = _z20[1] +_z21[1] +_z22[1] +_z23[1] +_z24[1];

  // z3 = x0y3 + x1y2 + x2y1 + x3y0 + 19x4y4
  u_int64_t _z10[2]; // {high, low}
  u_int64_t _z11[2];
  u_int64_t _z12[2];
  u_int64_t _z13[2];
  u_int64_t _z14[2];
  limb_t _z1[2] = {0,0}; //{high, low}
  _umul64wide(x[4], y[1], &_z10[0], &_z10[1]);
  c = _addc64(_z10[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z10[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(x[3], y[2], &_z11[0], &_z11[1]);
  c = _addc64(_z11[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z11[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(x[2], y[3], &_z12[0], &_z12[1]);
  c = _addc64(_z12[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z12[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(x[1], y[4], &_z13[0], &_z13[1]);
  c = _addc64(_z13[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z13[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(19*x[0], y[0], &_z14[0], &_z14[1]);
  c = _addc64(_z14[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z14[0], _z1[0], &_z1[0], c);
  c=0;

  // _z1[0] = _z10[0] +_z11[0] +_z12[0] +_z13[0] +_z14[0];
  // _z1[1] = _z10[1] +_z11[1] +_z12[1] +_z13[1] +_z14[1];

  // z4 = x0y4 + x1y3 + x2y2 + x3y1 + x4y0
  u_int64_t _z00[2]; // {high, low}
  u_int64_t _z01[2];
  u_int64_t _z02[2];
  u_int64_t _z03[2];
  u_int64_t _z04[2];
  limb_t _z0[2] = {0,0}; //{high, low}
  _umul64wide(x[4], y[0], &_z00[0], &_z00[1]);
  c = _addc64(_z00[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z00[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(x[3], y[1], &_z01[0], &_z01[1]);
  c = _addc64(_z01[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z01[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(x[2], y[2], &_z02[0], &_z02[1]);
  c = _addc64(_z02[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z02[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(x[1], y[3], &_z03[0], &_z03[1]);
  c = _addc64(_z03[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z03[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(x[0], y[4], &_z04[0], &_z04[1]);
  c = _addc64(_z04[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z04[0], _z0[0], &_z0[0], c);
  c=0;
  // _z0[0] = _z00[0] +_z01[0] +_z02[0] +_z03[0] +_z04[0];
  // _z0[1] = _z00[1] +_z01[1] +_z02[1] +_z03[1] +_z04[1];

  //**********************
  // reduce limbs to 64bit (limb_t)
  _ureduce128wide(_z0, _z1, _z2, _z3, _z4, z);
  //**********************
  // carry over and normalize
  _carry_and_normalize(z);
  //**********************
  // fix degenerate values
  _fix_degenerate(z);
}

void _r51_sqr(const limb_t* x, limb_t* z){
  // _r51_mul(x,x,z);
  //**********************
  // compute each 128bit limb
  u_int64_t _z40[2]; // {high, low}
  u_int64_t _z41[2];
  u_int64_t _z42[2];
  limb_t _z4[2] = {0,0}; //{high, low}
  u_int64_t c = 0;
  u_int64_t trash;
  _umul64wide(x[4], x[4], &_z40[0], &_z40[1]);
  c = _addc64(_z40[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z40[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(2*19*x[3], x[0], &_z41[0], &_z41[1]);
  c = _addc64(_z41[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z41[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(2*19*x[2], x[1], &_z42[0], &_z42[1]);
  c = _addc64(_z42[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z42[0], _z4[0], &_z4[0], c);
  c=0;
  // _z4[0] = _z40[0] +_z41[0] +_z42[0];
  // _z4[1] = _z40[1] +_z41[1] +_z42[1];

  u_int64_t _z30[2]; // {high, low}
  u_int64_t _z31[2];
  u_int64_t _z32[2];
  limb_t _z3[2] = {0,0}; //{high, low}
  _umul64wide(2*x[4], x[3], &_z30[0], &_z30[1]);
  c = _addc64(_z30[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z30[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(2*19*x[2], x[0], &_z31[0], &_z31[1]);
  c = _addc64(_z31[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z31[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(19*x[1], x[1], &_z32[0], &_z32[1]);
  c = _addc64(_z32[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z32[0], _z3[0], &_z3[0], c);
  c=0;
  // _z3[0] = _z30[0] +_z31[0] +_z32[0];
  // _z3[1] = _z30[1] +_z31[1] +_z32[1];

  u_int64_t _z20[2]; // {high, low}
  u_int64_t _z21[2];
  u_int64_t _z22[2];
  limb_t _z2[2] ={0,0}; //{high, low}
  _umul64wide(2*x[4], x[2], &_z20[0], &_z20[1]);
  c = _addc64(_z20[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z20[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(x[3], x[3], &_z21[0], &_z21[1]);
  c = _addc64(_z21[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z21[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(2*19*x[1], x[0], &_z22[0], &_z22[1]);
  c = _addc64(_z22[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z22[0], _z2[0], &_z2[0], c);
  c=0;
  // _z2[0] = _z20[0] +_z21[0] +_z22[0];
  // _z2[1] = _z20[1] +_z21[1] +_z22[1];

  u_int64_t _z10[2]; // {high, low}
  u_int64_t _z11[2];
  u_int64_t _z12[2];
  limb_t _z1[2] = {0,0}; //{high, low}
  _umul64wide(2*x[4], x[1], &_z10[0], &_z10[1]);
  c = _addc64(_z10[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z10[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(2*x[3], x[2], &_z11[0], &_z11[1]);
  c = _addc64(_z11[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z11[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(19*x[0], x[0], &_z12[0], &_z12[1]);
  c = _addc64(_z12[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z12[0], _z1[0], &_z1[0], c);
  c=0;
  // _z1[0] = _z10[0] +_z11[0] +_z12[0];
  // _z1[1] = _z10[1] +_z11[1] +_z12[1];

  u_int64_t _z00[2]; // {high, low}
  u_int64_t _z01[2];
  u_int64_t _z02[2];
  limb_t _z0[2] = {0,0}; //{high, low}
  _umul64wide(2*x[4], x[0], &_z00[0], &_z00[1]);
  c = _addc64(_z00[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z00[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(2*x[3], x[1], &_z01[0], &_z01[1]);
  c = _addc64(_z01[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z01[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(x[2], x[2], &_z02[0], &_z02[1]);
  c = _addc64(_z02[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z02[0], _z0[0], &_z0[0], c);
  c=0;
  // _z0[0] = _z00[0] +_z01[0] +_z02[0];
  // _z0[1] = _z00[1] +_z01[1] +_z02[1];

  //**********************
  // reduce limbs to 64bit (limb_t)
  _ureduce128wide(_z0, _z1, _z2, _z3, _z4, z);
  //**********************
  // carry over and normalize
  _carry_and_normalize(z);
}

void _r51_sub(const limb_t* x, const limb_t* y, limb_t* z){
  z[0] = (x[0] + 0xFFFFFFFFFFFFE) - y[0];
  z[1] = (x[1] + 0xFFFFFFFFFFFFE) - y[1];
  z[2] = (x[2] + 0xFFFFFFFFFFFFE) - y[2];
  z[3] = (x[3] + 0xFFFFFFFFFFFFE) - y[3];
  z[4] = (x[4] + 0xFFFFFFFFFFFDA) - y[4];

  _carry_and_normalize(z);
}

void _r51_neg(const limb_t* x, limb_t* z){
  z[0] = CURVE_Q.limbs[0] - x[0];
  z[1] = CURVE_Q.limbs[1] - x[1];
  z[2] = CURVE_Q.limbs[2] - x[2];
  z[3] = CURVE_Q.limbs[3] - x[3];
  z[4] = CURVE_Q.limbs[4] - x[4];

  _carry_and_normalize(z);
}

/* this is a port of the original code from Bernstein*/
void _r51_inv(const limb_t* x, limb_t* z){
  limb_t x2[] = {0,0,0,0,0};
  limb_t x9[] = {0,0,0,0,0};
  limb_t x11[] = {0,0,0,0,0};
  limb_t x2_5_0[] = {0,0,0,0,0};
  limb_t x2_10_0[] = {0,0,0,0,0};
  limb_t x2_20_0[] = {0,0,0,0,0};
  limb_t x2_50_0[] = {0,0,0,0,0};
  limb_t x2_100_0[] = {0,0,0,0,0};
  limb_t t0[] = {0,0,0,0,0};
  limb_t t1[] = {0,0,0,0,0};
  int i;

  /* 2 */ _r51_sqr(x, x2);
  /* 4 */ _r51_sqr(x2, t0);
  /* 8 */ _r51_sqr(t0, t1);
  /* 9 */ _r51_mul(t1, x, x9);
  /* 11 */ _r51_mul(x9, x2, x11);
  /* 22 */ _r51_sqr(x11, t0);
  /* 2^5 - 2^0 = 31 */ _r51_mul(t0, x9, x2_5_0);

  /* 2^6 - 2^1 */ _r51_sqr(x2_5_0, t0);
  /* 2^7 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^8 - 2^3 */ _r51_sqr(t1, t0);
  /* 2^9 - 2^4 */ _r51_sqr(t0, t1);
  /* 2^10 - 2^5 */ _r51_sqr(t1, t0);
  /* 2^10 - 2^0 */ _r51_mul(t0, x2_5_0, x2_10_0);

  /* 2^11 - 2^1 */ _r51_sqr(x2_10_0, t0);
  /* 2^12 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^20 - 2^10 */ for (i = 2;i < 10;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^20 - 2^0 */ _r51_mul(t1, x2_10_0, x2_20_0);

  /* 2^21 - 2^1 */ _r51_sqr(x2_20_0, t0);
  /* 2^22 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^40 - 2^20 */ for (i = 2;i < 20;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^40 - 2^0 */ _r51_mul(t1, x2_20_0, t0);

  /* 2^41 - 2^1 */ _r51_sqr(t0, t1);
  /* 2^42 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^50 - 2^10 */ for (i = 2;i < 10;i += 2) { _r51_sqr(t0, t1); _r51_sqr(t1, t0); }
  /* 2^50 - 2^0 */ _r51_mul(t0, x2_10_0, x2_50_0);

  /* 2^51 - 2^1 */ _r51_sqr(x2_50_0, t0);
  /* 2^52 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^100 - 2^50 */ for (i = 2;i < 50;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^100 - 2^0 */ _r51_mul(t1, x2_50_0, x2_100_0);

  /* 2^101 - 2^1 */ _r51_sqr(x2_100_0, t1);
  /* 2^102 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^200 - 2^100 */ for (i = 2;i < 100;i += 2) { _r51_sqr(t0, t1); _r51_sqr(t1, t0); }
  /* 2^200 - 2^0 */ _r51_mul(t0, x2_100_0, t1);

  /* 2^201 - 2^1 */ _r51_sqr(t1, t0);
  /* 2^202 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^250 - 2^50 */ for (i = 2;i < 50;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^250 - 2^0 */ _r51_mul(t1, x2_50_0, t0);

  /* 2^251 - 2^1 */ _r51_sqr(t0, t1);
  /* 2^252 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^253 - 2^3 */ _r51_sqr(t0, t1);
  /* 2^254 - 2^4 */ _r51_sqr(t1, t0);
  /* 2^255 - 2^5 */ _r51_sqr(t0, t1);
  /* 2^255 - 21 */ _r51_mul(t1, x11, z);
}

void _r51_div(const limb_t* x, const limb_t* y, limb_t* z){
  limb_t inv_y[] = {0, 0, 0, 0, 0};
  _r51_inv(y, inv_y);
  _r51_mul(x, inv_y, z);
}

/* expmod(xx,(q+3)//8,q) */
Radix51 expmod1(const limb_t* x){
  limb_t tmp[] = {0,0,0,0,0};
  Radix51 z(tmp);
  limb_t x2[] = {0,0,0,0,0};
  limb_t x9[] = {0,0,0,0,0};
  limb_t x11[] = {0,0,0,0,0};
  limb_t x2_5_0[] = {0,0,0,0,0};
  limb_t x2_10_0[] = {0,0,0,0,0};
  limb_t x2_20_0[] = {0,0,0,0,0};
  limb_t x2_50_0[] = {0,0,0,0,0};
  limb_t x2_100_0[] = {0,0,0,0,0};
  limb_t t0[] = {0,0,0,0,0};
  limb_t t1[] = {0,0,0,0,0};
  int i;

  /* 2 */ _r51_sqr(x, x2);
  /* 4 */ _r51_sqr(x2, t0);
  /* 8 */ _r51_sqr(t0, t1);
  /* 9 */ _r51_mul(t1, x, x9);
  /* 11 */ _r51_mul(x9, x2, x11);
  /* 22 */ _r51_sqr(x11, t0);
  /* 2^5 - 2^0 = 31 */ _r51_mul(t0, x9, x2_5_0);

  /* 2^6 - 2^1 */ _r51_sqr(x2_5_0, t0);
  /* 2^7 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^8 - 2^3 */ _r51_sqr(t1, t0);
  /* 2^9 - 2^4 */ _r51_sqr(t0, t1);
  /* 2^10 - 2^5 */ _r51_sqr(t1, t0);
  /* 2^10 - 2^0 */ _r51_mul(t0, x2_5_0, x2_10_0);

  /* 2^11 - 2^1 */ _r51_sqr(x2_10_0, t0);
  /* 2^12 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^20 - 2^10 */ for (i = 2;i < 10;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^20 - 2^0 */ _r51_mul(t1, x2_10_0, x2_20_0);

  /* 2^21 - 2^1 */ _r51_sqr(x2_20_0, t0);
  /* 2^22 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^40 - 2^20 */ for (i = 2;i < 20;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^40 - 2^0 */ _r51_mul(t1, x2_20_0, t0);

  /* 2^41 - 2^1 */ _r51_sqr(t0, t1);
  /* 2^42 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^50 - 2^10 */ for (i = 2;i < 10;i += 2) { _r51_sqr(t0, t1); _r51_sqr(t1, t0); }
  /* 2^50 - 2^0 */ _r51_mul(t0, x2_10_0, x2_50_0);

  /* 2^51 - 2^1 */ _r51_sqr(x2_50_0, t0);
  /* 2^52 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^100 - 2^50 */ for (i = 2;i < 50;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^100 - 2^0 */ _r51_mul(t1, x2_50_0, x2_100_0);

  /* 2^101 - 2^1 */ _r51_sqr(x2_100_0, t1);
  /* 2^102 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^200 - 2^100 */ for (i = 2;i < 100;i += 2) { _r51_sqr(t0, t1); _r51_sqr(t1, t0); }
  /* 2^200 - 2^0 */ _r51_mul(t0, x2_100_0, t1);

  /* 2^201 - 2^1 */ _r51_sqr(t1, t0);
  /* 2^202 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^250 - 2^50 */ for (i = 2;i < 50;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^250 - 2^0 */ _r51_mul(t1, x2_50_0, t0);

  /* 2^251 - 2^1 */ _r51_sqr(t0, t1);
  /* 2^252 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^252 - 2^1 */ _r51_mul(t0, x2, z.limbs);
  
  return z;
}


void _r51_shl(const limb_t *x, const int _a, limb_t *z) {
  int between = _a / 51;
  int within = _a % 51;
  limb_t msbs = 0;
  
  for (int i = 4 - between ; i >= 0 ; i--) { // Starting at the LSB
    z[i] = (x[i + between] << within);
    msbs = x[i + between] >> (50 - _a);
    z[i] &= RADIXMASK51;
  }
}

void _r51_shr(const limb_t *x, const int _a, limb_t *z) {
  int between = _a / 51;
  int within = _a % 51;
  
  int mask = 0;
  for (int i = 0 ; i < within ; i++) {
    mask |= (1 << i);
  }
  limb_t lsbs = 0;
  
  for (int i = 0 ; i < between ; i++) {
    z[i] = 0;
  }
  for (int i = between ; i < 5 ; i++) {
    z[i] = x[i - between] >> within;
    z[i] |= lsbs << (51 - within);
    lsbs = x[i - between] & mask;
  }
}

/* Translates a radix51 interpretation of a number into
* an array of 4 uint64.
*/
ulli *radixtoulli(const r51_t *r, ulli *dst) {
  dst[3] = (r->limbs[4])      | (r->limbs[3]<<51);
  dst[2] = (r->limbs[3]>>13)  | (r->limbs[2]<<38);
  dst[1] = (r->limbs[2]>>26)  | (r->limbs[1]<<25);
  dst[0] = (r->limbs[1]>>39)  | (r->limbs[0]<<12);
  return dst;
}
