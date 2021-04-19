#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "radix51_common.h"
// #include "../const.h"

/* From the Go language:
*  Add64 returns the sum with carry of x, y and carry: sum = x + y + carry.
*  The carry input must be 0 or 1; otherwise the behavior is undefined.
*  The carryOut output is guaranteed to be 0 or 1.
*/
u_int64_t _addc64(const u_int64_t x, const u_int64_t y, u_int64_t* z, u_int64_t c){

  *z = x + y + c;

  return ((x & y) | ((x | y) & ~(*z))) >> 63;
};

// multiplies two 64wide integers and stores result in two registers chi and clo
void _umul64wide (const u_int64_t a,const u_int64_t b, u_int64_t *chi, u_int64_t *clo)
{
  u_int64_t a_lo = (u_int64_t)(u_int32_t)a;
  // u_int64_t a_lo = a & LOW32OF64; // This is slower
  u_int64_t a_hi = a >> 32;
  u_int64_t b_lo = (u_int64_t)(u_int32_t)b;
  // u_int64_t b_lo = b & LOW32OF64;
  u_int64_t b_hi = b >> 32;

  u_int64_t p0 = a_lo * b_lo;
  u_int64_t p1 = a_lo * b_hi;
  u_int64_t p2 = a_hi * b_lo;
  u_int64_t p3 = a_hi * b_hi;

  u_int32_t cy = (u_int32_t)(((p0 >> 32) + (u_int32_t)p1 + (u_int32_t)p2) >> 32);
  // u_int64_t _t0 = p1 & LOW32OF64;
  // u_int64_t _t1 = p2 & LOW32OF64;
  // u_int64_t _t2 = p0 >> 32;
  // u_int64_t _t3 = (_t2 + _t0 + _t1) >> 32;
  // u_int64_t cy = _t3 & LOW32OF64;

  *clo = p0 + (p1 << 32) + (p2 << 32);
  *chi = p3 + (p1 >> 32) + (p2 >> 32) + cy;
}

/* Reduce an integer 128bit wide, divided in 10 registers
*  z0_hi, z0_lo, z1_hi, z1_lo, ..., z4_hi, z4_lo
*  into a single limb_t[5] z
 */
void _ureduce128wide(limb_t* _z0, limb_t* _z1, limb_t* _z2, limb_t* _z3, limb_t* _z4, limb_t*z){
  // prep
  // _z4[0] = (_z4[0] << 13) | (_z4[1] >> 51);
  // _z4[1] &= RADIXMASK51;
  limb_t _t40 = _z4[0]; limb_t _t41 = _z4[1];
  limb_t _r40 = _t40 << 13; limb_t _r41 = _t41 >> 51;
  limb_t _z40 = _r40 | _r41;
  limb_t _z41 = _t41 & RADIXMASK51;
  // z3
  // _z3[0] = (_z3[0] << 13) | (_z3[1] >>51);
  // _z3[1] &= RADIXMASK51;
  // z[3] = _z3[1] + _z4[0];
  limb_t _t30 = _z3[0]; limb_t _t31 = _z3[1];
  limb_t _r30 = _t30 << 13; limb_t _r31 = _t31 >> 51;
  limb_t _z30 = _r30 | _r31;
  limb_t _z31 = _t31 & RADIXMASK51;
  limb_t z3 = _z31 + _z40;
  // z2
  // _z2[0] = (_z2[0] << 13) | (_z2[1] >>51);
  // _z2[1] &= RADIXMASK51;
  // z[2] = _z2[1] + _z3[0];
  limb_t _t20 = _z2[0]; limb_t _t21 = _z2[1];
  limb_t _r20 = _t20 << 13; limb_t _r21 = _t21 >> 51;
  limb_t _z20 = _r20 | _r21;
  limb_t _z21 = _t21 & RADIXMASK51;
  limb_t z2 = _z21 + _z30;
  // z1
  // _z1[0] = (_z1[0] << 13) | (_z1[1] >>51);
  // _z1[1] &= RADIXMASK51;
  // z[1] = _z1[1] + _z2[0];
  limb_t _t10 = _z1[0]; limb_t _t11 = _z1[1];
  limb_t _r10 = _t10 << 13; limb_t _r11 = _t11 >> 51;
  limb_t _z10 = _r10 | _r11;
  limb_t _z11 = _t11 & RADIXMASK51;
  limb_t z1 = _z11 + _z20;
  // z0
  // _z0[0] = (_z0[0] << 13) | (_z0[1] >>51);
  // _z0[1] &= RADIXMASK51;
  // z[0] = _z0[1] + _z1[0];
  limb_t _t00 = _z0[0]; limb_t _t01 = _z0[1];
  limb_t _r00 = _t00 << 13; limb_t _r01 = _t01 >> 51;
  limb_t _z00 = _r00 | _r01;
  limb_t _z01 = _t01 & RADIXMASK51;
  limb_t z0 = _z01 + _z10;
  // z4
  // z[4] = _z4[1] + 19*_z0[0];
  limb_t z4 = _z41 + 19*_z00;

  // store
  z[0] = z0;
  z[1] = z1;
  z[2] = z2;
  z[3] = z3;
  z[4] = z4;
}

void _carry_and_normalize(limb_t* z){
  // ulli tmp;
  // tmp = z[4];
  // tmp = tmp >>51; // isolate carries
  // z[3]+=tmp; // carry over
  // z[4]&=RADIXMASK51; // zero the carries
  limb_t _t4 = z[4]; limb_t _t3 = z[3];
  limb_t _r40 = _t4 >> 51;
  _t3 = _t3 + _r40; //z3
  _t4 = _t4 & RADIXMASK51; //z4

  // tmp = z[3];
  // tmp = tmp >>51;
  // z[2]+=tmp;
  // z[3]&=RADIXMASK51;
  limb_t _t2 = z[2];
  limb_t _r30 = _t3 >> 51;
  _t2 = _t2 + _r30;
  limb_t _z3 = _t3 & RADIXMASK51;

  // tmp = z[2];
  // tmp = tmp >>51;
  // z[1]+=tmp;
  // z[2]&=RADIXMASK51;
  limb_t _t1 = z[1];
  limb_t _r20 = _t2 >> 51;
  _t1 = _t1 + _r20;
  limb_t _z2 = _t2 & RADIXMASK51;

  // tmp = z[1];
  // tmp = tmp >>51;
  // z[0]+=tmp;
  // z[1]&=RADIXMASK51;
  limb_t _t0 = z[0];
  limb_t _r10 = _t1 >> 51;
  _t0 = _t0 + _r10;
  limb_t _z1 = _t1 & RADIXMASK51;

  // tmp = z[0];
  // tmp = tmp >> 51;
  // z[4]+= tmp * 19;
  // z[0]&=RADIXMASK51;
  limb_t _r00 = _t0 >> 51;
  limb_t _z4 = _t4 + _r00*19;
  limb_t _z0 = _t0 & RADIXMASK51;

  //store
  z[0] = _z0;
  z[1] = _z1;
  z[2] = _z2;
  z[3] = _z3;
  z[4] = _z4;

}

void _fix_degenerate(limb_t* z){
  if (z[0]>= 0x0007ffffffffffff &&
      z[1]>= 0x0007ffffffffffff &&
      z[2]>= 0x0007ffffffffffff &&
      z[3]>= 0x0007ffffffffffff &&
      z[4]>= 0x0007ffffffffffed)
  {
    z[4] = z[4] - 0x0007ffffffffffed;
    z[0] = z[1] = z[2] = z[3] = 0;
    // z[4] += 18;
  }
}
