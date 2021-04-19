#include <immintrin.h>
#include "radix51_mul.h"

// z = x*y mod (2^255 -19)
inline void _r51_mul(const limb_t* x, const limb_t* y, limb_t* z){
  limb_t _x0 = x[0]; limb_t _19x0 = 19*_x0;
  limb_t _x1 = x[1]; limb_t _19x1 = 19*_x1;
  limb_t _x2 = x[2]; limb_t _19x2 = 19*_x2;
  limb_t _x3 = x[3]; limb_t _19x3 = 19*_x3;
  limb_t _x4 = x[4];
  limb_t _y0 = y[0];
  limb_t _y1 = y[1];
  limb_t _y2 = y[2];
  limb_t _y3 = y[3];
  limb_t _y4 = y[4];
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
  unsigned char c_in = 0;
  _z40[1] = _mulx_u64(_x4, _y4, &_z40[0]);
  // c = _addc64(_z40[1], _z4[1], &_z4[1], 0);
  // trash = _addc64(_z40[0], _z4[0], &_z4[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z40[1], _z4[1], &_z4[1]);
  _addcarry_u64(c_in, _z40[0], _z4[0], &_z4[0]);
  _z41[1] = _mulx_u64(_19x3, _y0, &_z41[0]);
  // c = _addc64(_z41[1], _z4[1], &_z4[1], 0);
  // trash = _addc64(_z41[0], _z4[0], &_z4[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z41[1], _z4[1], &_z4[1]);
  _addcarry_u64(c_in, _z41[0], _z4[0], &_z4[0]);
  // _umul64wide(_19x2, _y1, &_z42[0], &_z42[1]);
  _z42[1] = _mulx_u64(_19x2, _y1, &_z42[0]);
  // c = _addc64(_z42[1], _z4[1], &_z4[1], 0);
  // trash = _addc64(_z42[0], _z4[0], &_z4[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z42[1], _z4[1], &_z4[1]);
  _addcarry_u64(c_in, _z42[0], _z4[0], &_z4[0]);
  // _umul64wide(_19x1, _y2, &_z43[0], &_z43[1]);
  _z43[1] = _mulx_u64(_19x1, _y2, &_z43[0]);
  // c = _addc64(_z43[1], _z4[1], &_z4[1], 0);
  // trash = _addc64(_z43[0], _z4[0], &_z4[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z43[1], _z4[1], &_z4[1]);
  _addcarry_u64(c_in, _z43[0], _z4[0], &_z4[0]);
  // _umul64wide(_19x0, _y3, &_z44[0], &_z44[1]);
  _z44[1] = _mulx_u64(_19x0, _y3, &_z44[0]);
  // c = _addc64(_z44[1], _z4[1], &_z4[1], 0);
  // trash = _addc64(_z44[0], _z4[0], &_z4[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z44[1], _z4[1], &_z4[1]);
  _addcarry_u64(c_in, _z44[0], _z4[0], &_z4[0]);
  // _z4[0] = _z40[0] +_z41[0] +_z42[0] +_z43[0] +_z44[0];
  // _z4[1] = _z40[1] +_z41[1] +_z42[1] +_z43[1] +_z44[1];

  // z1 = x0y1 + x1y0 + 19x2y4 + 19x3y3 + 19x4y2
  u_int64_t _z30[2]; // {high, low}
  u_int64_t _z31[2];
  u_int64_t _z32[2];
  u_int64_t _z33[2];
  u_int64_t _z34[2];
  limb_t _z3[2] = {0,0}; //{high, low}
  // _umul64wide(_x4, _y3, &_z30[0], &_z30[1]);
  _z30[1] = _mulx_u64(_x4, _y3, &_z30[0]);
  // c = _addc64(_z30[1], _z3[1], &_z3[1], 0);
  // trash = _addc64(_z30[0], _z3[0], &_z3[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z30[1], _z3[1], &_z3[1]);
  _addcarry_u64(c_in, _z30[0], _z3[0], &_z3[0]);
  // _umul64wide(_x3, _y4, &_z31[0], &_z31[1]);
  _z31[1] = _mulx_u64(_x3, _y4, &_z31[0]);
  // c = _addc64(_z31[1], _z3[1], &_z3[1], 0);
  // trash = _addc64(_z31[0], _z3[0], &_z3[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z31[1], _z3[1], &_z3[1]);
  _addcarry_u64(c_in, _z31[0], _z3[0], &_z3[0]);
  // _umul64wide(_19x2, _y0, &_z32[0], &_z32[1]);
  _z32[1] = _mulx_u64(_19x2, _y0, &_z32[0]);
  // c = _addc64(_z32[1], _z3[1], &_z3[1], 0);
  // trash = _addc64(_z32[0], _z3[0], &_z3[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z32[1], _z3[1], &_z3[1]);
  _addcarry_u64(c_in, _z32[0], _z3[0], &_z3[0]);
  // _umul64wide(_19x1, _y1, &_z33[0], &_z33[1]);
  _z33[1] = _mulx_u64(_19x1, _y1, &_z33[0]);
  // c = _addc64(_z33[1], _z3[1], &_z3[1], 0);
  // trash = _addc64(_z33[0], _z3[0], &_z3[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z33[1], _z3[1], &_z3[1]);
  _addcarry_u64(c_in, _z33[0], _z3[0], &_z3[0]);
  // _umul64wide(_19x0, _y2, &_z34[0], &_z34[1]);
  _z34[1] = _mulx_u64(_19x0, _y2, &_z34[0]);
  // c = _addc64(_z34[1], _z3[1], &_z3[1], 0);
  // trash = _addc64(_z34[0], _z3[0], &_z3[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z34[1], _z3[1], &_z3[1]);
  _addcarry_u64(c_in, _z34[0], _z3[0], &_z3[0]);
  // _z3[0] = _z30[0] +_z31[0] +_z32[0] +_z33[0] +_z34[0];
  // _z3[1] = _z30[1] +_z31[1] +_z32[1] +_z33[1] +_z34[1];

  // z2 = x0y2 + x1y1 + x2y0 +19x3y4 + 19x4y3
  u_int64_t _z20[2]; // {high, low}
  u_int64_t _z21[2];
  u_int64_t _z22[2];
  u_int64_t _z23[2];
  u_int64_t _z24[2];
  limb_t _z2[2] = {0,0}; //{high, low}
  // _umul64wide(_x4, _y2, &_z20[0], &_z20[1]);
  _z20[1] = _mulx_u64(_x4, _y2, &_z20[0]);
  // c = _addc64(_z20[1], _z2[1], &_z2[1], 0);
  // trash = _addc64(_z20[0], _z2[0], &_z2[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z20[1], _z2[1], &_z2[1]);
  _addcarry_u64(c_in, _z20[0], _z2[0], &_z2[0]);
  // _umul64wide(_x3, _y3, &_z21[0], &_z21[1]);
  _z21[1] = _mulx_u64(_x3, _y3, &_z21[0]);
  // c = _addc64(_z21[1], _z2[1], &_z2[1], 0);
  // trash = _addc64(_z21[0], _z2[0], &_z2[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z21[1], _z2[1], &_z2[1]);
  _addcarry_u64(c_in, _z21[0], _z2[0], &_z2[0]);
  // _umul64wide(_x2, _y4, &_z22[0], &_z22[1]);
  _z22[1] = _mulx_u64(_x2, _y4, &_z22[0]);
  // c = _addc64(_z22[1], _z2[1], &_z2[1], 0);
  // trash = _addc64(_z22[0], _z2[0], &_z2[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z22[1], _z2[1], &_z2[1]);
  _addcarry_u64(c_in, _z22[0], _z2[0], &_z2[0]);
  // _umul64wide(_19x1, _y0, &_z23[0], &_z23[1]);
  _z23[1] = _mulx_u64(_19x1, _y0, &_z23[0]);
  // c = _addc64(_z23[1], _z2[1], &_z2[1], 0);
  // trash = _addc64(_z23[0], _z2[0], &_z2[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z23[1], _z2[1], &_z2[1]);
  _addcarry_u64(c_in, _z23[0], _z2[0], &_z2[0]);
  // _umul64wide(_19x0, _y1, &_z24[0], &_z24[1]);
  _z24[1] = _mulx_u64(_19x0, _y1, &_z24[0]);
  // c = _addc64(_z24[1], _z2[1], &_z2[1], 0);
  // trash = _addc64(_z24[0], _z2[0], &_z2[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z24[1], _z2[1], &_z2[1]);
  _addcarry_u64(c_in, _z24[0], _z2[0], &_z2[0]);
  // _z2[0] = _z20[0] +_z21[0] +_z22[0] +_z23[0] +_z24[0];
  // _z2[1] = _z20[1] +_z21[1] +_z22[1] +_z23[1] +_z24[1];

  // z3 = x0y3 + x1y2 + x2y1 + x3y0 + 19x4y4
  u_int64_t _z10[2]; // {high, low}
  u_int64_t _z11[2];
  u_int64_t _z12[2];
  u_int64_t _z13[2];
  u_int64_t _z14[2];
  limb_t _z1[2] = {0,0}; //{high, low}
  // _umul64wide(_x4, _y1, &_z10[0], &_z10[1]);
  _z10[1] = _mulx_u64(_x4, _y1, &_z10[0]);
  // c = _addc64(_z10[1], _z1[1], &_z1[1], 0);
  // trash = _addc64(_z10[0], _z1[0], &_z1[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z10[1], _z1[1], &_z1[1]);
  _addcarry_u64(c_in, _z10[0], _z1[0], &_z1[0]);
  // _umul64wide(_x3, _y2, &_z11[0], &_z11[1]);
  _z11[1] = _mulx_u64(_x3, _y2, &_z11[0]);
  // c = _addc64(_z11[1], _z1[1], &_z1[1], 0);
  // trash = _addc64(_z11[0], _z1[0], &_z1[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z11[1], _z1[1], &_z1[1]);
  _addcarry_u64(c_in, _z11[0], _z1[0], &_z1[0]);
  // _umul64wide(_x2, _y3, &_z12[0], &_z12[1]);
  _z12[1] = _mulx_u64(_x2, _y3, &_z12[0]);
  // c = _addc64(_z12[1], _z1[1], &_z1[1], 0);
  // trash = _addc64(_z12[0], _z1[0], &_z1[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z12[1], _z1[1], &_z1[1]);
  _addcarry_u64(c_in, _z12[0], _z1[0], &_z1[0]);
  // _umul64wide(_x1, _y4, &_z13[0], &_z13[1]);
  _z13[1] = _mulx_u64(_x1, _y4, &_z13[0]);
  // c = _addc64(_z13[1], _z1[1], &_z1[1], 0);
  // trash = _addc64(_z13[0], _z1[0], &_z1[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z13[1], _z1[1], &_z1[1]);
  _addcarry_u64(c_in, _z13[0], _z1[0], &_z1[0]);
  // _umul64wide(_19x0, _y0, &_z14[0], &_z14[1]);
  _z14[1] = _mulx_u64(_19x0, _y0, &_z14[0]);
  // c = _addc64(_z14[1], _z1[1], &_z1[1], 0);
  // trash = _addc64(_z14[0], _z1[0], &_z1[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z14[1], _z1[1], &_z1[1]);
  _addcarry_u64(c_in, _z14[0], _z1[0], &_z1[0]);
  // _z1[0] = _z10[0] +_z11[0] +_z12[0] +_z13[0] +_z14[0];
  // _z1[1] = _z10[1] +_z11[1] +_z12[1] +_z13[1] +_z14[1];

  // z4 = x0y4 + x1y3 + x2y2 + x3y1 + x4y0
  u_int64_t _z00[2]; // {high, low}
  u_int64_t _z01[2];
  u_int64_t _z02[2];
  u_int64_t _z03[2];
  u_int64_t _z04[2];
  limb_t _z0[2] = {0,0}; //{high, low}
  // _umul64wide(_x4, _y0, &_z00[0], &_z00[1]);
  _z00[1] = _mulx_u64(_x4, _y0, &_z00[0]);
  // c = _addc64(_z00[1], _z0[1], &_z0[1], 0);
  // trash = _addc64(_z00[0], _z0[0], &_z0[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z00[1], _z0[1], &_z0[1]);
  _addcarry_u64(c_in, _z00[0], _z0[0], &_z0[0]);
  // _umul64wide(_x3, _y1, &_z01[0], &_z01[1]);
  _z01[1] = _mulx_u64(_x3, _y1, &_z01[0]);
  // c = _addc64(_z01[1], _z0[1], &_z0[1], 0);
  // trash = _addc64(_z01[0], _z0[0], &_z0[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z01[1], _z0[1], &_z0[1]);
  _addcarry_u64(c_in, _z01[0], _z0[0], &_z0[0]);
  // _umul64wide(_x2, _y2, &_z02[0], &_z02[1]);
  _z02[1] = _mulx_u64(_x2, _y2, &_z02[0]);
  // c = _addc64(_z02[1], _z0[1], &_z0[1], 0);
  // trash = _addc64(_z02[0], _z0[0], &_z0[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z02[1], _z0[1], &_z0[1]);
  _addcarry_u64(c_in, _z02[0], _z0[0], &_z0[0]);
  // _umul64wide(_x1, _y3, &_z03[0], &_z03[1]);
  _z03[1] = _mulx_u64(_x1, _y3, &_z03[0]);
  // c = _addc64(_z03[1], _z0[1], &_z0[1], 0);
  // trash = _addc64(_z03[0], _z0[0], &_z0[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z03[1], _z0[1], &_z0[1]);
  _addcarry_u64(c_in, _z03[0], _z0[0], &_z0[0]);
  // _umul64wide(_x0, _y4, &_z04[0], &_z04[1]);
  _z04[1] = _mulx_u64(_x0, _y4, &_z04[0]);
  // c = _addc64(_z04[1], _z0[1], &_z0[1], 0);
  // trash = _addc64(_z04[0], _z0[0], &_z0[0], c);
  // c=0;
  c_in = _addcarry_u64(0, _z04[1], _z0[1], &_z0[1]);
  _addcarry_u64(c_in, _z04[0], _z0[0], &_z0[0]);
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
  // _fix_degenerate(z);
}