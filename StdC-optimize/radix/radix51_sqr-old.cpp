#include "radix51_sqr.h"

inline void _r51_sqr(const limb_t* x, limb_t* z){
  limb_t _x0 = x[0]; limb_t _19x0 = 19*_x0;
  limb_t _x1 = x[1]; limb_t _19x1 = 19*_x1;
  limb_t _x2 = x[2];
  limb_t _x3 = x[3];
  limb_t _x4 = x[4];
                     limb_t _38x1 = 38*_x1;
                     limb_t _38x2 = 38*_x2;
  limb_t _2x3 = 2*_x3; limb_t _38x3 = 38*_x3;
  limb_t _2x4 = 2*_x4;
  // _r51_mul(x,x,z);
  //**********************
  // compute each 128bit limb
  u_int64_t _z40[2]; // {high, low}
  u_int64_t _z41[2];
  u_int64_t _z42[2];
  limb_t _z4[2] = {0,0}; //{high, low}
  u_int64_t c = 0;
  u_int64_t trash;
  _umul64wide(_x4, _x4, &_z40[0], &_z40[1]);
  c = _addc64(_z40[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z40[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(_38x3, _x0, &_z41[0], &_z41[1]);
  c = _addc64(_z41[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z41[0], _z4[0], &_z4[0], c);
  c=0;
  _umul64wide(_38x2, _x1, &_z42[0], &_z42[1]);
  c = _addc64(_z42[1], _z4[1], &_z4[1], 0);
  trash = _addc64(_z42[0], _z4[0], &_z4[0], c);
  c=0;
  // _z4[0] = _z40[0] +_z41[0] +_z42[0];
  // _z4[1] = _z40[1] +_z41[1] +_z42[1];

  u_int64_t _z30[2]; // {high, low}
  u_int64_t _z31[2];
  u_int64_t _z32[2];
  limb_t _z3[2] = {0,0}; //{high, low}
  _umul64wide(_2x4, _x3, &_z30[0], &_z30[1]);
  c = _addc64(_z30[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z30[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(_38x2, _x0, &_z31[0], &_z31[1]);
  c = _addc64(_z31[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z31[0], _z3[0], &_z3[0], c);
  c=0;
  _umul64wide(_19x1, _x1, &_z32[0], &_z32[1]);
  c = _addc64(_z32[1], _z3[1], &_z3[1], 0);
  trash = _addc64(_z32[0], _z3[0], &_z3[0], c);
  c=0;
  // _z3[0] = _z30[0] +_z31[0] +_z32[0];
  // _z3[1] = _z30[1] +_z31[1] +_z32[1];

  u_int64_t _z20[2]; // {high, low}
  u_int64_t _z21[2];
  u_int64_t _z22[2];
  limb_t _z2[2] ={0,0}; //{high, low}
  _umul64wide(_2x4, _x2, &_z20[0], &_z20[1]);
  c = _addc64(_z20[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z20[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(_x3, _x3, &_z21[0], &_z21[1]);
  c = _addc64(_z21[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z21[0], _z2[0], &_z2[0], c);
  c=0;
  _umul64wide(_38x1, _x0, &_z22[0], &_z22[1]);
  c = _addc64(_z22[1], _z2[1], &_z2[1], 0);
  trash = _addc64(_z22[0], _z2[0], &_z2[0], c);
  c=0;
  // _z2[0] = _z20[0] +_z21[0] +_z22[0];
  // _z2[1] = _z20[1] +_z21[1] +_z22[1];

  u_int64_t _z10[2]; // {high, low}
  u_int64_t _z11[2];
  u_int64_t _z12[2];
  limb_t _z1[2] = {0,0}; //{high, low}
  _umul64wide(_2x4, _x1, &_z10[0], &_z10[1]);
  c = _addc64(_z10[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z10[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(_2x3, _x2, &_z11[0], &_z11[1]);
  c = _addc64(_z11[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z11[0], _z1[0], &_z1[0], c);
  c=0;
  _umul64wide(_19x0, _x0, &_z12[0], &_z12[1]);
  c = _addc64(_z12[1], _z1[1], &_z1[1], 0);
  trash = _addc64(_z12[0], _z1[0], &_z1[0], c);
  c=0;
  // _z1[0] = _z10[0] +_z11[0] +_z12[0];
  // _z1[1] = _z10[1] +_z11[1] +_z12[1];

  u_int64_t _z00[2]; // {high, low}
  u_int64_t _z01[2];
  u_int64_t _z02[2];
  limb_t _z0[2] = {0,0}; //{high, low}
  _umul64wide(_2x4, _x0, &_z00[0], &_z00[1]);
  c = _addc64(_z00[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z00[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(_2x3, _x1, &_z01[0], &_z01[1]);
  c = _addc64(_z01[1], _z0[1], &_z0[1], 0);
  trash = _addc64(_z01[0], _z0[0], &_z0[0], c);
  c=0;
  _umul64wide(_x2, _x2, &_z02[0], &_z02[1]);
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