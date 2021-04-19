#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "radix25_addsub.h"

/* z = x + y */
void _r25_add(const limb_t* x, const limb_t* y, limb_t* z){
  // add limbs
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
  z[3] = x[3] + y[3];
  z[4] = x[4] + y[4];
  z[5] = x[5] + y[5];
  z[6] = x[6] + y[6];
  z[7] = x[7] + y[7];
  z[8] = x[8] + y[8];
  z[9] = x[9] + y[9];

  // carry over and normalize
  // _carry_and_normalize(z);
  // _fix_degenerate(z);
}

void _r25_sub(const limb_t* x, const limb_t* y, limb_t* z){
  /* First add 2*q to avoid underflow
  */
  z[0] = (x[0] + 0x3fffffe) - y[0];
  z[1] = (x[1] + 0x7fffffe) - y[1];
  z[2] = (x[2] + 0x3fffffe) - y[2];
  z[3] = (x[3] + 0x7fffffe) - y[3];
  z[4] = (x[4] + 0x3fffffe) - y[4];
  z[5] = (x[5] + 0x7fffffe) - y[5];
  z[6] = (x[6] + 0x3fffffe) - y[6];
  z[7] = (x[7] + 0x7fffffe) - y[7];
  z[8] = (x[8] + 0x3fffffe) - y[8];
  z[9] = (x[9] + 0x7ffffda) - y[9];

  // _carry_and_normalize(z);
  // _fix_degenerate(z);
}

void _r25_neg(const limb_t* x, limb_t* z){
  z[0] = 0x1ffffff - x[0];
  z[1] = 0x3ffffff - x[1];
  z[2] = 0x1ffffff - x[2];
  z[3] = 0x3ffffff - x[3];
  z[4] = 0x1ffffff - x[4];
  z[5] = 0x3ffffff - x[5];
  z[6] = 0x1ffffff - x[6];
  z[7] = 0x3ffffff - x[7];
  z[8] = 0x1ffffff - x[8];
  z[9] = 0x3ffffed - x[9];

  // _carry_and_normalize(z);
  // _fix_degenerate(z);
}
