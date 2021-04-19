#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "radix51_addsub.h"

/* z = x + y */
void _r51_add(const limb_t* x, const limb_t* y, limb_t* z){
  // add limbs
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
  z[3] = x[3] + y[3];
  z[4] = x[4] + y[4];

  // carry over and normalize
  // _carry_and_normalize(z);
}

void _r51_sub(const limb_t* x, const limb_t* y, limb_t* z){
  z[0] = (x[0] + 0xFFFFFFFFFFFFE) - y[0];
  z[1] = (x[1] + 0xFFFFFFFFFFFFE) - y[1];
  z[2] = (x[2] + 0xFFFFFFFFFFFFE) - y[2];
  z[3] = (x[3] + 0xFFFFFFFFFFFFE) - y[3];
  z[4] = (x[4] + 0xFFFFFFFFFFFDA) - y[4];

  // _carry_and_normalize(z);
}

void _r51_neg(const limb_t* x, limb_t* z){
  z[0] = 0x0007ffffffffffff - x[0];
  z[1] = 0x0007ffffffffffff - x[1];
  z[2] = 0x0007ffffffffffff - x[2];
  z[3] = 0x0007ffffffffffff - x[3];
  z[4] = 0x0007ffffffffffed - x[4];

  _carry_and_normalize(z);
}