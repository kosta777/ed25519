#include "radix25_mul.h"
#include <immintrin.h>

void _r25_mul(const limb_t* x, const limb_t* y, limb_t* z){
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

  // first limb
  z[9] = _x9*_y9 + _38x8*_y0   + _19x7*_y1 + _38x6*_y2 + _19x5*_y3 + _38x4*_y4 + _19x3*_y5 + _38x2*_y6 + _19x1*_y7 + _38x0*_y8;
  // second
  z[8] = _x9*_y8 + _x8*_y9     + _19x7*_y0 + _19x6*_y1 + _19x5*_y2 +  _19x4*_y3 + _19x3*_y4 + _19x2*_y5 + _19x1*_y6 + _19x0*_y7;
  // third
  z[7] = _x9*_y7 + _2x8*_y8    + _x7*_y9   + _38x6*_y0 + _19x5*_y1 + _38x4*_y2 + _19x3*_y3 + _38x2*_y4 + _19x1*_y5 + _38x0*_y6;
  // fourth
  z[6] = _x9*_y6 + _x8*_y7     + _x7*_y8   + _x6*_y9   + _19x5*_y0 + _19x4*_y1 + _19x3*_y2 + _19x2*_y3 + _19x1*_y4 + _19x0*_y5;
  // fifth
  z[5] = _x9*_y5 + _2x8*_y6    + _x7*_y7   + _2x6*_y8  + _x5*_y9 + _38x4*_y0 + _19x3*_y1 + _38x2*_y2 + _19x1*_y3 + _38x0*_y4;
  // sixth
  z[4] = _x9*_y4 + _x8*_y5     + _x7*_y6   + _x6*_y7   + _x5*_y8 + _x4*_y9   + _19x3*_y0 + _19x2*_y1 + _19x1*_y2 + _19x0*_y3;
  // seventh
  z[3] = _x9*_y3 + _2x8*_y4    + _x7*_y5   + _2x6*_y6  + _x5*_y7 + _2x4*_y8  + _x3*_y9   + _38x2*_y0 + _19x1*_y1 + _38x0*_y2;
  // eighth
  z[2] = _x9*_y2 + _x8*_y3     + _x7*_y4   + _x6*_y5   + _x5*_y6 + _x4*_y7   + _x3*_y8   + _x2*_y9   + _19x1*_y0 + _19x0*_y1;
  // nineth
  z[1] = _x9*_y1 + _2x8*_y2    + _x7*_y3   + _2x6*_y4  + _x5*_y5 + _2x4*_y6  + _x3*_y7   + _2x2*_y8  + _x1*_y9   + _38x0*_y0;
  // tenth
  z[0] = _x9*_y0 + _x8*_y1     + _x7*_y2   + _x6*_y3   + _x5*_y4 + _x4*_y5   + _x3*_y6   + _x2*_y7   + _x1*_y8   + _x0*_y9;

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