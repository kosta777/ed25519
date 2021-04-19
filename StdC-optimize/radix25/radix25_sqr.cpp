#include "radix25_sqr.h"

void _r25_sqr(const limb_t* x, limb_t* z){
  //                uses add to balance ports
  limb_t _x0 = x[0]; limb_t _2x0 = _x0+_x0; limb_t _38x0 = 38*_x0;
  limb_t _x1 = x[1]; limb_t _2x1 = _x1+_x1; limb_t _19x1 = 19*_x1;
  limb_t _x2 = x[2]; limb_t _2x2 = _x2+_x2; limb_t _19x2 = 19*_x2; limb_t _38x2 = 38*_x2;
  limb_t _x3 = x[3]; limb_t _2x3 = _x3+_x3; limb_t _19x3 = 19*_x3;
  limb_t _x4 = x[4]; limb_t _2x4 = _x4+_x4; limb_t _19x4 = 19*_x4; limb_t _38x4 = 38*_x4;
  limb_t _x5 = x[5]; limb_t _2x5 = _x5+_x5; limb_t _19x5 = 19*_x5;
  limb_t _x6 = x[6]; limb_t _2x6 = _x6+_x6; limb_t _19x6 = 19*_x6; limb_t _38x6 = 38*_x6;
  limb_t _x7 = x[7]; limb_t _2x7 = _x7+_x7; limb_t _19x7 = 19*_x7;
  limb_t _x8 = x[8]; limb_t _2x8 = _x8+_x8; limb_t _38x8 = 38*_x8;
  limb_t _x9 = x[9]; limb_t _2x9 = _x9+_x9;

/* Each multiplication produces 58bits numbers.
*  Each limb is a summation of 5-6 products.
*  No nee to split the sum this time as it should always fit in 63 bits
*/
  // first limb
  z[9] = _x9*_x9  + _38x8*_2x0  + _19x7*_2x1 + _38x6*_2x2 + _19x5*_2x3 + _38x4*_x4;
  // second
  z[8] = _2x9*_x8 + _19x7*_2x0  + _19x6*_2x1 + _19x5*_2x2 + _19x4*_2x3;
  // third
  z[7] = _2x9*_x7 + _2x8*_x8    + _38x6*_2x0 + _19x5*_2x1 + _38x4*_2x2 + _19x3*_x3;
  // fourth
  z[6] = _2x9*_x6 + _2x8*_x7    + _19x5*_2x0 + _19x4*_2x1 + _19x3*_2x2;
  // fifth
  z[5] = _2x9*_x5 + _2x8*_2x6   + _x7*_x7    + _38x4*_2x0 + _19x3*_2x1 + _38x2*_x2;
  // sixth
  z[4] = _2x9*_x4 + _2x8*_x5    + _2x7*_x6   + _19x3*_2x0 + _19x2*_2x1;
  // seventh
  z[3] = _2x9*_x3 + _2x8*_2x4   + _2x7*_x5   + _2x6*_x6   + _38x2*_2x0 + _19x1*_x1;
  // eighth
  z[2] = _2x9*_x2 + _2x8*_x3    + _2x7*_x4   + _2x6*_x5   + _19x1*_2x0;
  // nineth
  z[1] = _2x9*_x1 + _2x8*_2x2   + _2x7*_x3   + _2x6*_2x4  + _x5*_x5    + _38x0*_x0;
  // tenth
  z[0] = _2x9*_x0 + _2x8*_x1    + _2x7*_x2   + _2x6*_x3   + _2x5*_x4;
  
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