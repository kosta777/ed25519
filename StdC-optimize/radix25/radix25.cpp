#include "radix25.h"

#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include "radix25_t.h"
#include "radix25_common.h"

std::string str_rep(const r25_t *r) {
  ulli r64[4];
  radixtoulli(r, &r64[0]);
  std::stringstream ss;
  for (int i = 0 ; i < 4 ; i++) {
    ss << std::hex << r64[i];
  }
  free(r64);
  return ss.str();
}

/* Translates number in radix64 to radix25 
*  a[0] should always be 0 as a consequence
*/
void _init_uint64a(const ulli* a, limb_t* limbs){
  if(a[0] > 0) {
    std::cerr << "Tried to convert to radix25 a number larger than 256 bits!" << std::endl;
    std::cerr << std::showbase << std::internal << std::setfill('0') << std::hex;
    std::cerr << std::setw(18) << a[0] << " "
              << std::setw(18) << a[1] << " "
              << std::setw(18) << a[2] << " "
              << std::setw(18) << a[3] << " "
              << std::setw(18) << a[4] << std::endl;
    exit(1);
  }
  limbs[9] = a[4]         & RADIXMASK26;
  limbs[8] = (a[4] >> 26) & RADIXMASK25;
  limbs[7] = ((a[3] << 13) | (a[4] >> 51)) & RADIXMASK26;
  limbs[6] = (a[3] >> 13) & RADIXMASK25;
  limbs[5] = (a[3] >> 38) & RADIXMASK26;
  limbs[4] = a[2]         & RADIXMASK25;
  limbs[3] = (a[2] >> 25) & RADIXMASK26;
  limbs[2] = ((a[1] << 13) | (a[2] >> 51)) & RADIXMASK25;
  limbs[1] = (a[1] >> 12) & RADIXMASK26;
  limbs[0] = (a[1] >> 38) & RADIXMASK25;
  _fix_degenerate(limbs);
}

/* Translates a radix25 interpretation of a number into
* an array of 4 uint64.
*/
ulli *radixtoulli(const r25_t *r, ulli *dst) {
  limb_t l9 = r->limbs[9];
  limb_t l8 = r->limbs[8];
  limb_t l7 = r->limbs[7];
  limb_t l6 = r->limbs[6];
  limb_t l5 = r->limbs[5];
  limb_t l4 = r->limbs[4];
  limb_t l3 = r->limbs[3];
  limb_t l2 = r->limbs[2];
  limb_t l1 = r->limbs[1];
  limb_t l0 = r->limbs[0];
  dst[3] = (l7 << 51) | (l8 << 26) | l9;
  dst[2] = (l5 << 38) | (l6 << 13) | (l7 >> 13);
  dst[1] = (l2 << 51) | (l3 << 25) | l4;
  dst[0] = (l0 << 38) | (l1 << 12) | (l2 >> 13);
  return dst;
}
