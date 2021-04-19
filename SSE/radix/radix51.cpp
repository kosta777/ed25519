#include "radix51.h"
#include "radix51_common.h"
// #include "../const.h"
#include <stdio.h>
#include <iomanip>
#include <iostream>

#include <string>
#include <sstream>

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
