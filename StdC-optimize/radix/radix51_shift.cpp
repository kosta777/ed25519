#include "radix51_shift.h"

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
