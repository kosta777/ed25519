#include "radix25_shift.h"

// WRONG
void _r25_shl(const limb_t *x, const int _a, limb_t *z) {
  int between = _a / 26;
  int within = _a % 26;
  
  for (int i = 9 - between ; i >= 0 ; i--) { // Starting at the LSB
    int size = 26 - ((i+1)%2);
    between = _a /size;
    within = _a % size;
    z[i] = (x[i + between] << within);
    // msbs = x[i + between] >> (50 - _a);
    long long int mask = size == 26 ? RADIXMASK26 : RADIXMASK25;
    z[i] &= mask;
  }
}

// void _r25_shr(const limb_t *x, const int _a, limb_t *z) {
//   int between = _a / 26;
//   int within = _a % 26;
  
//   int mask = 0;
//   for (int i = 0 ; i < within ; i++) {
//     mask |= (1 << i);
//   }
//   limb_t lsbs = 0;
  
//   for (int i = 0 ; i < between ; i++) {
//     z[i] = 0;
//   }
//   for (int i = between ; i < 5 ; i++) {
//     z[i] = x[i - between] >> within;
//     z[i] |= lsbs << (51 - within);
//     lsbs = x[i - between] & mask;
//   }
// }
