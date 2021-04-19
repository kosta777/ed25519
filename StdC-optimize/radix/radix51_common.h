#ifndef RADIX51_COMMON_H
#define RADIX51_COMMON_H

#include "radix51_t.h"
// #include "../const.h"
#include <stdint.h>
#include <stdio.h>


#define RADIXMASK51 0x0007FFFFFFFFFFFF // 2^51 - 1; this is the largest value a (normalized) limb can store
#define RADIXMASK52 0x000FFFFFFFFFFFFF // 2^52 - 1;

#define LOW32OF64   0X00000000FFFFFFFF
#define HIGH32OF64  0xFFFFFFFF00000000


u_int64_t _addc64(const u_int64_t x, const u_int64_t y, u_int64_t* z, u_int64_t c);
void _umul64wide (const u_int64_t a,const u_int64_t b, u_int64_t *chi, u_int64_t *clo);
void _ureduce128wide(limb_t* _z0, limb_t* _z1, limb_t* _z2, limb_t* _z3, limb_t* _z4, limb_t*z);
void _carry_and_normalize(limb_t* z);
void _fix_degenerate(limb_t* z);
#endif