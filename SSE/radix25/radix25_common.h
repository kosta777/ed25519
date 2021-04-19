#ifndef RADIX25_COMMON_H
#define RADIX25_COMMON_H

#include "radix25_t.h"
#include <stdint.h>
#include <stdio.h>

#define RADIXMASK25 0x1ffffff // 2^25 - 1;
#define RADIXMASK26 0x3ffffff // 2^26 - 1;

#define LOW32OF64   0X00000000FFFFFFFF
#define HIGH32OF64  0xFFFFFFFF00000000


u_int64_t _addc64(const u_int64_t x, const u_int64_t y, u_int64_t* z, u_int64_t c);
// void _umul64wide (const u_int64_t a,const u_int64_t b, u_int64_t *chi, u_int64_t *clo);
// void _ureduce128wide(limb_t* _z0, limb_t* _z1, limb_t* _z2, limb_t* _z3, limb_t* _z4, limb_t*z);
void _carry_and_normalize(limb_t* z);
void _v_carry_and_normalize(v_limb_t* z);
void _fix_degenerate(limb_t* z);
void _v_fix_degenerate(v_limb_t *z);
#endif
