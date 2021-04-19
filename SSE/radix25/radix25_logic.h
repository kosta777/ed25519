#ifndef RADIX25_LOGIC_H
#define RADIX25_LOGIC_H

#include "radix25_common.h"

void _r25_and(const limb_t *x, const limb_t *y, limb_t *z);
void _r25_xor(const limb_t *x, const limb_t *y, limb_t *z);
void _r25_not(const limb_t *x, limb_t *z);
void _r25_mask(const limb_t *x, const int64_t m, limb_t *z);

void _v_r25_and(const v_limb_t *x, const v_limb_t *y, v_limb_t *z);
void _v_r25_xor(const v_limb_t *x, const v_limb_t *y, v_limb_t *z);
void _v_r25_not(const v_limb_t *x, v_limb_t *z);
void _v_r25_mask(const v_limb_t *x, const __m256i m, v_limb_t *z);

#endif
