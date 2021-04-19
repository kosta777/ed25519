#ifndef RADIX51_LOGIC_H
#define RADIX51_LOGIC_H

#include "radix51_common.h"

void _r51_and(const limb_t *x, const limb_t *y, limb_t *z);
void _r51_xor(const limb_t *x, const limb_t *y, limb_t *z);
void _r51_not(const limb_t *x, limb_t *z);
void _r51_mask(const limb_t *x, const int64_t m, limb_t *z);

#endif