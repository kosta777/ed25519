#ifndef RADIX51_ADDSUB_H
#define RADIX51_ADDSUB_H

// #include "radix51.h"
#include "radix51_common.h"
// #include "../const.h"

void _r51_add(const limb_t* x, const limb_t* y, limb_t* z);
void _r51_sub(const limb_t* x, const limb_t* y, limb_t* z);
void _r51_neg(const limb_t* x, limb_t* z);

#endif