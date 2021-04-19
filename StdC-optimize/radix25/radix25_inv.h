#ifndef RADIX25_INV_H
#define RADIX25_INV_H

// #include "radix25.h"
#include "radix25_common.h"
// #include "../const.h"

void _r25_inv(const limb_t* x, limb_t* z);
void _r25_div(const limb_t* x, const limb_t* y, limb_t* z);

#endif