#ifndef RADIX25_ADDSUB_H
#define RADIX25_ADDSUB_H

// #include "radix25.h"
#include "radix25_common.h"
// #include "../const.h"

void _r25_add(const limb_t* x, const limb_t* y, limb_t* z);
void _r25_sub(const limb_t* x, const limb_t* y, limb_t* z);
void _r25_neg(const limb_t* x, limb_t* z);

#endif