#ifndef RADIX25_ADDSUB_H
#define RADIX25_ADDSUB_H

// #include "radix25.h"
#include "radix25_common.h"
// #include "../const.h"

void _r25_add(const limb_t* x, const limb_t* y, limb_t* z);
void _r25_sub(const limb_t* x, const limb_t* y, limb_t* z);
void _r25_neg(const limb_t* x, limb_t* z);

void _v_r25_add(const v_limb_t *x, const v_limb_t *y, v_limb_t *z);
void _v_r25_sub(const v_limb_t *x, const v_limb_t *y, v_limb_t *z);
void _v_r25_neg(const v_limb_t *x, v_limb_t *z);

#endif
