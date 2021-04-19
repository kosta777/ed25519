#include "radix51_logic.h"

void _r51_and(const limb_t *x, const limb_t *y, limb_t *z) {
    z[0] = x[0] & y[0];
    z[1] = x[1] & y[1];
    z[2] = x[2] & y[2];
    z[3] = x[3] & y[3];
    z[4] = x[4] & y[4];
}

void _r51_xor(const limb_t *x, const limb_t *y, limb_t *z) {
    z[0] = x[0] ^ y[0];
    z[1] = x[1] ^ y[1];
    z[2] = x[2] ^ y[2];
    z[3] = x[3] ^ y[3];
    z[4] = x[4] ^ y[4];
}

void _r51_not(const limb_t *x, limb_t *z) {
    z[0] = ~x[0];
    z[1] = ~x[1];
    z[2] = ~x[2];
    z[3] = ~x[3];
    z[4] = ~x[4];
}

void _r51_mask(const limb_t *x, const int64_t m, limb_t *z) {
    z[0] = m & x[0];
    z[1] = m & x[1];
    z[2] = m & x[2];
    z[3] = m & x[3];
    z[4] = m & x[4];
}