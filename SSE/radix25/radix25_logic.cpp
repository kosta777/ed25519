#include "radix25_logic.h"

void _r25_and(const limb_t *x, const limb_t *y, limb_t *z) {
    z[0] = x[0] & y[0];
    z[1] = x[1] & y[1];
    z[2] = x[2] & y[2];
    z[3] = x[3] & y[3];
    z[4] = x[4] & y[4];
    z[5] = x[5] & y[5];
    z[6] = x[6] & y[6];
    z[7] = x[7] & y[7];
    z[8] = x[8] & y[8];
    z[9] = x[9] & y[9];
}

void _r25_xor(const limb_t *x, const limb_t *y, limb_t *z) {
    z[0] = x[0] ^ y[0];
    z[1] = x[1] ^ y[1];
    z[2] = x[2] ^ y[2];
    z[3] = x[3] ^ y[3];
    z[4] = x[4] ^ y[4];
    z[5] = x[5] ^ y[5];
    z[6] = x[6] ^ y[6];
    z[7] = x[7] ^ y[7];
    z[8] = x[8] ^ y[8];
    z[9] = x[9] ^ y[9];
}

void _r25_not(const limb_t *x, limb_t *z) {
    z[0] = ~x[0];
    z[1] = ~x[1];
    z[2] = ~x[2];
    z[3] = ~x[3];
    z[4] = ~x[4];
    z[5] = ~x[5];
    z[6] = ~x[6];
    z[7] = ~x[7];
    z[8] = ~x[8];
    z[9] = ~x[9];
}

void _r25_mask(const limb_t *x, const int64_t m, limb_t *z) {
    z[0] = m & x[0];
    z[1] = m & x[1];
    z[2] = m & x[2];
    z[3] = m & x[3];
    z[4] = m & x[4];
    z[5] = m & x[5];
    z[6] = m & x[6];
    z[7] = m & x[7];
    z[8] = m & x[8];
    z[9] = m & x[9];
}
