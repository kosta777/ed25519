#ifndef SCALARMUL_H
#define SCALARMUL_H

#if R25
#include "../radix25/radix25.h"
#else
#include "../radix/radix51.h"
#endif
// #include "../ed25519.h"

/* Module for elliptic curve scalar mult */
/*
 int shifts: 1
 int ands: 1
 int subs: 4
 */
short get_bit(const radix_t& e, int ind);

/*
 r51 adds: 3
 r51 subs: 1
 r51 muls: 14
 r51 invs: 2
*/
void edwards(const radix_t *P, const radix_t *Q, radix_t *outX, radix_t *outY);

/*
 257 * cost of get_bit()
 + 514 * cost of curve_add()
 */
radix_t* scalar_mult(const radix_t* P,const radix_t e, radix_t *P0);
radix_t* scalar_mult_branchless(const radix_t *P,const radix_t e, radix_t *P0);

/**
    P = 8 (4*2) radixes - actual points
    e = 4 radixes - exponents
    P00, P01, P02, P03 = outputs
 */
void v_scalar_mult(const radix_t P[4][2], const radix_t e[4], radix_t *P00, radix_t *P01, radix_t *P02, radix_t *P03);


#endif
