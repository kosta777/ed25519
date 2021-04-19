#ifndef SCALARMUL_H
#define SCALARMUL_H

#include "radix51.h"

/* Module for elliptic curve scalar mult */

/*
 int shifts: 1
 int ands: 1
 int subs: 4
 */
short get_bit(const r51_t& e, int ind);

/*
 r51 adds: 3
 r51 subs: 1
 r51 muls: 14
 r51 invs: 2
*/
void edwards(const r51_t *P, const r51_t *Q, r51_t *outX, r51_t *outY);

/*
 r51 adds: 3
 r51 subs: 1
 r51 muls: 14
 r51 invs: 2
 */
r51_t* curve_add(const r51_t& x1, const r51_t& x2, const r51_t& y1, const r51_t& y2);

/*
 257 * cost of get_bit()
 + 514 * cost of curve_add()
 */
r51_t* scalar_mult(const r51_t* P,const r51_t e, r51_t *P0);

#endif
