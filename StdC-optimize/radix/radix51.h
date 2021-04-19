#ifndef RADIX51_H
#define RADIX51_H

#include <stdint.h>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sstream>

#include "radix51_t.h"
#include "radix51_common.h"
#include "radix51_addsub.h" // _r51_add _r51_sub _r51_neg
#include "radix51_mul.h" // _r51_mul
#include "radix51_sqr.h" // _r51_sqr
#include "radix51_inv.h" // _r51_inv _r51_div
#include "radix51_shift.h" // _r51_shl
#include "radix51_logic.h" // _r51_and _r51_xor _r51_not _r51_mask
// #include "radix51_expmod.h" // _r51_shr _r51_shl

typedef u_int64_t ulli;
typedef ulli limb_t;

void _init_uint64a(const ulli* a, limb_t* limbs);
/** 
 * Implements elements of F_{2^255-19} using
 * data type unsigned int256 with radix-51 representation
 * Implements  x+y: _r51_add
 *             x*y: _r51_mul
 *             x-y: _r51_sub
 *              -x: _r51_neg
 *             x/y: _r51_div
 *           x.sqr: _r51_sqr
 *           x.inv: _r51_inv
 **/

typedef struct Radix51
{
  limb_t limbs[5];

  Radix51() {
    limbs[0] = limbs[1] = limbs[2] = limbs[3] = limbs[4] = 0;
  }
  Radix51(const ulli* _n) {
    _init_uint64a(_n, &limbs[0]);
  }
  Radix51 operator+(const Radix51& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_add(limbs, _a.limbs, z.limbs);
    return z;
  }
  Radix51 operator-(const Radix51& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_sub(limbs, _a.limbs, z.limbs);
    return z;
  }
  Radix51 operator-() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_neg(limbs, z.limbs);
    return z;
  }
  Radix51 operator*(const Radix51& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_mul(limbs, _a.limbs, z.limbs);
    return z;
  }
  Radix51 operator/(const Radix51& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_div(limbs, _a.limbs, z.limbs);
    return z;
  }
  
  Radix51 inv() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_inv(limbs, z.limbs);
    return z;
  }

  Radix51 sqr() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_sqr(limbs, z.limbs);
    return z;
  }

  Radix51 operator<<(const int& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_shl(limbs, _a, z.limbs);
    return z;
  }
  Radix51 operator>>(const int& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_shr(limbs, _a, z.limbs);
    return z;
  }
  Radix51 operator&(const Radix51& _y) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_and(limbs, _y.limbs, z.limbs);
    return z;
  }
  Radix51 operator^(const Radix51& _y) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_xor(limbs, _y.limbs, z.limbs);
    return z;
  }
  Radix51 operator~() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_not(limbs, z.limbs);
    return z;
  }
  Radix51 operator&(const int64_t& _m) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix51 z(temp);
    _r51_mask(limbs, _m, z.limbs);
    return z;
  }

  
  // Comparison operators
  bool operator==(const Radix51& _a) const{
    for (int i = 4 ; i >= 0 ; i--) {
      if (limbs[i] != _a.limbs[i]) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const Radix51& _a) const{
    for (int i = 4 ; i >= 0 ; i--) {
      if (limbs[i] != _a.limbs[i]) {
        return true;
      }
    }
    return false;
  }

} r51_t;

/* Converts a r51_t in an ulli[4] */
ulli *radixtoulli(const r51_t *r, ulli *dst);
std::string str_rep(const r51_t *r);

/* DEFINES THE RADIX TYPE UNIVERSALLY */
typedef r51_t radix_t;


#endif
