#ifndef RADIX25_H
#define RADIX25_H
  #ifndef R25
  #define R25
  #endif

#include <stdint.h>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sstream>

#include "radix25_t.h"
#include "radix25_common.h"
#include "radix25_addsub.h" // _r25_add _r25_sub _r25_neg
#include "radix25_mul.h" // _r25_mul
#include "radix25_sqr.h" // _r25_sqr
#include "radix25_inv.h" // _r25_inv _r25_div
#include "radix25_shift.h" // _r25_shl
#include "radix25_logic.h" // _r25_and _r25_xor _r25_not _r25_mask

typedef u_int64_t ulli;
typedef ulli limb_t;

void _init_uint64a(const ulli* a, limb_t* limbs);

/** 
 * Implements elements of F_{2^255-19} using
 * data type unsigned int256 with radix-25 representation
 * Implements  x+y: _r25_add
 *             x*y: _r25_mul
 *             x-y: _r25_sub
 *              -x: _r25_neg
 *             x/y: _r25_div
 *           x.sqr: _r25_sqr
 *           x.inv: _r25_inv
 **/
typedef struct Radix25
{
  limb_t limbs[10];

  Radix25(){
    // _init_uint64a(ulli_ZERO, &limbs[0]);
    limbs[0] = limbs[1] = limbs[2] = limbs[3] = limbs[4] =
    limbs[5] = limbs[6] = limbs[7] = limbs[8] = limbs[9] = 0;
  }
  Radix25(const ulli* _n) {
    _init_uint64a(_n, &limbs[0]);
  }
  Radix25 operator+(const Radix25& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_add(limbs, _a.limbs, z.limbs);
    return z;
  }
  Radix25 operator-(const Radix25& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_sub(limbs, _a.limbs, z.limbs);
    return z;
  }
  Radix25 operator-() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_neg(limbs, z.limbs);
    return z;
  }
  Radix25 operator*(const Radix25& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_mul(limbs, _a.limbs, z.limbs);
    return z;
  }
  Radix25 operator/(const Radix25& _a) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_div(limbs, _a.limbs, z.limbs);
    return z;
  }
  
  Radix25 inv() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_inv(limbs, z.limbs);
    return z;
  }

  Radix25 sqr() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_sqr(limbs, z.limbs);
    return z;
  }
  Radix25 operator&(const Radix25& _y) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_and(limbs, _y.limbs, z.limbs);
    return z;
  }
  Radix25 operator^(const Radix25& _y) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_xor(limbs, _y.limbs, z.limbs);
    return z;
  }
  Radix25 operator~() const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_not(limbs, z.limbs);
    return z;
  }
  Radix25 operator&(const int64_t& _m) const{
    limb_t temp[] = {0, 0, 0, 0, 0};
    Radix25 z(temp);
    _r25_mask(limbs, _m, z.limbs);
    return z;
  }

  
  // Comparison operators
  bool operator==(const Radix25& _a) const{
    for (int i = 9 ; i >= 0 ; i--) {
      if (limbs[i] != _a.limbs[i]) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const Radix25& _a) const{
    for (int i = 9 ; i >= 0 ; i--) {
      if (limbs[i] != _a.limbs[i]) {
        return true;
      }
    }
    return false;
  }

} r25_t;

/* Converts a r25_t in an ulli[4] */
ulli *radixtoulli(const r25_t *r, ulli *dst);
std::string str_rep(const r25_t *r);

/* DEFINES THE RADIX TYPE UNIVERSALLY */
typedef r25_t radix_t;

#endif
