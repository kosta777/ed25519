#ifndef RADIX51_H
#define RADIX51_H

#include <stdint.h>
#include <stdio.h>
#include <iostream>

#include <string>
#include <iomanip>
#include <sstream>

#define RADIXMASK51 0x0007FFFFFFFFFFFF // 2^51 - 1; this is the largest value a (normalized) limb can store
#define RADIXMASK52 0x000FFFFFFFFFFFFF // 2^52 - 1;

typedef u_int64_t ulli;
typedef ulli limb_t;
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
typedef struct Radix51 r51_t;

/* Adds two r51_t mod q and stores result in third argument */
void _r51_add(const limb_t*, const limb_t*, limb_t*);
/* Multiplies two r51_t mod q and stores result in third argument */
void _r51_mul(const limb_t*, const limb_t*, limb_t*);
void _r51_div(const limb_t*, const limb_t*, limb_t*);
/* Substract second argument from first mod q and stores result in third argument */
void _r51_sub(const limb_t*, const limb_t*, limb_t*);
/* negate first argument mod q and stores result in second argument */
void _r51_neg(const limb_t*, limb_t*);
/* Squares first argument modulo q and stores result in second argument */
void _r51_sqr(const limb_t* x, limb_t* z);
/* Stores the multiplicative inverse of x in z*/
void _r51_inv(const limb_t* x, limb_t* z);

void _r51_shl(const limb_t *x, const int _a, limb_t *z);
void _r51_shr(const limb_t *x, const int _a, limb_t *z);

/* Translates a radix51 interpretation of a number into
* an array of 4 uint64.
*/
ulli *radixtoulli(const r51_t *r, ulli *dst);

std::string str_rep(const r51_t *r);

/* expmod(xx,(q+3)//8,q) */
/* used by xrecover */
Radix51 expmod1(const limb_t* x);

void _init_uint64a(const ulli* a, limb_t* limbs);

struct Radix51
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
  bool operator<(const Radix51& _a) const{
    for (int i = 4 ; i >= 0 ; i--) { // Starting with the MSB
      if (limbs[i] > _a.limbs[i]) {
        // No way to get smaller than _a now
        return false;
      }
      if (limbs[i] < _a.limbs[i]) {
        // No way to get bigger than _a now
        return true;
      }
    }
    return false; // Reached only if this == _a
  }
  bool operator>(const Radix51& _a) const{
    for (int i = 4 ; i >= 0 ; i--) { // Starting with the MSB
      if (limbs[i] > _a.limbs[i]) {
        // No way to get smaller than _a now
        return true;
      }
      if (limbs[i] < _a.limbs[i]) {
        // No way to get bigger than _a now
        return false;
      }
    }
    return false; // Reached only if this == _a
  }

};

/* DEFINES THE RADIX TYPE UNIVERSALLY */
typedef r51_t radix_t;

#endif
