#include "radix51_inv.h"
#include "radix51_mul.h"
#include "radix51_sqr.h"

/* this is a port of the original code from Bernstein*/
/* Long chain of operations. All operations depend on the one before for arguments */
void _r51_inv(const limb_t* x, limb_t* z){
  limb_t x2[] = {0,0,0,0,0};
  limb_t x9[] = {0,0,0,0,0};
  limb_t x11[] = {0,0,0,0,0};
  limb_t x2_5_0[] = {0,0,0,0,0};
  limb_t x2_10_0[] = {0,0,0,0,0};
  limb_t x2_20_0[] = {0,0,0,0,0};
  limb_t x2_50_0[] = {0,0,0,0,0};
  limb_t x2_100_0[] = {0,0,0,0,0};
  limb_t t0[] = {0,0,0,0,0};
  limb_t t1[] = {0,0,0,0,0};
  int i;

  /* 2 */ _r51_sqr(x, x2);
  /* 4 */ _r51_sqr(x2, t0);
  /* 8 */ _r51_sqr(t0, t1);
  /* 9 */ _r51_mul(t1, x, x9);
  /* 11 */ _r51_mul(x9, x2, x11);
  /* 22 */ _r51_sqr(x11, t0);
  /* 2^5 - 2^0 = 31 */ _r51_mul(t0, x9, x2_5_0);

  /* 2^6 - 2^1 */ _r51_sqr(x2_5_0, t0);
  /* 2^7 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^8 - 2^3 */ _r51_sqr(t1, t0);
  /* 2^9 - 2^4 */ _r51_sqr(t0, t1);
  /* 2^10 - 2^5 */ _r51_sqr(t1, t0);
  /* 2^10 - 2^0 */ _r51_mul(t0, x2_5_0, x2_10_0);

  /* 2^11 - 2^1 */ _r51_sqr(x2_10_0, t0);
  /* 2^12 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^20 - 2^10 */ for (i = 2;i < 10;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^20 - 2^0 */ _r51_mul(t1, x2_10_0, x2_20_0);

  /* 2^21 - 2^1 */ _r51_sqr(x2_20_0, t0);
  /* 2^22 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^40 - 2^20 */ for (i = 2;i < 20;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^40 - 2^0 */ _r51_mul(t1, x2_20_0, t0);

  /* 2^41 - 2^1 */ _r51_sqr(t0, t1);
  /* 2^42 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^50 - 2^10 */ for (i = 2;i < 10;i += 2) { _r51_sqr(t0, t1); _r51_sqr(t1, t0); }
  /* 2^50 - 2^0 */ _r51_mul(t0, x2_10_0, x2_50_0);

  /* 2^51 - 2^1 */ _r51_sqr(x2_50_0, t0);
  /* 2^52 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^100 - 2^50 */ for (i = 2;i < 50;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^100 - 2^0 */ _r51_mul(t1, x2_50_0, x2_100_0);

  /* 2^101 - 2^1 */ _r51_sqr(x2_100_0, t1);
  /* 2^102 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^200 - 2^100 */ for (i = 2;i < 100;i += 2) { _r51_sqr(t0, t1); _r51_sqr(t1, t0); }
  /* 2^200 - 2^0 */ _r51_mul(t0, x2_100_0, t1);

  /* 2^201 - 2^1 */ _r51_sqr(t1, t0);
  /* 2^202 - 2^2 */ _r51_sqr(t0, t1);
  /* 2^250 - 2^50 */ for (i = 2;i < 50;i += 2) { _r51_sqr(t1, t0); _r51_sqr(t0, t1); }
  /* 2^250 - 2^0 */ _r51_mul(t1, x2_50_0, t0);

  /* 2^251 - 2^1 */ _r51_sqr(t0, t1);
  /* 2^252 - 2^2 */ _r51_sqr(t1, t0);
  /* 2^253 - 2^3 */ _r51_sqr(t0, t1);
  /* 2^254 - 2^4 */ _r51_sqr(t1, t0);
  /* 2^255 - 2^5 */ _r51_sqr(t0, t1);
  /* 2^255 - 21 */ _r51_mul(t1, x11, z);
}

void _r51_div(const limb_t* x, const limb_t* y, limb_t* z){
  limb_t inv_y[] = {0, 0, 0, 0, 0};
  _r51_inv(y, inv_y);
  _r51_mul(x, inv_y, z);
}