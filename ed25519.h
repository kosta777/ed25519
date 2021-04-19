#ifndef ED25519_H
#define ED25519_H
#include <openssl/sha.h>
#include <math.h>
// #include "ed25519.h"
#if STDC
#include "StdC-optimize/scalarmult/scalarmult.h"
  #if R25
  #include "StdC-optimize/radix25/radix25.h"
  #include "StdC-optimize/radix25/radix25_expmod.h"
  #else
  #include "StdC-optimize/radix/radix51.h"
  #include "StdC-optimize/radix/radix51_expmod.h"
  #endif
#elif BASE
#include "straight/scalarmult.h"
#include "straight/radix51.h"
#elif SSE
#include "SSE/scalarmult/scalarmult.h"
  #if R25
  #include "SSE/radix25/radix25.h"
  #include "SSE/radix25/radix25_expmod.h"
  #else
  #include "SSE/radix/radix51.h"
  #include "SSE/radix/radix51_expmod.h"
  #endif
#endif
#include "const.h"

/*
 int divs: 1
 int mods: 1
 int shifts: 1
 int ands: 1
*/
int bit_at(char *m, int i);

/*
 r51 adds: 1
 r51 subs: 3
 r51 muls: 6
 + cost of expmod1()
 int ands: 1
*/
radix_t x_recover(radix_t y);

/*
 int shifts: 512
 int ands: 256
 int adds: 256
*/
void encode_int(radix_t y, char *dst);

/*
 int shifts: 511
 int ands: 256
 int adds: 256
 */
inline void encode_point(radix_t *P, char *dst);


/*
 cost of sha512()
 */
inline void H(const unsigned char *m, int n, char *dst);

/*
 r51 adds: 512
 r51 muls: 512
 r51 shifts: 512
 + 512 * cost of bit_at()
 */
radix_t H_int(const unsigned char *m);


/*
 r51 adds: 252
 r51 subs: 1 (Or as many as the wonky modulo costs)
 r51 muls: 252
 r51 shifts: 251
 + 2 * cost of H_int()
 + cost of scalar_mult()
 + 2 * cost of encode_point()
 + cost of encode_int()
 + cost of unhexlify()
 + cost of radixify()
 */
void signature(char *m, char *sk, char *pk, char *sig);

/*
 r51 adds: 251
 r51 muls: 251
 r51 shifts: 251
 + cost of scalar_mult()
 + cost of encode_point()
 + cost of unhexlify()
 + cost of radixify()
 */
void publickey(char *sk, char *pk);

void publickeys(int n, char sk[][32], char pk[][64]);

#if SSE && R25
void v_publickeys(int n, char sk[][32], char pk[][64]);
#endif














/*
 r51 adds: 1
 r51 subs: 2
 r51_negs: 1
 r51 muls: 6
 */
bool is_on_curve(radix_t *P);

/*
 r51 adds: 256
 r51 muls: 256
 r51 shifts: 256
 */
radix_t decode_int(char *s);

/*
 r51 adds: 255
 r51 subs: 1
 r51 muls: 255
 r51 shifts: 255
 int ands: 1
 + cost of x_recover()
 + cost of bit_at()
 + cost of is_on_curve()
 */
radix_t *decode_point(char *s);


/*
 2 * cost of decode_point()
 + cost of decode_int()
 + cost of H_int()
 + cost of scalar_mult()
 + cost of edwards()
 */
int verify(char *sig, char *m, char *pk);

#endif
