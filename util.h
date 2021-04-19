#ifndef _UTIL_H_
#define _UTIL_H_

#if STDC
#if R25
  #include "StdC-optimize/radix25/radix25.h"
#else
  #include "StdC-optimize/radix/radix51.h"
#endif
#elif BASE
#include "straight/radix51.h"
#elif SSE
#if R25
  #include "SSE/radix25/radix25.h"
#else
  #include "SSE/radix/radix51.h"
#endif
#endif

int a2v(char c);
char v2a(int c);
char *hexlify(char *bstr, int n);
char *unhexlify(const char *hstr, int offset);
radix_t radixify(ulli *num_in_radix64);

void tmp_print_bytes(char *bstr, int n);

#endif
