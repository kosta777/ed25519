#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <fstream>
// #include "StdC-optimize/radix/radix51.h"

int a2v(char c) {
  if ((c >= '0') && (c <= '9')) return c - '0';
  if ((c >= 'a') && (c <= 'f')) return c - 'a' + 10;
  if ((c >= 'A') && (c <= 'F')) return c - 'A' + 10;
  else return 0;
}

char v2a(int c) {
  const char hex[] = "0123456789abcdef";
  return hex[c];
}

//num_in_radix64 is ulli array with 4 elements
radix_t radixify(ulli *num_in_radix64)
{
  ulli* prep_for_init = (ulli *)malloc(5 * sizeof(ulli));
  prep_for_init[0] = 0;
  memcpy(&prep_for_init[1], num_in_radix64, 4*sizeof(ulli));
  radix_t retval(prep_for_init);
  free(prep_for_init);
  return retval;
}


char *hexlify(char *bstr, int n)
{
  char *hstr = (char *)malloc((n * 2) + 1);
  char *phstr = hstr;
  for (int i = 0; i < n; i++)
  {
      *phstr++ = v2a((bstr[i] >> 4) & 0x0F);
      *phstr++ = v2a((bstr[i]) & 0x0F);
  }
  *phstr++ = '\0';
  return hstr;
}

//Returns ulli array with 4 elements
char *unhexlify(const char *hstr, int offset)
{
  char *bstr = (char *)malloc((strlen(hstr) / 2) + 1);
  char *pbstr = bstr;
  for (int i = 0; i < strlen(hstr); i += 2)
  {
    char c = (a2v(hstr[i]) << 4) + a2v(hstr[i + 1]);
      *pbstr++ = c;
  }
  *pbstr++ = '\0';
  return bstr;
}

void tmp_print_bytes(char *bstr, int n) {
  for (int i = 0 ; i < n ; i++) {
    fprintf(stderr, "%02x", (int)bstr[i] & 0xff);
  }
  fprintf(stderr, "\n");
}
