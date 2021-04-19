#include "../radix51.h"
#include "../const.h"
#include <iostream>
#include <stdio.h>

int main(int argc, char const *argv[])
{
  ulli tempQ[] = {0x0, 0x7fffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffed};
  ulli tempb[] = {0x0,0x0,0x0,0x0,0x100};
  ulli tempL[] = {0x0, 0x1000000000000000, 0x0000000000000000, 0x14def9dea2f79cd6, 0x5812631a5cf5d3ed};

  ulli tempD[] = {0x0, 0x52036cee2b6ffe73, 0x8cc740797779e898, 0x00700a4d4141d8ab, 0x75eb4dca135978a3};
  ulli tempI[] = {0x0, 0x2b8324804fc1df0b, 0x2b4d00993dfbd7a7, 0x2f431806ad2fe478, 0xc4ee1b274a0ea0b0};
  ulli tempBx[] = {0x0, 0x216936d3cd6e53fe, 0xc0a4e231fdd6dc5c, 0x692cc7609525a7b2, 0xc9562d608f25d51a};
  ulli tempBy[] = {0x0, 0x6666666666666666, 0x6666666666666666, 0x6666666666666666, 0x6666666666666658};

  r51_t Q(tempQ);
  r51_t b(tempb);
  r51_t L(tempL);

  r51_t D(tempD);
  r51_t I(tempI);
  r51_t Bx(tempBx);
  r51_t By(tempBy);

  printf("CURVE Q:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", Q.limbs[0], Q.limbs[1], Q.limbs[2], Q.limbs[3], Q.limbs[4]);
  printf("CURVE b:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", b.limbs[0], b.limbs[1], b.limbs[2], b.limbs[3], b.limbs[4]);
  printf("CURVE L:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", L.limbs[0], L.limbs[1], L.limbs[2], L.limbs[3], L.limbs[4]);
  printf("CURVE D:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", D.limbs[0], D.limbs[1], D.limbs[2], D.limbs[3], D.limbs[4]);
  printf("CURVE I:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", I.limbs[0], I.limbs[1], I.limbs[2], I.limbs[3], I.limbs[4]);
  printf("CURVE Bx:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", Bx.limbs[0], Bx.limbs[1], Bx.limbs[2], Bx.limbs[3], Bx.limbs[4]);
  printf("CURVE By:\t{0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx, 0x%016llx};\n", By.limbs[0], By.limbs[1], By.limbs[2], By.limbs[3], By.limbs[4]);
  
  return 0;
}
