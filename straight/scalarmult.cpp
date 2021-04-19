#include "scalarmult.h"
#include "../const.h"
#include "radix51.h"
#include <stdlib.h>
// #include "ed25519.h"

short get_bit(const r51_t& e, int ind)
{
	if(ind < 51)
		return e.limbs[4]>>ind & 1;

	ind -= 51;
	if(ind < 51)
		return e.limbs[3]>>ind & 1;

	ind -= 51;
	if(ind < 51)
		return e.limbs[2]>>ind & 1;

	ind -= 51;
	if(ind < 51)
		return e.limbs[1] >> ind & 1;
	ind -= 51;
	return e.limbs[0]>>ind & 1;
	
}

void edwards(const Radix51 *P, const Radix51 *Q, Radix51 *outX, Radix51 *outY) {
  Radix51 x1 = P[0];
  Radix51 y1 = P[1];
  Radix51 x2 = Q[0];
  Radix51 y2 = Q[1];
  
  Radix51 x3 = x1*y2 + x2*y1;
  x3 = x3 / (ONE + (CURVE_D * x1 * x2 * y1 * y2));
  Radix51 y3 = y1*y2 + x1*x2;
  y3 = y3 / (ONE - (CURVE_D * x1 * x2 * y1 * y2));
  
  // Modulo taken care of when doing Radix51 ops
  *outX = x3;
  *outY = y3;
}

//Montgomery's Ladder
r51_t* scalar_mult(const r51_t *P,const r51_t e, r51_t *P0)
{
	//Should initialize to 0s
  P0[1].limbs[4] = 1;
  //Create a copy of the passed point P because P1 will be edited
//  r51_t *P1 = (r51_t *) malloc(10 * sizeof(ulli));
  radix_t P1[2];
  P1[0].limbs[0] = P[0].limbs[0];
  P1[0].limbs[1] = P[0].limbs[1];
  P1[0].limbs[2] = P[0].limbs[2];
  P1[0].limbs[3] = P[0].limbs[3];
  P1[0].limbs[4] = P[0].limbs[4];

  P1[1].limbs[0] = P[1].limbs[0];
  P1[1].limbs[1] = P[1].limbs[1];
  P1[1].limbs[2] = P[1].limbs[2];
  P1[1].limbs[3] = P[1].limbs[3];
  P1[1].limbs[4] = P[1].limbs[4];


	for(int i=256;i>=0;i--)
	{
		if(get_bit(e, i) == 0)
		{
			edwards(P0, P1, &P1[0], &P1[1]); // P1 <- P0 + P1
			edwards(P0, P0, &P0[0], &P0[1]); // P0 <- 2 * P0
		}
		else
		{
			edwards(P0, P1, &P0[0], &P0[1]); // P0 <- P0 + P1
			edwards(P1, P1, &P1[0], &P1[1]); // P1 <- 2 * P1
		}
	}
	
//  free(P1);
	return P0;
}
