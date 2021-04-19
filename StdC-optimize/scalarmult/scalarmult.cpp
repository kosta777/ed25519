#include "scalarmult.h"
#include "../../const.h"
#if R25
#include "../radix25/radix25.h"
#else
#include "../radix/radix51.h"
#endif
#include <stdlib.h>
#include <iostream>

typedef radix_t* point_t;

#if R25
short get_bit(const r25_t& e, int ind)
{
  // printf("get_bit::ind: %d\n", ind);
	if(ind < 26)
		return e.limbs[9]>>ind & 1;
	// printf("skip limb9\n");
  ind -= 26;
	if(ind < 25)
		return e.limbs[8]>>ind & 1;
	// printf("skip limb8\n");
	ind -= 25;
	if(ind < 26)
		return e.limbs[7]>>ind & 1;
	// printf("skip limb7\n");
	ind -= 26;
	if(ind < 25)
		return e.limbs[6] >> ind & 1;
	// printf("skip limb6\n");
	ind -= 25;
	if(ind < 26)
		return e.limbs[5] >> ind & 1;
	// printf("skip limb5\n");
	ind -= 26;
	if(ind < 25)
		return e.limbs[4] >> ind & 1;
	// printf("skip limb4\n");
	ind -= 25;
	if(ind < 26)
		return e.limbs[3] >> ind & 1;
	// printf("skip limb3\n");
	ind -= 26;
	if(ind < 25)
		return e.limbs[2] >> ind & 1;
	// printf("skip limb2\n");
	ind -= 25;
	if(ind < 26)
		return e.limbs[1] >> ind & 1;
	// printf("skip limb1\n");
	ind -= 26;
	if(ind <= 26)
  	return e.limbs[0] >> ind & 1;
  printf("ind: %d\n", ind);
  printf("oopsie doopsie in get_bit radix25\nIf you got here that means something is wrong with ind\n");
  exit(1);
}
#else
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
#endif


/*
* Conditional move
* Sets P<-Q if b != 0, otherwise P remains unchanged
*/
void CMOV_point(short b, const point_t P, const point_t Q){

  int64_t m = (int64_t)(0 - b);
  // std::cerr << std::showbase << std::internal << std::setfill('0') << std::hex;
  // std::cerr<< "b: " << b << "  m: " << std::setw(18) << m << std::endl;
  // std::cerr<< "P:\t" << std::setw(18) << P[0].limbs[1] << "\tQ:\t" << std::setw(18) << Q[0].limbs[1] << std::endl;
  // std::cerr<< "\t" << std::setw(18) << P[1].limbs[1] << "\t\t" << std::setw(18) << Q[1].limbs[1] << std::endl;

  radix_t R[2];
  P[0] = (P[0] & ~m)^(Q[0] & m);
  P[1] = (P[1] & ~m)^(Q[1] & m);

  // std::cerr<< "\nP:\t" << std::setw(18) << P[0].limbs[1] << std::endl;
  // std::cerr<< "\t" << std::setw(18) << P[1].limbs[1] << std::endl; 
}

void edwards(const radix_t *P, const radix_t *Q, radix_t *outX, radix_t *outY) {
  radix_t x1 = P[0];
  radix_t y1 = P[1];
  radix_t x2 = Q[0];
  radix_t y2 = Q[1];
  
  radix_t x3 = x1*y2 + x2*y1;
  x3 = x3 / (ONE + (CURVE_D * x1 * x2 * y1 * y2));
  radix_t y3 = y1*y2 + x1*x2;
  y3 = y3 / (ONE - (CURVE_D * x1 * x2 * y1 * y2));
  
  // Modulo taken care of when doing Radix51 ops
  *outX = x3;
  *outY = y3;
}

//Montgomery's Ladder
radix_t* scalar_mult_unopt(const radix_t *P,const radix_t e, radix_t *P0)
{
	//Should initialize to 0s
//  radix_t *P0 = (radix_t *) calloc(10, sizeof(ulli)); // Need 0 init // initialised in front and freed after
#if R25
  P0[1].limbs[9] = 1;
#else
  P0[1].limbs[4] = 1;
#endif
  //Create a copy of the passed point P because P1 will be edited
#if R25
  radix_t *P1 = (radix_t *) malloc(20 * sizeof(ulli));
#else
  radix_t *P1 = (radix_t *) malloc(10 * sizeof(ulli));
#endif
  P1[0].limbs[0] = P[0].limbs[0];
  P1[0].limbs[1] = P[0].limbs[1];
  P1[0].limbs[2] = P[0].limbs[2];
  P1[0].limbs[3] = P[0].limbs[3];
  P1[0].limbs[4] = P[0].limbs[4];
#if R25
  P1[0].limbs[5] = P[0].limbs[5];
  P1[0].limbs[6] = P[0].limbs[6];
  P1[0].limbs[7] = P[0].limbs[7];
  P1[0].limbs[8] = P[0].limbs[8];
  P1[0].limbs[9] = P[0].limbs[9];
#endif
  P1[1].limbs[0] = P[1].limbs[0];
  P1[1].limbs[1] = P[1].limbs[1];
  P1[1].limbs[2] = P[1].limbs[2];
  P1[1].limbs[3] = P[1].limbs[3];
  P1[1].limbs[4] = P[1].limbs[4];
#if R25
  P1[1].limbs[5] = P[1].limbs[5];
  P1[1].limbs[6] = P[1].limbs[6];
  P1[1].limbs[7] = P[1].limbs[7];
  P1[1].limbs[8] = P[1].limbs[8];
  P1[1].limbs[9] = P[1].limbs[9];
#endif


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
	
  free(P1);
	return P0;
}

//Montgomery's Ladder
radix_t* scalar_mult_branchless(const radix_t *P,const radix_t e, radix_t *P0)
{
  // prepare accumulators
  radix_t _P0x_acc(ulli_ZERO);
  radix_t _P0y_acc(ulli_ONE);
  radix_t _P1x;
  radix_t _P1y;
  radix_t _P2x;
  radix_t _P2y;

  //Create a copy of the passed point P because P1 will be edited
  _P1x = P[0];
  _P1y = P[1];

  short b = 0;
  short ki1 = 0;
  short ki, ki_1;
	for(int i=255;i>=1;i-=2)
	{
    ki = get_bit(e,i);
    b = ki^ki1;
    /**********************************/
    // edwards(P0, P1, &P2[0], &P2[1]); // P2 <- P0 + P1
    // grouped by independent operations
    radix_t _t0 = _P0x_acc * _P1x;
    radix_t _t1 = _P0y_acc * _P1y;
    radix_t _t4 = _P0x_acc * _P1y;
    radix_t _t5 = _P1x * _P0y_acc;
    
    radix_t _x3 = _t4 + _t5;
    radix_t _y3 = _t1 + _t0;
    radix_t _t2 = _t0 * _t1;

    radix_t _t3 = (CURVE_D * _t2);

    _P2x = _x3 / (ONE + _t3);
    _P2y = _y3 / (ONE - _t3);
    /*********************************/

    /*********************************/
    // CMOV_point(b, P0, P1); // P0 <- CMOV(b, P0, P1)
    int64_t m = (int64_t)(0-b);
    _P0x_acc = (_P0x_acc & ~m)^(_P1x & m);
    _P0y_acc = (_P0y_acc & ~m)^(_P1y & m);
    /**********************************/

    /**********************************/
    // edwards(P0, P0, &P0[0], &P0[1]); // P0 <- 2*P0
    // grouped by dependencies
    radix_t _t6 = _P0x_acc.sqr();
    radix_t _t7 = _P0y_acc.sqr();
    radix_t _t10 = _P0x_acc*_P0y_acc;

    radix_t _t9 = (CURVE_D * _t6 * _t7);
    radix_t _x4 = _t10 + _t10;
    radix_t _y4 = _t7 + _t6;

    _P0x_acc = _x4 / (ONE + _t9);
    _P0y_acc = _y4 / (ONE - _t9);
    /**********************************/

    /**********************************/
    _P1x = _P2x; // P1 <- P2
    _P1y = _P2y; // P1 <- P2
    /**********************************/
    
    /**********************************************************/
    ki_1 = get_bit(e,i-1);
    b = ki_1^ki;
    /**********************************/
    // edwards(P0, P1, &P2[0], &P2[1]); // P2 <- P0 + P1
    // grouped by independent operations
    _t0 = _P0x_acc * _P1x;
    _t1 = _P0y_acc * _P1y;
    _t4 = _P0x_acc * _P1y;
    _t5 = _P1x * _P0y_acc;
    
    _x3 = _t4 + _t5;
    _y3 = _t1 + _t0;
    _t2 = _t0 * _t1;

    _t3 = (CURVE_D * _t2);

    _P2x = _x3 / (ONE + _t3);
    _P2y = _y3 / (ONE - _t3);
    /*********************************/

    /*********************************/
    // CMOV_point(b, P0, P1); // P0 <- CMOV(b, P0, P1)
    m = (int64_t)(0-b);
    _P0x_acc = (_P0x_acc & ~m)^(_P1x & m);
    _P0y_acc = (_P0y_acc & ~m)^(_P1y & m);
    /**********************************/

    /**********************************/
    // edwards(P0, P0, &P0[0], &P0[1]); // P0 <- 2*P0
    // grouped by dependencies
    _t6 = _P0x_acc.sqr();
    _t7 = _P0y_acc.sqr();
    _t10 = _P0x_acc*_P0y_acc;

    _t9 = (CURVE_D * _t6 * _t7);
    _x4 = _t10 + _t10;
    _y4 = _t7 + _t6;

    _P0x_acc = _x4 / (ONE + _t9);
    _P0y_acc = _y4 / (ONE - _t9);
    /**********************************/

    /**********************************/
    _P1x = _P2x; // P1 <- P2
    _P1y = _P2y; // P1 <- P2
    /**********************************/
    
    ki1 = ki_1; // k_i+1 <- k_i
  }
  /**********************************/
  // CMOV_point(ki, P0, P1); // return Q <- CMOV(k_0, P0, P1)
  int64_t m = (int64_t)(0-ki);
  _P0x_acc = (_P0x_acc & ~m)^(_P1x & m);
  _P0y_acc = (_P0y_acc & ~m)^(_P1y & m);
  /**********************************/
  // store
  P0[0] = _P0x_acc;
  P0[1] = _P0y_acc;
	return P0;
}

//Montgomery's Ladder
radix_t* scalar_mult(const radix_t *P,const radix_t e, radix_t *P0)
{
  //Should initialize to 0s
//  radix_t *P0 = (radix_t *) calloc(10, sizeof(ulli)); // Need 0 init // initialised in front and freed after
#if R25
  P0[1].limbs[9] = 1;
#else
  P0[1].limbs[4] = 1;
#endif  //Create a copy of the passed point P because P1 will be edited
  
  limb_t temp[] = {0,0,0,0,0};
  radix_t P00acc(temp);
  radix_t P01acc(temp);
#if R25
  P01acc.limbs[9] = 1;
#else
  P01acc.limbs[4] = 1;
#endif
  radix_t P10acc(temp);
  radix_t P11acc(temp);

  P10acc.limbs[0] = P[0].limbs[0];
  P10acc.limbs[1] = P[0].limbs[1];
  P10acc.limbs[2] = P[0].limbs[2];
  P10acc.limbs[3] = P[0].limbs[3];
  P10acc.limbs[4] = P[0].limbs[4];
#if R25
  P10acc.limbs[5] = P[0].limbs[5];
  P10acc.limbs[6] = P[0].limbs[6];
  P10acc.limbs[7] = P[0].limbs[7];
  P10acc.limbs[8] = P[0].limbs[8];
  P10acc.limbs[9] = P[0].limbs[9];
#endif
  P11acc.limbs[0] = P[1].limbs[0];
  P11acc.limbs[1] = P[1].limbs[1];
  P11acc.limbs[2] = P[1].limbs[2];
  P11acc.limbs[3] = P[1].limbs[3];
  P11acc.limbs[4] = P[1].limbs[4];
#if R25
  P11acc.limbs[5] = P[1].limbs[5];
  P11acc.limbs[6] = P[1].limbs[6];
  P11acc.limbs[7] = P[1].limbs[7];
  P11acc.limbs[8] = P[1].limbs[8];
  P11acc.limbs[9] = P[1].limbs[9];
#endif


  for(int i=256;i>=2;i-=2)
  {
    //edwards(p0, p1) > x3, y3
    radix_t p00mulp11 = P00acc * P11acc;
    radix_t p10mulp01 = P10acc * P01acc;
    radix_t allmul = p00mulp11*p10mulp01;
    radix_t Dallmul = CURVE_D * allmul;
    radix_t x3 = (p00mulp11 + p10mulp01) / (ONE + Dallmul);
    radix_t y3 = (P01acc*P11acc + P00acc*P10acc) / (ONE - Dallmul);

    if(get_bit(e, i) == 0)
    {
      //edwards(P0, P1, &P1[0], &P1[1]); // P1 <- P0 + P1
      P10acc = x3;
      P11acc = y3;

      //edwards(P0, P0, &P0[0], &P0[1]); // P0 <- 2 * P0
      radix_t p00sqr = P00acc.sqr();
      radix_t p01sqr = P01acc.sqr();
      radix_t p00mulp01 = P00acc * P01acc;
      Dallmul = CURVE_D * p00sqr * p01sqr;
      
      x3 = (p00mulp01 + p00mulp01) / (ONE + Dallmul);
      y3 = (p01sqr + p00sqr) / (ONE - Dallmul);
      P00acc = x3;
      P01acc = y3;
    }
    else
    {
      //edwards(P0, P1, &P0[0], &P0[1]); // P0 <- P0 + P1
      P00acc = x3;
      P01acc = y3;
      //edwards(P1, P1, &P1[0], &P1[1]); // P1 <- 2 * P1
      radix_t p10sqr = P10acc.sqr();
      radix_t p11sqr = P11acc.sqr();
      radix_t p10mulp11 = P10acc * P11acc;
      Dallmul = CURVE_D * p10sqr * p11sqr;
      
      x3 = (p10mulp11 + p10mulp11) / (ONE + Dallmul);
      y3 = (p11sqr + p10sqr) / (ONE - Dallmul);

      P10acc = x3;
      P11acc = y3;
    }
  
    p00mulp11 = P00acc * P11acc;
    p10mulp01 = P10acc * P01acc;
    allmul = p00mulp11*p10mulp01;
    Dallmul = CURVE_D * allmul;
    x3 = (p00mulp11 + p10mulp01) / (ONE + Dallmul);
    y3 = (P01acc*P11acc + P00acc*P10acc) / (ONE - Dallmul);

    if(get_bit(e, i-1) == 0)
    {
      //edwards(P0, P1, &P1[0], &P1[1]); // P1 <- P0 + P1
      P10acc = x3;
      P11acc = y3;

      //edwards(P0, P0, &P0[0], &P0[1]); // P0 <- 2 * P0
      radix_t p00sqr = P00acc.sqr();
      radix_t p01sqr = P01acc.sqr();
      radix_t p00mulp01 = P00acc * P01acc;
      Dallmul = CURVE_D * p00sqr * p01sqr;
      
      x3 = (p00mulp01 + p00mulp01) / (ONE + Dallmul);
      y3 = (p01sqr + p00sqr) / (ONE - Dallmul);
      P00acc = x3;
      P01acc = y3;
    }
    else
    {
      //edwards(P0, P1, &P0[0], &P0[1]); // P0 <- P0 + P1
      P00acc = x3;
      P01acc = y3;
      //edwards(P1, P1, &P1[0], &P1[1]); // P1 <- 2 * P1
      radix_t p10sqr = P10acc.sqr();
      radix_t p11sqr = P11acc.sqr();
      radix_t p10mulp11 = P10acc * P11acc;
      Dallmul = CURVE_D * p10sqr * p11sqr;
      
      x3 = (p10mulp11 + p10mulp11) / (ONE + Dallmul);
      y3 = (p11sqr + p10sqr) / (ONE - Dallmul);

      P10acc = x3;
      P11acc = y3;
    }


  }
  //last peeled iteration
  radix_t p00mulp11 = P00acc * P11acc;
  radix_t p10mulp01 = P10acc * P01acc;
  radix_t allmul = p00mulp11*p10mulp01;
  radix_t Dallmul = CURVE_D * allmul;
  radix_t x3 = (p00mulp11 + p10mulp01) / (ONE + Dallmul);
  radix_t y3 = (P01acc*P11acc + P00acc*P10acc) / (ONE - Dallmul);

  if(get_bit(e, 0) == 0)
  {
    //edwards(P0, P1, &P1[0], &P1[1]); // P1 <- P0 + P1
    P10acc = x3;
    P11acc = y3;

    //edwards(P0, P0, &P0[0], &P0[1]); // P0 <- 2 * P0
    radix_t p00sqr = P00acc.sqr();
    radix_t p01sqr = P01acc.sqr();
    radix_t p00mulp01 = P00acc * P01acc;
    Dallmul = CURVE_D * p00sqr * p01sqr;
    
    x3 = (p00mulp01 + p00mulp01) / (ONE + Dallmul);
    y3 = (p01sqr + p00sqr) / (ONE - Dallmul);
    P00acc = x3;
    P01acc = y3;
  }
  else
  {
    //edwards(P0, P1, &P0[0], &P0[1]); // P0 <- P0 + P1
    P00acc = x3;
    P01acc = y3;
    //edwards(P1, P1, &P1[0], &P1[1]); // P1 <- 2 * P1
    radix_t p10sqr = P10acc.sqr();
    radix_t p11sqr = P11acc.sqr();
    radix_t p10mulp11 = P10acc * P11acc;
    Dallmul = CURVE_D * p10sqr * p11sqr;
    
    x3 = (p10mulp11 + p10mulp11) / (ONE + Dallmul);
    y3 = (p11sqr + p10sqr) / (ONE - Dallmul);

    P10acc = x3;
    P11acc = y3;
  }


  P0[0].limbs[0] = P00acc.limbs[0];
  P0[0].limbs[1] = P00acc.limbs[1];
  P0[0].limbs[2] = P00acc.limbs[2];
  P0[0].limbs[3] = P00acc.limbs[3];
  P0[0].limbs[4] = P00acc.limbs[4];
#if R25
  P0[0].limbs[5] = P00acc.limbs[5];
  P0[0].limbs[6] = P00acc.limbs[6];
  P0[0].limbs[7] = P00acc.limbs[7];
  P0[0].limbs[8] = P00acc.limbs[8];
  P0[0].limbs[9] = P00acc.limbs[9];
#endif
  P0[1].limbs[0] = P01acc.limbs[0];
  P0[1].limbs[1] = P01acc.limbs[1];
  P0[1].limbs[2] = P01acc.limbs[2];
  P0[1].limbs[3] = P01acc.limbs[3];
  P0[1].limbs[4] = P01acc.limbs[4];
#if R25
  P0[1].limbs[5] = P01acc.limbs[5];
  P0[1].limbs[6] = P01acc.limbs[6];
  P0[1].limbs[7] = P01acc.limbs[7];
  P0[1].limbs[8] = P01acc.limbs[8];
  P0[1].limbs[9] = P01acc.limbs[9];
#endif
  return P0;
}

