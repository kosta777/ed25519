#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "../../catch2/catch.hpp"
#include "../radix25/radix25.h"
#include "../const.h"
// #include "../../util.h"
#include <iostream>
#include <stdio.h>
#include "../scalarmult/scalarmult.h"
#include <immintrin.h>
#include "../radix25/radix25_mul.h"
#include "../radix25/radix25_sqr.h"
#include "../radix25/radix25_addsub.h"

using namespace std;

TEST_CASE( "Radix25 - conversion", "[Radix25]" )
{
  // printf("Q:\t\t%016lx \n%016lx \n%016lx \n%016lx \n%016lx \n%016lx \n%016lx \n%016lx \n%016lx \n%016lx\n",
  //             2*CURVE_Q.limbs[0], 2*CURVE_Q.limbs[1], 2*CURVE_Q.limbs[2], 2*CURVE_Q.limbs[3], 2*CURVE_Q.limbs[4],
  //             2*CURVE_Q.limbs[5], 2*CURVE_Q.limbs[6], 2*CURVE_Q.limbs[7], 2*CURVE_Q.limbs[8], 2*CURVE_Q.limbs[9]);
  SECTION("identity"){
		r25_t x(ulli_Q);
		REQUIRE(x == ZERO);
	}

  SECTION("zero"){
		limb_t ulli_x[] = {0x0, 0x0, 0x0, 0x0, 0x0};
		r25_t x(ulli_x);
    ulli c[4];
		radixtoulli(&x, &c[0]);

		REQUIRE(ulli_x[1] == c[0]);
		REQUIRE(ulli_x[2] == c[1]);
		REQUIRE(ulli_x[3] == c[2]);
		REQUIRE(ulli_x[4] == c[3]);
	}

	SECTION("small numbers"){
		limb_t ulli_x[] = {0x0, 0x0, 0x0, 0x0, 0xC77DE711E};

		r25_t x(ulli_x);

    ulli c[4];
    radixtoulli(&x, &c[0]);

		REQUIRE(ulli_x[1] == c[0]);
		REQUIRE(ulli_x[2] == c[1]);
		REQUIRE(ulli_x[3] == c[2]);
		REQUIRE(ulli_x[4] == c[3]);
	}

	SECTION("large numbers"){
		limb_t ulli_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};

		r25_t x(ulli_x);
    
    ulli c[4];
    radixtoulli(&x, &c[0]);

		REQUIRE(ulli_x[1] == c[0]);
		REQUIRE(ulli_x[2] == c[1]);
		REQUIRE(ulli_x[3] == c[2]);
		REQUIRE(ulli_x[4] == c[3]);
	}

	SECTION("large numbers with square"){
		limb_t ulli_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		ulli expected[] = {0x3e0dae68d87559d9, 0x6d583f3989782ebb, 0xdc1a212518888180, 0xc05f04242dbfea82};
		r25_t x(ulli_x);

		x = x.sqr();
    
    ulli c[4];
    radixtoulli(&x, &c[0]);

		REQUIRE(expected[0] == c[0]);
		REQUIRE(expected[1] == c[1]);
		REQUIRE(expected[2] == c[2]);
		REQUIRE(expected[3] == c[3]);
	}

	SECTION("large with mul")
	{
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x1E9F9DF516146CCA, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x78031505f132f09c, 0xc2575c55d0f365af, 0xed10b08c2a30da40, 0x3e67765a710d0ccc};
		r25_t x(temp_x);
		r25_t y(temp_y);
		// r25_t ex(expected);

		r25_t z = x*y;
    
    ulli c[4];
    radixtoulli(&z, &c[0]);

		REQUIRE(expected[0] == c[0]);
		REQUIRE(expected[1] == c[1]);
		REQUIRE(expected[2] == c[2]);
		REQUIRE(expected[3] == c[3]);
	}
}

#if 0
TEST_CASE( "Radix25 - Addition", "[Radix25]" )
{
	SECTION("small numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_add1[] = {0x0, 0x0, 0x0, 0x0, 0x7FFFFFFFFD9DA5FF };
	    limb_t temp_add2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0xBD3F3BEA2AEA8BFE};
	    r25_t x(temp_add1);
	    r25_t y(temp_add2);
	    r25_t ex(expected);
	    r25_t z = x+y;

	    REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_add1[] = {0x0, 0x0, 0x0, 0xB1A2BC4AC4647A, 0x76B0177693348A1F };
	    limb_t temp_add2[] =  {0x0, 0x0, 0x0, 0x15C2620766D3, 0x1B300382F71DE6FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0xB1B87EACCBCB4D, 0x91E01AF98A52711E};
	    r25_t x(temp_add1);
	    r25_t y(temp_add2);
	    r25_t ex(expected);
	    r25_t z = x+y;
	 
	    REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{ 
	     	//values before normalization -> constructor will call carry and normalize
	    limb_t temp_add1[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E };
	    limb_t temp_add2[] =  {0x0, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};
	    limb_t expected[] = {0x0, 0x4473F8C, 0x08AACFC6D9754A98, 0x95E2A83A2F628698, 0xA09B30AB79872291};
	    r25_t x(temp_add1);
	    r25_t y(temp_add2);
	    r25_t ex(expected);
	    r25_t z = x+y;
	   
	    REQUIRE(z == ex);
	}

  SECTION("wrap around")
	{ 
	     	//values before normalization -> constructor will call carry and normalize
	    limb_t temp_add2[] =  {0x0, 0x0, 0x0, 0x0, 0x14};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x14};
	    r25_t x(ulli_Q);
      // printf("x:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
      //         x.limbs[0], x.limbs[1], x.limbs[2], x.limbs[3], x.limbs[4],
      //         x.limbs[5], x.limbs[6], x.limbs[7], x.limbs[8], x.limbs[9]);
      
	    r25_t y(temp_add2);
	    r25_t ex(expected);
	    r25_t z = x+y;
	   
	    REQUIRE(z == ex);
	}

}

TEST_CASE( "Radix25 - SUBTRACTION", "[Radix25]" ) 
{
    SECTION("neutral element")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x1 };
	    limb_t temp_sub2[] =  {0x0, 0x0, 0x0, 0x0, 0x0};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x1};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
	    r25_t ex(expected);
	    r25_t z = x-y;

	    REQUIRE(z == ex);
	}

	SECTION("small numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x7FFFFFFFFD9DA5FF };
	    limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x42C0C415D050C000};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
	    r25_t ex(expected);
	    r25_t z = x-y;
	    REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0xB1A2BC4AC4647A, 0x76B0177693348A1F };
	    limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x15C2620766D3, 0x1B300382F71DE6FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0xB18CF9E8BCFDA7, 0x5B8013F39C16A320};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
	    r25_t ex(expected);
	    r25_t z = x-y;
	    REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{ 
	     	//values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
	    limb_t temp_sub2[] =  {0x0000000000000, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};
	    limb_t expected[] = {0x0, 0x7FFFFFFFFBE59D62, 0x087C9658A9E2D3EE, 0x345E2D914FDA7343, 0xEFD64E6D7635BF98};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
	    r25_t ex(expected);
	    r25_t z = x-y;
	    REQUIRE(z == ex);
	}

    SECTION("negative numbers")
	{ 
	     	//values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
	    limb_t temp_sub2[] =  {0x0000000000000, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};
	    limb_t expected[] = {0x0, 0x41a629d, 0xf78369a7561d2c11, 0xcba1d26eb0258cbc, 0x1029b19289ca4055};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
	    r25_t ex(expected);
	    r25_t z = y-x;
	    REQUIRE(z == ex);
	}
}

#endif
TEST_CASE( "Radix25 - MULTIPLICATION", "[Radix25]" ) 
{
	SECTION("SANITY")
	{
		u_int64_t foo = 0x1;
		u_int64_t foo2 = 0xFFFFFFFFFFFFFFFF;
		u_int64_t bar = 0;
		u_int64_t c = _addc64(foo, foo2, &bar, 0);
		REQUIRE(c == 1);
	}

	SECTION("small numbers")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x2 };
		limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x04000000000000};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x8000000000000};
		r25_t x(temp_sub1);
		r25_t y(temp_sub2);
		r25_t ex(expected);
		
		r25_t z = x*y;
		REQUIRE(z == ex);
	}

	SECTION("aligned numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x1000, 0x0000000002000000, 0x0000004000000000, 0x0008000000000001};
	    limb_t temp_sub2[] = {0x0, 0x1000, 0x0000000002000000, 0x0000004000000000, 0x0008000000000001};
		// x = 1 + 2**51 + 2**102 + 2**153 + 2**204
		limb_t temp_exp[] = {0x0, 0x1000, 0x0000000002000000, 0x0000004000000000, 0x0008000000000001};
	    limb_t expected[] = {0x0, 0x5000, 0x000000002E000000, 0x00000A4000000000, 0x01D800000000004D};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
		r25_t xexp(temp_exp);
	    r25_t ex(expected);

	    r25_t z = x*y;

	    REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x7FFFFFFFFD9DA5FF};
	    limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x1E9F9DF516146CCA, 0x43B07EECC0397401};
	    r25_t x(temp_sub1);
	    r25_t y(temp_sub2);
	    r25_t ex(expected);
	    r25_t z = x*y;
		
	    REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x1E9F9DF516146CCA, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x78031505f132f09c, 0xc2575c55d0f365af, 0xed10b08c2a30da40, 0x3e67765a710d0ccc};
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y;
		// printf("result:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
		// 				z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4],
		// 				z.limbs[5], z.limbs[6], z.limbs[7], z.limbs[8], z.limbs[9]);
		// printf("expected:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
		// 				ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4],
		// 				ex.limbs[5], ex.limbs[6], ex.limbs[7], ex.limbs[8], ex.limbs[9]);
		

		REQUIRE(z == ex);
	}

	SECTION("First limb")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x11475c82c0eda21f,  0x146c737d3caa82e2};
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Second limb")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x4838BF8C77DE711E, 0x0};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x0, 0x11475c82c0eda21f, 0x0146c737d3caa82e2, 0x00000000000000000};
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Third limb")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x4838BF8C77DE711E, 0x0, 0x0};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x11475c82c0eda21f, 0x146c737d3caa82e2, 0x0000000000000000, 0x0000000000000000};
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Fourth limb")
	{
		limb_t temp_x[] = {0x0, 0x4838BF8C77DE711E, 0x0, 0x0, 0x0};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x146c737d3caa82e2, 0x0000000000000000, 0x0000000000000002, 0x9097bb68a346109a};
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y;
		// printf("result:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
		// 				z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4],
		// 				z.limbs[5], z.limbs[6], z.limbs[7], z.limbs[8], z.limbs[9]);
		// printf("expected:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
		// 				ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4],
		// 				ex.limbs[5], ex.limbs[6], ex.limbs[7], ex.limbs[8], ex.limbs[9]);
		
		REQUIRE(z == ex);
	}

	SECTION("Identity")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x1};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x4838BF8C77DE711E};
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y;

		REQUIRE(z == ex);
	}
}

TEST_CASE( "Radix25 - SQUARING", "[Radix25]" ) 
{
	SECTION("small numbers")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x2 };
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x4};
		r25_t x(temp_x);
		r25_t ex(expected);

		r25_t z = x.sqr();

		r25_t te(ulli_ZERO);
		// _r25_sqr(x.limbs, te.limbs);
		// printf("blbl:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
		// 				te.limbs[0], te.limbs[1], te.limbs[2], te.limbs[3], te.limbs[4],
		// 				te.limbs[5], te.limbs[6], te.limbs[7], te.limbs[8], te.limbs[9]);
		
		REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
		limb_t temp_x[] = {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0ea7322c24deb88d, 0x041e2021040a3401};
		r25_t x(temp_x);
		r25_t ex(expected);

		r25_t z = x.sqr();
		REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t expected[] = {0x0, 0x3e0dae68d87559d9, 0x6d583f3989782ebb, 0xdc1a212518888180, 0xc05f04242dbfea82};
		r25_t x(temp_x);
		r25_t ex(expected);

		r25_t z = x.sqr();

		REQUIRE(z == ex);
	}

	SECTION("very large numbers")
	{
		limb_t expected[] = {0x0, 0x7fffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffec};
		r25_t x = CURVE_I;
		r25_t ex(expected);

		r25_t z = x.sqr();

		REQUIRE(z == ex);
	}

}

TEST_CASE( "Radix25 - INVERSE", "[Radix25]" ) 
{
// Assumes multiplication works
  SECTION("small numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x2};
    limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x0};
		limb_t exp_inv[] = {0x0, 0x3fffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffffffffff7};
    
    r25_t x(temp_x);
    r25_t y(temp_y);
		r25_t invy(exp_inv);

    y = x.inv(); // y == x^{-1}
    
    REQUIRE(y == invy);
  }
  
  SECTION("medium numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_x[] = {0x0, 0x0, 0x0, 0xB1A2BC4AC4647A, 0x76B0177693348A1F};
    limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x0};
		limb_t exp_inv[] = {0x0, 0x2fcb000fc38ef656, 0x8b606c29d512053a, 0xf12318a9aa3102e6, 0x9a097c2e37edea77};
    
    r25_t x(temp_x);
    r25_t y(temp_y);
		r25_t invy(exp_inv);
    
    y = x.inv(); // y == x^{-1}
    
    REQUIRE(y == invy);
  }
  
  SECTION("large numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
    limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x0};
    limb_t expected[] = {0x0, 0x6b2c6ce8f04949ac, 0x6bb0f89dbbd5c91d, 0x7057401e6e09dc38, 0xb14cde306209ac72};
    
    r25_t x(temp_x);
    r25_t y(temp_y);
    r25_t ex(expected);
    
    y = x.inv(); // y == x^{-1}
    
    
    REQUIRE(y == ex);
  }
}

TEST_CASE( "Radix51 - DIVISION", "[Radix51]" )
{
	// FOR THE GLORY OF CYCLES
	// SECTION("identity division")
	// {
	// 	//values before normalization -> constructor will call carry and normalize
	// 	limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x2};
	// 	limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x1};
	// 	// _r25_inv(&temp_x[0], &temp_y[0]);
	// 	// _r25_mul(&temp_x[0], &temp_y[0], &temp_z[0]);
		
	// 	r25_t x(temp_x);
	// 	// r25_t z(temp_z);
	// 	r25_t ex(expected);

	// 		r25_t z = x/x; // neutral element
	// 		// printf("result:\t\t%016lx  %016lx  %016lx  %016lx  %016lx\n", z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4]);
	// 		// printf("expected:\t%016lx  %016lx  %016lx  %016lx  %016lx\n", ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4]);

	// 	REQUIRE(z == ex);
	// }

	SECTION("medium division")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x0c691b46e8200b53};
		limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x7057401e6e09dc38};
		limb_t expected[] = {0x0, 0x38f7bcfc6091d2d8, 0x48853e72b223d3e0, 0xb706a4f7014f5be0, 0x3b90fbbbca33c30c};
		// _r25_inv(&temp_x[0], &temp_y[0]);
		// _r25_mul(&temp_x[0], &temp_y[0], &temp_z[0]);
		
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y.inv(); // same as x/y
		REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t temp_y[] = {0x0, 0x6b2c6ce8f04949ac, 0x6bb0f89dbbd5c91d, 0x7057401e6e09dc38, 0xb14cde306209ac72};
		limb_t expected[] = {0x0, 0x3e0dae68d87559d9, 0x6d583f3989782ebb, 0xdc1a212518888180, 0xc05f04242dbfea82};
		
		r25_t x(temp_x);
		r25_t y(temp_y);
		r25_t ex(expected);

		r25_t z = x*y.inv(); // same as x/y
		// printf("result:\t\t%016lx  %016lx  %016lx  %016lx  %016lx\n", z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4]);
		// printf("expected:\t%016lx  %016lx  %016lx  %016lx  %016lx\n", ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4]);

		REQUIRE(z == ex);
	}
}

// TEST_CASE( "Radix25 - Shift left", "[Radix25]" )
// {
//   SECTION("small numbers")
//   {
//     //values before normalization -> constructor will call carry and normalize
//     limb_t temp_shl1[] = {0x0, 0x0, 0x0, 0x0, 0x1};
//     limb_t temp_shl2[] = {0x0, 0x0, 0x0, 0x0, 0x8000000000000000};
//     r25_t x(temp_shl1);
//     r25_t y(temp_shl2);
//    
//    
//     limb_t expected_x0[] = {0x0, 0x0, 0x0, 0x0, 0x1};
//     limb_t expected_x1[] = {0x0, 0x0, 0x0, 0x0, 0x2};
//     limb_t expected_y1[] = {0x0, 0x0, 0x0, 0x1, 0x0};
//     r25_t ex0(expected_x0);
//     r25_t ex1(expected_x1);
//     r25_t ex2(expected_y1);
//    
//     REQUIRE((x << 0) == ex0);
//     REQUIRE((x << 1) == ex1);
//     REQUIRE((y << 1) == ex2);
//   }
//
// 	SECTION("all the bits")
// 	{
// 		limb_t temp1[] = {0,0,0,0,1};
// 		r25_t x(temp1);
// 		for (int i = 0; i < 64; i++)
// 		{
// 			printf("Shift %d\n", i);
// 			r25_t z = x << i;
// 			limb_t expected[5] = {0,0,0,0,0};
// 			expected[4] = 1 << i;
// 			r25_t ex(expected);
// 			printf("result:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
// 						z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4],
// 						z.limbs[5], z.limbs[6], z.limbs[7], z.limbs[8], z.limbs[9]);
// 			printf("expected:\n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx \n\t%016lx\n",
// 						ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4],
// 						ex.limbs[5], ex.limbs[6], ex.limbs[7], ex.limbs[8], ex.limbs[9]);
//
// 			REQUIRE(z == ex);
// 		}
//
//
// 	}
// }

TEST_CASE("Scalar Multiplication", "[ScalarMult]")
{
	SECTION("identity")
    {
        ulli a[5];
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;
        a[3] = 0;
        a[4] = 1;
        Radix25 ar(a);

        ulli A0EXPECTED[5];
        A0EXPECTED[0] = 0;
        A0EXPECTED[1] = 0x216936d3cd6e53fe;
        A0EXPECTED[2] = 0xc0a4e231fdd6dc5c;
        A0EXPECTED[3] = 0x692cc7609525a7b2;
        A0EXPECTED[4] = 0xc9562d608f25d51a;
        Radix25 a0exp(A0EXPECTED);

        ulli A1EXPECTED[5];
        A1EXPECTED[0] = 0;
        A1EXPECTED[1] = 0x6666666666666666;
        A1EXPECTED[2] = 0x6666666666666666;
        A1EXPECTED[3] = 0x6666666666666666;
        A1EXPECTED[4] = 0x6666666666666658;
        Radix25 a1exp(A1EXPECTED);
        r25_t A[2] = {Radix25(ulli_ZERO), Radix25(ulli_ZERO)};
        scalar_mult(R25_CURVE_B2, ar, A);
        REQUIRE(A[0] == a0exp);
        REQUIRE(A[1] == a1exp);

    }

	
	SECTION("test1")
    {
        ulli a[5];
        a[0] = 0;
        a[1] = 0x4fe94d9006f020a5;
        a[2] = 0xa3c080d96827fffd;
        a[3] = 0x3c010ac0f12e7a42;
        a[4] = 0xcb33284f86837c30;
        Radix25 ar(a);

        ulli A0EXPECTED[5];
        A0EXPECTED[0] = 0;
        A0EXPECTED[1] = 0x55d0e09a2b9d3429;
        A0EXPECTED[2] = 0x2297e08d60d0f620;
        A0EXPECTED[3] = 0xc513d47253187c24;
        A0EXPECTED[4] = 0xb12786bd777645ce;
        Radix25 a0exp(A0EXPECTED);
        
        ulli A1EXPECTED[5];
        A1EXPECTED[0] = 0;
        A1EXPECTED[1] = 0x1a5107f7681a02af;
        A1EXPECTED[2] = 0x2523a6daf372e10e;
        A1EXPECTED[3] = 0x3a0764c9d3fe4bd5;
        A1EXPECTED[4] = 0xb70ab18201985ad7;
        Radix25 a1exp(A1EXPECTED);
        r25_t A[2] = {Radix25(ulli_ZERO), Radix25(ulli_ZERO)};
        scalar_mult(R25_CURVE_B2, ar, A);
        REQUIRE(A[0].limbs[0] == a0exp.limbs[0]);
        REQUIRE(A[0] == a0exp);
        REQUIRE(A[1] == a1exp);

    }
	

	
	SECTION("test2")
    {
        ulli a[5];
        a[0] = 0;
        a[1] = 0x512e502eb0249a25;
        a[2] = 0x5e1c827f3b6b6c7f;
        a[3] = 0x0a79f4ca8575a915;
        a[4] = 0x28d58258d79ebd68;
        Radix25 ar(a);

        ulli A0EXPECTED[5];
        A0EXPECTED[0] = 0;
        A0EXPECTED[1] = 0x74ad28205b4f384b;
        A0EXPECTED[2] = 0xc0813e6585864e52;
        A0EXPECTED[3] = 0x8085f91fb6a5096f;
        A0EXPECTED[4] = 0x244ae01e57de43ae;
        Radix25 a0exp(A0EXPECTED);
        
        ulli A1EXPECTED[5];
        A1EXPECTED[0] = 0;
        A1EXPECTED[1] = 0x0c66f42af155cdc0;
        A1EXPECTED[2] = 0x8c96c42ecf2c989c;
        A1EXPECTED[3] = 0xbc7e1b4da70ab792;
        A1EXPECTED[4] = 0x5a8943e8c317403d;
        Radix25 a1exp(A1EXPECTED);
        r25_t A[2] = {Radix25(ulli_ZERO), Radix25(ulli_ZERO)};
        scalar_mult(R25_CURVE_B2, ar, A);
        REQUIRE(A[0].limbs[0] == a0exp.limbs[0]);

        REQUIRE(A[0] == a0exp);
        REQUIRE(A[1] == a1exp);
    }
}


TEST_CASE("Branchless Scalar Multiplication", "[ScalarMult]")
{
	// FOR THE GLORY OF CYCLES
    // SECTION("identity")
    // {
    //     ulli a[5];
    //     a[0] = 0;
    //     a[1] = 0;
    //     a[2] = 0;
    //     a[3] = 0;
    //     a[4] = 1;
    //     Radix25 ar(a);

    //     ulli A0EXPECTED[5];
    //     A0EXPECTED[0] = 0;
    //     A0EXPECTED[1] = 0x216936d3cd6e53fe;
    //     A0EXPECTED[2] = 0xc0a4e231fdd6dc5c;
    //     A0EXPECTED[3] = 0x692cc7609525a7b2;
    //     A0EXPECTED[4] = 0xc9562d608f25d51a;
    //     Radix25 a0exp(A0EXPECTED);

    //     ulli A1EXPECTED[5];
    //     A1EXPECTED[0] = 0;
    //     A1EXPECTED[1] = 0x6666666666666666;
    //     A1EXPECTED[2] = 0x6666666666666666;
    //     A1EXPECTED[3] = 0x6666666666666666;
    //     A1EXPECTED[4] = 0x6666666666666658;
    //     Radix25 a1exp(A1EXPECTED);
    //     r25_t A[2] = {Radix25(ulli_ZERO), Radix25(ulli_ZERO)};
    //     scalar_mult_branchless(R25_CURVE_B2, ar, A);
    //     REQUIRE(A[0] == a0exp);
    //     REQUIRE(A[1] == a1exp);
    // }

	SECTION("test1")
        {
			ulli *a = (ulli*) calloc(5, sizeof(ulli));
			a[0] = 0;
			a[1] = 0x4fe94d9006f020a5;
			a[2] = 0xa3c080d96827fffd;
			a[3] = 0x3c010ac0f12e7a42;
			a[4] = 0xcb33284f86837c30;
			Radix25 ar(a);

			ulli A0EXPECTED[5];
			A0EXPECTED[0] = 0;
			A0EXPECTED[1] = 0x55d0e09a2b9d3429;
			A0EXPECTED[2] = 0x2297e08d60d0f620;
			A0EXPECTED[3] = 0xc513d47253187c24;
			A0EXPECTED[4] = 0xb12786bd777645ce;
			Radix25 a0exp((ulli*)&A0EXPECTED);
			
			ulli A1EXPECTED[5];
			A1EXPECTED[0] = 0;
			A1EXPECTED[1] = 0x1a5107f7681a02af;
			A1EXPECTED[2] = 0x2523a6daf372e10e;
			A1EXPECTED[3] = 0x3a0764c9d3fe4bd5;
			A1EXPECTED[4] = 0xb70ab18201985ad7;
			Radix25 a1exp((ulli*)&A1EXPECTED);
			r25_t A[2] = {Radix25(ulli_ZERO), Radix25(ulli_ZERO)};
			scalar_mult_branchless(R25_CURVE_B2, ar, A);
			REQUIRE(A[0].limbs[0] == a0exp.limbs[0]);
			REQUIRE(A[0].limbs[1] == a0exp.limbs[1]);
			REQUIRE(A[0].limbs[2] == a0exp.limbs[2]);
			REQUIRE(A[0].limbs[3] == a0exp.limbs[3]);
			REQUIRE(A[0].limbs[4] == a0exp.limbs[4]);

			REQUIRE(A[1].limbs[0] == a1exp.limbs[0]);
			REQUIRE(A[1].limbs[1] == a1exp.limbs[1]);
			REQUIRE(A[1].limbs[2] == a1exp.limbs[2]);
			REQUIRE(A[1].limbs[3] == a1exp.limbs[3]);
			REQUIRE(A[1].limbs[4] == a1exp.limbs[4]);
        }
}

TEST_CASE( "Radix25 - MULTIPLICATIONx4", "[Radix25]" ) 
{
	SECTION("test1")
	{
	    limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
	    limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
	    limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
 	    limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};
	    limb_t temp5[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
	    limb_t temp6[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
	    limb_t temp7[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
 	    limb_t temp8[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

	    limb_t temp0[] = {0, 0, 0, 0, 0};

	    r25_t x1(temp1);
	    r25_t x2(temp2);
	    r25_t x3(temp3);
	    r25_t x4(temp4);
	    r25_t y1(temp5);
	    r25_t y2(temp6);
	    r25_t y3(temp7);
	    r25_t y4(temp8);

	    r25_t z1_ex = x1 * y1;
	    r25_t z2_ex = x2 * y2;
	    r25_t z3_ex = x3 * y3;
	    r25_t z4_ex = x4 * y4;


	    __m256i converted_x [10];
	    __m256i converted_y [10];

	    
	    for(int i=0;i<10;i++)
	    {
		converted_x[i] = _mm256_set_epi64x(x4.limbs[i], x3.limbs[i], x2.limbs[i], x1.limbs[i]);
		converted_y[i] = _mm256_set_epi64x(y4.limbs[i], y3.limbs[i], y2.limbs[i], y1.limbs[i]);
	    }

	    __m256i converted_z [10];

	    _v_r25_mul(converted_x, converted_y, converted_z);

	    unsigned long long int calc[10][4] ={
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0},
		    {0, 0, 0, 0}};

	    for(int i=0;i<10;i++)
	    {
	    	_mm256_storeu_si256((__m256i*)&calc[i], converted_z[i]);
	    }
	    
	    r25_t z1(temp0);
	    r25_t z2(temp0);
	    r25_t z3(temp0);
	    r25_t z4(temp0);

	    for(int i=0;i<10;i++)
	    {
		    z1.limbs[i] = calc[i][0];
		    z2.limbs[i] = calc[i][1];
		    z3.limbs[i] = calc[i][2];
		    z4.limbs[i] = calc[i][3];
	    }

	    REQUIRE(z1 == z1_ex);
	    REQUIRE(z2 == z2_ex);
	    REQUIRE(z3 == z3_ex);
	    REQUIRE(z4 == z4_ex);

	}

}

TEST_CASE( "Radix25 - SQRx4", "[Radix25]" )
{

  SECTION("test1")
  {
      limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
      limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
      limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
      limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

      limb_t temp0[] = {0, 0, 0, 0, 0};

      r25_t x1(temp1);
      r25_t x2(temp2);
      r25_t x3(temp3);
      r25_t x4(temp4);

      r25_t z1_ex = x1 * x1;
      r25_t z2_ex = x2 * x2;
      r25_t z3_ex = x3 * x3;
      r25_t z4_ex = x4 * x4;


      __m256i converted_x [10];

      
      for(int i = 0 ; i < 10 ; i++) {
        converted_x[i] = _mm256_set_epi64x(x4.limbs[i], x3.limbs[i], x2.limbs[i], x1.limbs[i]);
      }

      __m256i converted_z [10];

      _v_r25_sqr(converted_x, converted_z);

      unsigned long long int calc[10][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}};

      for(int i = 0 ; i < 10 ; i++) {
        _mm256_storeu_si256((__m256i*)&calc[i], converted_z[i]);
      }
      
      r25_t z1(temp0);
      r25_t z2(temp0);
      r25_t z3(temp0);
      r25_t z4(temp0);

      for(int i = 0 ; i < 10 ; i++) {
        z1.limbs[i] = calc[i][0];
        z2.limbs[i] = calc[i][1];
        z3.limbs[i] = calc[i][2];
        z4.limbs[i] = calc[i][3];
      }

      REQUIRE(z1 == z1_ex);
      REQUIRE(z2 == z2_ex);
      REQUIRE(z3 == z3_ex);
      REQUIRE(z4 == z4_ex);

  }

}

TEST_CASE( "Radix25 - ADDSUBx4", "[Radix25]" )
{
  SECTION("addx4")
  {
      
      limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
      limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
      limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
      limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

      limb_t temp0[] = {0, 0, 0, 0, 0};

      r25_t x1(temp1);
      r25_t x2(temp2);
      r25_t x3(temp3);
      r25_t x4(temp4);
      r25_t y1(temp4);
      r25_t y2(temp3);
      r25_t y3(temp2);
      r25_t y4(temp1);

      r25_t z1_ex = x1 + y1;
      r25_t z2_ex = x2 + y2;
      r25_t z3_ex = x3 + y3;
      r25_t z4_ex = x4 + y4;


      __m256i converted_x[10];
      __m256i converted_y[10];

      
      for(int i = 0 ; i < 10 ; i++) {
        converted_x[i] = _mm256_set_epi64x(x4.limbs[i], x3.limbs[i], x2.limbs[i], x1.limbs[i]);
        converted_y[i] = _mm256_set_epi64x(y4.limbs[i], y3.limbs[i], y2.limbs[i], y1.limbs[i]);
      }

      __m256i converted_z[10];

      _v_r25_add(converted_x, converted_y, converted_z);

      unsigned long long int calc[10][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}};

      for(int i = 0 ; i < 10 ; i++) {
        _mm256_storeu_si256((__m256i*)&calc[i], converted_z[i]);
      }
      
      r25_t z1(temp0);
      r25_t z2(temp0);
      r25_t z3(temp0);
      r25_t z4(temp0);

      for(int i = 0 ; i < 10 ; i++) {
        z1.limbs[i] = calc[i][0];
        z2.limbs[i] = calc[i][1];
        z3.limbs[i] = calc[i][2];
        z4.limbs[i] = calc[i][3];
      }

      REQUIRE(z1 == z1_ex);
      REQUIRE(z2 == z2_ex);
      REQUIRE(z3 == z3_ex);
      REQUIRE(z4 == z4_ex);

  }
  
    SECTION("sub neutral elementx4")
    {
      limb_t temp1[] = {0x0, 0x0, 0x0, 0x0, 0x1};
      limb_t temp2[] = {0x0, 0x0, 0x0, 0x0, 0x0};
      limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
      limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

      limb_t temp0[] = {0, 0, 0, 0, 0};

      r25_t x1(temp1);
      r25_t x2(temp2);
      r25_t x3(temp3);
      r25_t x4(temp4);
      r25_t y1(temp0);
      r25_t y2(temp0);
      r25_t y3(temp0);
      r25_t y4(temp0);

      r25_t z1_ex = x1;
      r25_t z2_ex = x2;
      r25_t z3_ex = x3;
      r25_t z4_ex = x4;


      __m256i converted_x[10];
      __m256i converted_y[10];

      
      for(int i = 0 ; i < 10 ; i++) {
        converted_x[i] = _mm256_set_epi64x(x4.limbs[i], x3.limbs[i], x2.limbs[i], x1.limbs[i]);
        converted_y[i] = _mm256_set_epi64x(y4.limbs[i], y3.limbs[i], y2.limbs[i], y1.limbs[i]);
      }

      __m256i converted_z[10];

      _v_r25_sub(converted_x, converted_y, converted_z);

      unsigned long long int calc[10][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}};

      for(int i = 0 ; i < 10 ; i++) {
        _mm256_storeu_si256((__m256i*)&calc[i], converted_z[i]);
      }
      
      r25_t z1(temp0);
      r25_t z2(temp0);
      r25_t z3(temp0);
      r25_t z4(temp0);

      for(int i = 0 ; i < 10 ; i++) {
        z1.limbs[i] = calc[i][0];
        z2.limbs[i] = calc[i][1];
        z3.limbs[i] = calc[i][2];
        z4.limbs[i] = calc[i][3];
      }

      REQUIRE(z1 == z1_ex);
      REQUIRE(z2 == z2_ex);
      REQUIRE(z3 == z3_ex);
      REQUIRE(z4 == z4_ex);

    }

  SECTION("subx4")
  {
      
      limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
      limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
      limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
      limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

      limb_t temp0[] = {0, 0, 0, 0, 0};

      r25_t x1(temp1);
      r25_t x2(temp2);
      r25_t x3(temp3);
      r25_t x4(temp4);
      r25_t y1(temp4);
      r25_t y2(temp3);
      r25_t y3(temp2);
      r25_t y4(temp1);

      r25_t z1_ex = x1 - y1;
      r25_t z2_ex = x2 - y2;
      r25_t z3_ex = x3 - y3;
      r25_t z4_ex = x4 - y4;


      __m256i converted_x[10];
      __m256i converted_y[10];

      
      for(int i = 0 ; i < 10 ; i++) {
        converted_x[i] = _mm256_set_epi64x(x4.limbs[i], x3.limbs[i], x2.limbs[i], x1.limbs[i]);
        converted_y[i] = _mm256_set_epi64x(y4.limbs[i], y3.limbs[i], y2.limbs[i], y1.limbs[i]);
      }

      __m256i converted_z[10];

      _v_r25_sub(converted_x, converted_y, converted_z);

      unsigned long long int calc[10][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}};

      for(int i = 0 ; i < 10 ; i++) {
        _mm256_storeu_si256((__m256i*)&calc[i], converted_z[i]);
      }
      
      r25_t z1(temp0);
      r25_t z2(temp0);
      r25_t z3(temp0);
      r25_t z4(temp0);

      for(int i = 0 ; i < 10 ; i++) {
        z1.limbs[i] = calc[i][0];
        z2.limbs[i] = calc[i][1];
        z3.limbs[i] = calc[i][2];
        z4.limbs[i] = calc[i][3];
      }

      REQUIRE(z1 == z1_ex);
      REQUIRE(z2 == z2_ex);
      REQUIRE(z3 == z3_ex);
      REQUIRE(z4 == z4_ex);

  }
  
  SECTION("negx4")
  {
      
      limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
      limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
      limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
      limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

      limb_t temp0[] = {0, 0, 0, 0, 0};

      r25_t x1(temp1);
      r25_t x2(temp2);
      r25_t x3(temp3);
      r25_t x4(temp4);
      r25_t x0(temp0);

      r25_t z1_ex = x0 - x1;
      r25_t z2_ex = x0 - x2;
      r25_t z3_ex = x0 - x3;
      r25_t z4_ex = x0 - x4;


      __m256i converted_x[10];

      
      for(int i = 0 ; i < 10 ; i++) {
        converted_x[i] = _mm256_set_epi64x(x4.limbs[i], x3.limbs[i], x2.limbs[i], x1.limbs[i]);
      }

      __m256i converted_z[10];

      _v_r25_neg(converted_x, converted_z);

      unsigned long long int calc[10][4] = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}};

      for(int i = 0 ; i < 10 ; i++) {
        _mm256_storeu_si256((__m256i*)&calc[i], converted_z[i]);
      }
      
      r25_t z1(temp0);
      r25_t z2(temp0);
      r25_t z3(temp0);
      r25_t z4(temp0);

      for(int i = 0 ; i < 10 ; i++) {
        z1.limbs[i] = calc[i][0];
        z2.limbs[i] = calc[i][1];
        z3.limbs[i] = calc[i][2];
        z4.limbs[i] = calc[i][3];
      }

      REQUIRE(z1 == z1_ex);
      REQUIRE(z2 == z2_ex);
      REQUIRE(z3 == z3_ex);
      REQUIRE(z4 == z4_ex);

  }

}

TEST_CASE( "Radix25 - Scalar_multx4", "[Radix25]" )
{
	SECTION("identity")
    {
        ulli a[5];
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;
        a[3] = 0;
        a[4] = 1;
        Radix25 ar(a);

        ulli A0EXPECTED[5];
        A0EXPECTED[0] = 0;
        A0EXPECTED[1] = 0x216936d3cd6e53fe;
        A0EXPECTED[2] = 0xc0a4e231fdd6dc5c;
        A0EXPECTED[3] = 0x692cc7609525a7b2;
        A0EXPECTED[4] = 0xc9562d608f25d51a;
        Radix25 a0exp((ulli*)&A0EXPECTED);

        ulli A1EXPECTED[5];
        A1EXPECTED[0] = 0;
        A1EXPECTED[1] = 0x6666666666666666;
        A1EXPECTED[2] = 0x6666666666666666;
        A1EXPECTED[3] = 0x6666666666666666;
        A1EXPECTED[4] = 0x6666666666666658;
        Radix25 a1exp((ulli*)&A1EXPECTED);

        r25_t sink;
        r25_t xibalba[2];
        r25_t E[4];
		E[0] = ar; E[1] = ar; E[2] = ar; E[3] = ar;
        r25_t P[4][2];
		P[0][0] = R25_CURVE_B2[0]; P[1][0] = R25_ZERO; P[2][0] = R25_ZERO; P[3][0] = R25_ZERO;
		P[0][1] = R25_CURVE_B2[1]; P[1][1] = R25_ZERO; P[2][1] = R25_ZERO; P[3][1] = R25_ZERO;

        r25_t res0[2] = {Radix25(ulli_ZERO), Radix25(ulli_ZERO)};
        v_scalar_mult(P, E, res0, xibalba, xibalba, xibalba);


        REQUIRE(res0[0] == a0exp);
        REQUIRE(res0[1] == a1exp);
    }

	SECTION("Test1")
	{
		// const numbers
		limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
		limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
		limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
		limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};
		limb_t temp5[] = {0x0, 0x1d2f4909d7fa7539, 0x2053473460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
		limb_t temp6[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a241ee199b3eb, 0x5ce89754aa0cae48};
		limb_t temp7[] = {0x0, 0x41ac52533aed04c4, 0x32a675e52a3fc8e2, 0x517d3c1edeadbeef, 0x2a6598d95e75122d};
		limb_t temp8[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f77666e9e6cce4, 0x2262f8aca9c26644};
		limb_t temp9[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d6bb05d3e0a09b};
    	limb_t temp10[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa00ae48};
    	limb_t temp11[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e67118312, 0x2a6598d95eaa122d};
    	limb_t temp12[] = {0x0, 0x214c921a9c0f39cc, 0x22b5a094467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca98e6644};

		// setup all variables, inputs
		r25_t p0x(temp1);
		r25_t p0y(temp2);
		r25_t p1x(temp3);
		r25_t p1y(temp4);
		r25_t p2x(temp5);
		r25_t p2y(temp6);
		r25_t p3x(temp7);
		r25_t p3y(temp8);
		r25_t p0[2]; p0[0] = p0x; p0[1] = p0y;
		r25_t p1[2]; p1[0] = p1x; p1[1] = p1y;
		r25_t p2[2]; p2[0] = p2x; p2[1] = p2y;
		r25_t p3[2]; p3[0] = p3x; p3[1] = p3y;
		r25_t P[4][2];
		P[0][0] = p0x; P[1][0] = p1x; P[2][0] = p2x; P[3][0] = p3x;
		P[0][1] = p0y; P[1][1] = p1y; P[2][1] = p2y; P[3][1] = p3y;

		r25_t e0(temp9);
		r25_t e1(temp10);
		r25_t e2(temp11);
		r25_t e3(temp12);
		r25_t E[4];
		E[0] = e0; E[1] = e1; E[2] = e2; E[3] = e3;

		// setup expected results
		r25_t res0[2];
		r25_t res1[2];
		r25_t res2[2];
		r25_t res3[2];

		r25_t expected0[2];
		r25_t expected1[2];
		r25_t expected2[2];
		r25_t expected3[2];
		// expected results depend on previous (tested) scalarmult
		scalar_mult_branchless(p0, e0, expected0);
		scalar_mult(p1, e1, expected1);
		scalar_mult(p2, e2, expected2);
		scalar_mult(p3, e3, expected3);
		
		// actual testing
		v_scalar_mult(P, E, &res0[0], &res1[0], &res2[0], &res3[0]);

		// first point
		REQUIRE(res0[0] == expected0[0]);
		REQUIRE(res0[1] == expected0[1]);
		// second point
		REQUIRE(res1[0] == expected1[0]);
		REQUIRE(res1[1] == expected1[1]);
		// third point
		REQUIRE(res2[0] == expected2[0]);
		REQUIRE(res2[1] == expected2[1]);
		// fourth point
		REQUIRE(res3[0] == expected3[0]);
		REQUIRE(res3[1] == expected3[1]);

	}

}
