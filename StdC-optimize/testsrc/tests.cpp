#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "../../catch2/catch.hpp"
#include "../radix/radix51.h"
#include "../const.h"
// #include "../../util.h"
#include <iostream>
#include <stdio.h>
#include "../scalarmult/scalarmult.h"

using namespace std;

#if 0
TEST_CASE( "Radix51 - Addition", "[Radix51]" )
{
	SECTION("small numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_add1[] = {0x0, 0x0, 0x0, 0x0, 0x7FFFFFFFFD9DA5FF };
	    limb_t temp_add2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0xBD3F3BEA2AEA8BFE};
	    r51_t x(temp_add1);
	    r51_t y(temp_add2);
	    r51_t ex(expected);
	    r51_t z = x+y;

	    REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_add1[] = {0x0, 0x0, 0x0, 0xB1A2BC4AC4647A, 0x76B0177693348A1F };
	    limb_t temp_add2[] =  {0x0, 0x0, 0x0, 0x15C2620766D3, 0x1B300382F71DE6FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0xB1B87EACCBCB4D, 0x91E01AF98A52711E};
	    r51_t x(temp_add1);
	    r51_t y(temp_add2);
	    r51_t ex(expected);
	    r51_t z = x+y;
	 
	    REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{ 
	     	//values before normalization -> constructor will call carry and normalize
	    limb_t temp_add1[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E };
	    limb_t temp_add2[] =  {0x0, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};
	    limb_t expected[] = {0x0, 0x4473F8C, 0x08AACFC6D9754A98, 0x95E2A83A2F628698, 0xA09B30AB79872291};
	    r51_t x(temp_add1);
	    r51_t y(temp_add2);
	    r51_t ex(expected);
	    r51_t z = x+y;
	   
	    REQUIRE(z == ex);
	}

}

TEST_CASE( "Radix51 - SUBTRACTION", "[Radix51]" ) 
{
	SECTION("small numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x7FFFFFFFFD9DA5FF };
	    limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x42C0C415D050C000};
	    r51_t x(temp_sub1);
	    r51_t y(temp_sub2);
	    r51_t ex(expected);
	    r51_t z = x-y;

	    REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0xB1A2BC4AC4647A, 0x76B0177693348A1F };
	    limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x15C2620766D3, 0x1B300382F71DE6FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0xB18CF9E8BCFDA7, 0x5B8013F39C16A320};
	    r51_t x(temp_sub1);
	    r51_t y(temp_sub2);
	    r51_t ex(expected);
	    r51_t z = x-y;
	 
	    REQUIRE(z == ex);
	}

	/* TODO: this case actually never happens so not needed for now */
	SECTION("large numbers")
	{ 
	     	//values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
	    limb_t temp_sub2[] =  {0x0000000000000, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};
	    limb_t expected[] = {0x0, 0x7FFFFFFFFBE59D62, 0x087C9658A9E2D3EE, 0x345E2D914FDA7343, 0xEFD64E6D7635BF98};
	    r51_t x(temp_sub1);
	    r51_t y(temp_sub2);
	    r51_t ex(expected);
	    r51_t z = x-y;

	    REQUIRE(z == ex);
	}
}
#endif
TEST_CASE( "Radix51 - MULTIPLICATION", "[Radix51]" ) 
{
	SECTION("small numbers")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x2 };
		limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x04000000000000};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x8000000000000};
		r51_t x(temp_sub1);
		r51_t y(temp_sub2);
		r51_t ex(expected);
		
		r51_t z = x*y;
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
	    r51_t x(temp_sub1);
	    r51_t y(temp_sub2);
		r51_t xexp(temp_exp);
	    r51_t ex(expected);

	    r51_t z = x*y;

	    REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
	    //values before normalization -> constructor will call carry and normalize
	    limb_t temp_sub1[] = {0x0, 0x0, 0x0, 0x0, 0x7FFFFFFFFD9DA5FF};
	    limb_t temp_sub2[] =  {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
	    limb_t expected[] = {0x0, 0x0, 0x0, 0x1E9F9DF516146CCA, 0x43B07EECC0397401};
	    r51_t x(temp_sub1);
	    r51_t y(temp_sub2);
	    r51_t ex(expected);
	    r51_t z = x*y;
		
	    REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x1E9F9DF516146CCA, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x78031505f132f09c, 0xc2575c55d0f365af, 0xed10b08c2a30da40, 0x3e67765a710d0ccc};
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("First limb")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x11475c82c0eda21f,  0x146c737d3caa82e2};
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Second limb")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x4838BF8C77DE711E, 0x0};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x0, 0x11475c82c0eda21f, 0x0146c737d3caa82e2, 0x00000000000000000};
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Third limb")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x4838BF8C77DE711E, 0x0, 0x0};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x11475c82c0eda21f, 0x146c737d3caa82e2, 0x0000000000000000, 0x0000000000000000};
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Fourth limb")
	{
		limb_t temp_x[] = {0x0, 0x4838BF8C77DE711E, 0x0, 0x0, 0x0};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x146c737d3caa82e2, 0x0000000000000000, 0x0000000000000002, 0x9097bb68a346109a};
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y;

		REQUIRE(z == ex);
	}

	SECTION("Identity")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x4838BF8C77DE711E};
		limb_t temp_y[] =  {0x0, 0x0, 0x0, 0x0, 0x1};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x4838BF8C77DE711E};
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y;

		REQUIRE(z == ex);
	}

}

TEST_CASE( "Radix51 - SQUARING", "[Radix51]" ) 
{
	SECTION("small numbers")
	{
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x2 };
		limb_t temp_z[] =  {0x0, 0x0, 0x0, 0x0, 0x0};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x4};
		r51_t x(temp_x);
		r51_t ex(expected);

		r51_t z(temp_z);
		_r51_sqr(temp_x, temp_z);
		z.limbs[0] = temp_z[0];
		z.limbs[1] = temp_z[1];
		z.limbs[2] = temp_z[2];
		z.limbs[3] = temp_z[3];
		z.limbs[4] = temp_z[4];
		
		REQUIRE(z == ex);
	}

	SECTION("medium numbers")
	{
		limb_t temp_x[] = {0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0, 0x3D3F3BEA2D4CE5FF};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0ea7322c24deb88d, 0x041e2021040a3401};
		r51_t x(temp_x);
		r51_t ex(expected);

		r51_t z = x.sqr();
		REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t expected[] = {0x0, 0x3e0dae68d87559d9, 0x6d583f3989782ebb, 0xdc1a212518888180, 0xc05f04242dbfea82};
		r51_t x(temp_x);
		r51_t ex(expected);

		r51_t z = x.sqr();

		REQUIRE(z == ex);
	}

	SECTION("very large numbers")
	{
		limb_t expected[] = {0x0, 0x7fffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffec};
		r51_t x = CURVE_I;
		r51_t ex(expected);

		r51_t z = x.sqr();

		REQUIRE(z == ex);
	}

}

TEST_CASE( "Radix51 - Shift left", "[Radix51]" )
{
  SECTION("small numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_shl1[] = {0x0, 0x0, 0x0, 0x0, 0x1};
    limb_t temp_shl2[] =  {0x0, 0x0, 0x0, 0x0, 0x8000000000000000};
    r51_t x(temp_shl1);
    r51_t y(temp_shl2);
    
    
    limb_t expected_x0[] = {0x0, 0x0, 0x0, 0x0, 0x1};
    limb_t expected_x1[] = {0x0, 0x0, 0x0, 0x0, 0x2};
    limb_t expected_y1[] = {0x0, 0x0, 0x0, 0x1, 0x0};
    r51_t ex0(expected_x0);
    r51_t ex1(expected_x1);
    r51_t ex2(expected_y1);
    
    REQUIRE((x << 0) == ex0);
    REQUIRE((x << 1) == ex1);
    REQUIRE((y << 1) == ex2);
  }
}

TEST_CASE( "Radix51 - INVERSE", "[Radix51]" ) 
{
// Assumes multiplication works
  SECTION("small numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x2};
    limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x0};
		limb_t exp_inv[] = {0x0, 0x3fffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xfffffffffffffff7};
    
    r51_t x(temp_x);
    r51_t y(temp_y);
		r51_t invy(exp_inv);

    y = x.inv(); // y == x^{-1}
    
    REQUIRE(y == invy);
  }
  
  SECTION("medium numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_x[] = {0x0, 0x0, 0x0, 0xB1A2BC4AC4647A, 0x76B0177693348A1F};
    limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x0};
		limb_t exp_inv[] = {0x0, 0x2fcb000fc38ef656, 0x8b606c29d512053a, 0xf12318a9aa3102e6, 0x9a097c2e37edea77};
    
    r51_t x(temp_x);
    r51_t y(temp_y);
		r51_t invy(exp_inv);
    
    y = x.inv(); // y == x^{-1}
    
    REQUIRE(y == invy);
  }
  
  SECTION("large numbers")
  {
    //values before normalization -> constructor will call carry and normalize
    limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
    limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x0};
    limb_t expected[] = {0x0, 0x6b2c6ce8f04949ac, 0x6bb0f89dbbd5c91d, 0x7057401e6e09dc38, 0xb14cde306209ac72};
    
    r51_t x(temp_x);
    r51_t y(temp_y);
    r51_t ex(expected);
    
    y = x.inv(); // y == x^{-1}
    
    
    REQUIRE(y == ex);
  }
}

TEST_CASE( "Radix51 - DIVISION", "[Radix51]" )
{
	SECTION("identity division")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x2};
		limb_t expected[] = {0x0, 0x0, 0x0, 0x0, 0x1};
		// _r51_inv(&temp_x[0], &temp_y[0]);
		// _r51_mul(&temp_x[0], &temp_y[0], &temp_z[0]);
		
		r51_t x(temp_x);
		// r51_t z(temp_z);
		r51_t ex(expected);

			r51_t z = x/x; // neutral element
			// printf("result:\t\t%016lx  %016lx  %016lx  %016lx  %016lx\n", z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4]);
			// printf("expected:\t%016lx  %016lx  %016lx  %016lx  %016lx\n", ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4]);

		REQUIRE(z == ex);
	}

	SECTION("medium division")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_x[] = {0x0, 0x0, 0x0, 0x0, 0x0c691b46e8200b53};
		limb_t temp_y[] = {0x0, 0x0, 0x0, 0x0, 0x7057401e6e09dc38};
		limb_t expected[] = {0x0, 0x38f7bcfc6091d2d8, 0x48853e72b223d3e0, 0xb706a4f7014f5be0, 0x3b90fbbbca33c30c};
		// _r51_inv(&temp_x[0], &temp_y[0]);
		// _r51_mul(&temp_x[0], &temp_y[0], &temp_z[0]);
		
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y.inv(); // same as x/y
		REQUIRE(z == ex);
	}

	SECTION("large numbers")
	{
		//values before normalization -> constructor will call carry and normalize
		limb_t temp_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		limb_t temp_y[] = {0x0, 0x6b2c6ce8f04949ac, 0x6bb0f89dbbd5c91d, 0x7057401e6e09dc38, 0xb14cde306209ac72};
		limb_t expected[] = {0x0, 0x3e0dae68d87559d9, 0x6d583f3989782ebb, 0xdc1a212518888180, 0xc05f04242dbfea82};
		
		r51_t x(temp_x);
		r51_t y(temp_y);
		r51_t ex(expected);

		r51_t z = x*y.inv(); // same as x/y
		// printf("result:\t\t%016lx  %016lx  %016lx  %016lx  %016lx\n", z.limbs[0], z.limbs[1], z.limbs[2], z.limbs[3], z.limbs[4]);
		// printf("expected:\t%016lx  %016lx  %016lx  %016lx  %016lx\n", ex.limbs[0], ex.limbs[1], ex.limbs[2], ex.limbs[3], ex.limbs[4]);

		REQUIRE(z == ex);
	}
}

TEST_CASE( "Radix51 - conversion", "[Radix51]" ) {

	SECTION("small numbers"){
		limb_t ulli_x[] = {0x0, 0x0, 0x0, 0x0, 0xC77DE711E};

		r51_t x(ulli_x);

		ulli c[4];
		radixtoulli(&x, c);

		REQUIRE(ulli_x[1] == c[0]);
		REQUIRE(ulli_x[2] == c[1]);
		REQUIRE(ulli_x[3] == c[2]);
		REQUIRE(ulli_x[4] == c[3]);
	}

	SECTION("large numbers"){
		limb_t ulli_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};

		r51_t x(ulli_x);

		ulli c[4];
		radixtoulli(&x, c);

		REQUIRE(ulli_x[1] == c[0]);
		REQUIRE(ulli_x[2] == c[1]);
		REQUIRE(ulli_x[3] == c[2]);
		REQUIRE(ulli_x[4] == c[3]);
	}

	SECTION("large numbers with square"){
		limb_t ulli_x[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
		ulli expected[] = {0x3e0dae68d87559d9, 0x6d583f3989782ebb, 0xdc1a212518888180, 0xc05f04242dbfea82};
		r51_t x(ulli_x);

		x = x.sqr();

		ulli c[4];
		radixtoulli(&x, c);

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
		r51_t x(temp_x);
		r51_t y(temp_y);
		// r51_t ex(expected);

		r51_t z = x*y;

		ulli c[4];
		radixtoulli(&z, c);

		REQUIRE(expected[0] == c[0]);
		REQUIRE(expected[1] == c[1]);
		REQUIRE(expected[2] == c[2]);
		REQUIRE(expected[3] == c[3]);
	}
}

TEST_CASE( "Util - Unhexlify and Hexlify", "[Unhexlify/Hexlify]" )
{
	/*
	SECTION("test1 - 64char string")
	{
		std::string str_in_hex = "9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60";
		ulli *converted_to_ulli = (ulli*) unhexlify(str_in_hex);
		//std::cout<<std::hex<<converted_to_ulli[0]<<" "<<converted_to_ulli[1]<<" "<<converted_to_ulli[2]<<" "<<converted_to_ulli[3]<<std::endl;
	        REQUIRE(converted_to_ulli[0] == 11340540628589959776);
		REQUIRE(converted_to_ulli[1] == 13439949602277829828);
		REQUIRE(converted_to_ulli[2] == 4920681124708116761);
		REQUIRE(converted_to_ulli[3] == 8087246685239934816);

		std::string hexlified = hexlify(converted_to_ulli, 64);

		REQUIRE(hexlified =="9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60");
	}


	SECTION("test2 - 32char string")
	{
		std::string str_in_hex = "4449c5697b326919703bac031cae7f60";
		ulli *converted_to_ulli = (ulli*) unhexlify(str_in_hex);
		std::cout<<std::hex<<converted_to_ulli[0]<<" "<<converted_to_ulli[1]<<" "<<converted_to_ulli[2]<<" "<<converted_to_ulli[3]<<std::endl;
	        REQUIRE(converted_to_ulli[0] == 0);
		REQUIRE(converted_to_ulli[1] == 0);
		REQUIRE(converted_to_ulli[2] == 4920681124708116761);
		REQUIRE(converted_to_ulli[3] == 8087246685239934816);

		std::string hexlified = hexlify(&converted_to_ulli[2], 32);

                REQUIRE(hexlified =="4449c5697b326919703bac031cae7f60");
	}
	*/

}

TEST_CASE("Scalar Multiplication", "[ScalarMult]")
{
	SECTION("identity")
        {
                ulli *a = (ulli*) calloc(5, sizeof(ulli));
                a[0] = 0;
                a[1] = 0;
                a[2] = 0;
                a[3] = 0;
                a[4] = 1;
                Radix51 ar(a);

                ulli A0EXPECTED[5];
                A0EXPECTED[0] = 0;
                A0EXPECTED[1] = 0x216936d3cd6e53fe;
                A0EXPECTED[2] = 0xc0a4e231fdd6dc5c;
                A0EXPECTED[3] = 0x692cc7609525a7b2;
                A0EXPECTED[4] = 0xc9562d608f25d51a;
                Radix51 a0exp((ulli*)&A0EXPECTED);

                ulli A1EXPECTED[5];
                A1EXPECTED[0] = 0;
                A1EXPECTED[1] = 0x6666666666666666;
                A1EXPECTED[2] = 0x6666666666666666;
                A1EXPECTED[3] = 0x6666666666666666;
                A1EXPECTED[4] = 0x6666666666666658;
                Radix51 a1exp((ulli*)&A1EXPECTED);
                r51_t A[2] = {Radix51(ulli_ZERO), Radix51(ulli_ZERO)};
                scalar_mult((r51_t *)CURVE_B2, ar, A);
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

	
	SECTION("test1")
        {
                ulli *a = (ulli*) calloc(5, sizeof(ulli));
                a[0] = 0;
                a[1] = 0x4fe94d9006f020a5;
                a[2] = 0xa3c080d96827fffd;
                a[3] = 0x3c010ac0f12e7a42;
                a[4] = 0xcb33284f86837c30;
                Radix51 ar(a);

		ulli A0EXPECTED[5];
		A0EXPECTED[0] = 0;
		A0EXPECTED[1] = 0x55d0e09a2b9d3429;
		A0EXPECTED[2] = 0x2297e08d60d0f620;
		A0EXPECTED[3] = 0xc513d47253187c24;
		A0EXPECTED[4] = 0xb12786bd777645ce;
		Radix51 a0exp((ulli*)&A0EXPECTED);
		
		ulli A1EXPECTED[5];
		A1EXPECTED[0] = 0;
		A1EXPECTED[1] = 0x1a5107f7681a02af;
		A1EXPECTED[2] = 0x2523a6daf372e10e;
		A1EXPECTED[3] = 0x3a0764c9d3fe4bd5;
		A1EXPECTED[4] = 0xb70ab18201985ad7;
		Radix51 a1exp((ulli*)&A1EXPECTED);
                r51_t A[2] = {Radix51(ulli_ZERO), Radix51(ulli_ZERO)};
                scalar_mult((r51_t *)CURVE_B2, ar, A);
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
	

	
	SECTION("test2")
        {
                ulli *a = (ulli*) calloc(5, sizeof(ulli));
                a[0] = 0;
                a[1] = 0x512e502eb0249a25;
                a[2] = 0x5e1c827f3b6b6c7f;
                a[3] = 0x0a79f4ca8575a915;
                a[4] = 0x28d58258d79ebd68;
                Radix51 ar(a);

		ulli A0EXPECTED[5];
		A0EXPECTED[0] = 0;
		A0EXPECTED[1] = 0x74ad28205b4f384b;
		A0EXPECTED[2] = 0xc0813e6585864e52;
		A0EXPECTED[3] = 0x8085f91fb6a5096f;
		A0EXPECTED[4] = 0x244ae01e57de43ae;
		Radix51 a0exp((ulli*)&A0EXPECTED);
		
		ulli A1EXPECTED[5];
		A1EXPECTED[0] = 0;
		A1EXPECTED[1] = 0x0c66f42af155cdc0;
		A1EXPECTED[2] = 0x8c96c42ecf2c989c;
		A1EXPECTED[3] = 0xbc7e1b4da70ab792;
		A1EXPECTED[4] = 0x5a8943e8c317403d;
		Radix51 a1exp((ulli*)&A1EXPECTED);
                r51_t A[2] = {Radix51(ulli_ZERO), Radix51(ulli_ZERO)};
                scalar_mult((r51_t *)CURVE_B2, ar, A);
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

TEST_CASE("Branchless Scalar Multiplication", "[ScalarMult]")
{
	SECTION("test1")
        {
			ulli *a = (ulli*) calloc(5, sizeof(ulli));
			a[0] = 0;
			a[1] = 0x4fe94d9006f020a5;
			a[2] = 0xa3c080d96827fffd;
			a[3] = 0x3c010ac0f12e7a42;
			a[4] = 0xcb33284f86837c30;
			Radix51 ar(a);

			ulli A0EXPECTED[5];
			A0EXPECTED[0] = 0;
			A0EXPECTED[1] = 0x55d0e09a2b9d3429;
			A0EXPECTED[2] = 0x2297e08d60d0f620;
			A0EXPECTED[3] = 0xc513d47253187c24;
			A0EXPECTED[4] = 0xb12786bd777645ce;
			Radix51 a0exp((ulli*)&A0EXPECTED);
			
			ulli A1EXPECTED[5];
			A1EXPECTED[0] = 0;
			A1EXPECTED[1] = 0x1a5107f7681a02af;
			A1EXPECTED[2] = 0x2523a6daf372e10e;
			A1EXPECTED[3] = 0x3a0764c9d3fe4bd5;
			A1EXPECTED[4] = 0xb70ab18201985ad7;
			Radix51 a1exp((ulli*)&A1EXPECTED);
			r51_t A[2] = {Radix51(ulli_ZERO), Radix51(ulli_ZERO)};
			scalar_mult_branchless((r51_t *)CURVE_B2, ar, A);
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


