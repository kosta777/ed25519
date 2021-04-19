#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "const.h"
#include "tsc_x86.h"
#include "util.h"
#include "ed25519.h"
#include <locale>
#include <fstream>
#include <cstring>
#include <immintrin.h>
#include "rnd_numbers.hpp"
// #include "rnd_numbers_2M.hpp"


// this is done by makefile
#if STDC
  #include "StdC-optimize/scalarmult/scalarmult.h"
  #if R25
    #include "StdC-optimize/radix25/radix25.h"
  #else
    #include "StdC-optimize/radix/radix51.h"
  #endif
#elif BASE
  #include "straight/scalarmult.h"
  #include "straight/radix51.h"
#elif SSE
#include "SSE/scalarmult/scalarmult.h"
  #if R25
    #include "SSE/radix25/radix25.h"
  #else
    #include "SSE/radix/radix51.h"
  #endif
#endif

// The number of keys to sign is selected there
#define KEY_INPUT "keygen_input/input_10.in"
#define N 10
#define PK_LABEL "public key (10 inputs)"
char *sk_unhexed;
char *pk_computed;

limb_t limb_sink[] = {0x0, 0x0, 0x0, 0x0, 0x0};
limb_t temp_P0[] = {0x0, 0x166E77, 0x0893B30FC1AC0F43, 0x65206AE5BF9E7CEE, 0x4838BF8C77DE711E};
limb_t temp_P1[] = {0x0, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};
limb_t temp_e[] = {0x0, 0x430D115, 0x00171CB717C93B55, 0x30C23D546FC409AA, 0x5862711F01A8B173};


radix_t sink(limb_sink);
radix_t P0(temp_P0);
radix_t P1(temp_P1);
radix_t e(temp_e);
radix_t P[] = {P0, P1};
radix_t P_sink[] = {ONE, ONE};

std::vector<radix_t> rnd_nums;
std::vector<__m256i*> rnd_nums_converted;

__m256i converted_y[10];
__m256i converted_z[10];
unsigned long NUM_INPUTS = 4; // default value, does not need to be a multiple of 4, can be passed as command line argument
unsigned long NUM_VECTORS;
#define NUM_BENCHMARKS 9
const char *func_names[NUM_BENCHMARKS] = {"multiplication", "squaring", "inverse", "scalar_mult", "v_mul", "v_sqr", "v_inv", "v_scalar_mult", PK_LABEL};

const bool vectorized[NUM_BENCHMARKS] = {false, false, false, false, true, true, true, true, false};

#if SSE | STDC
  #if R25
  unsigned long ops_count[NUM_BENCHMARKS] = {241, 156, 42275, 44448164, 241, 156, 42275, 44448164, 0};
  #else //                                   mul  sqr   inv   scalmul  vmul vsqr vinv  vscalmul  pk
  unsigned long ops_count[NUM_BENCHMARKS] = {116, 89, 23882, 24682219, 116, 89, 23882, 24682219, 0};
  #endif
#else  // BASE
  #if R25
  unsigned long ops_count[NUM_BENCHMARKS] = {312, 157, 43310, 54747808, 312, 157, 43310, 54747808, 0};
  #else
  unsigned long ops_count[NUM_BENCHMARKS] = {917, 569, 154613, 162309635, 917, 569, 154613, 162309635, 0};
  #endif
#endif

#if SSE && R25 
  limb_t temp1[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
  limb_t temp2[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
  limb_t temp3[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
  limb_t temp4[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};
  limb_t temp5[] = {0x0, 0x1d2f4909d7fa7539, 0x205ec73460a744ca, 0x2959fdc8e114e1ea, 0x64d65005d3e0a09b};
  limb_t temp6[] = {0x0, 0x7df9b56e994884d3, 0x6f54b68e77ee106f, 0x402a242e2799b3eb, 0x5ce89754aa0cae48};
  limb_t temp7[] = {0x0, 0x41ac5253aaed04c4, 0x32a675e52a3fc8e2, 0x517d3c1e678c8312, 0x2a6598d95e75122d};
  limb_t temp8[] = {0x0, 0x214c921a9c0f39cc, 0x22b5e894467b7588, 0x12f7766ce9e6cce4, 0x2262f8aca9c26644};

  limb_t temp0[] = {0, 0, 0, 0, 0};

  r25_t _x1(temp1);
  r25_t _x2(temp2);
  r25_t _x3(temp3);
  r25_t _x4(temp4);
  r25_t _y1(temp1);
  r25_t _y2(temp2);
  r25_t _y3(temp3);
  r25_t _y4(temp4);

  r25_t vp_P[4][2];

#endif

void run_func(int index)
{
  switch (index)
  {
  case 0: // MUL
    sink = P0 * P1;
    break;
  case 1: // SQR
    sink = P0.sqr();
    break;
  case 2: // INV
    sink = P0.inv();
    break;
  case 3: { // SCALAR_MULT
#if BASE
    scalar_mult(P, e, P_sink);
#else
    scalar_mult_branchless(P, e, P_sink);
#endif
  }
    break;
  case 4: { // V_MUL
  #if SSE
    #if R25
      for (size_t i = 0; i < NUM_VECTORS; i++)
      {
        __m256i* converted_x = rnd_nums_converted[i];
        // __m256i* converted_y = rnd_nums_converted[2*i+1];
        _v_r25_mul(converted_x, converted_y, converted_z);
      }
    #else
      for (size_t i = 0; i < NUM_INPUTS; i++)
      {
        limb_t* x = rnd_nums[i].limbs;
        limb_t* y = P0.limbs;
        _r51_mul(x, y, sink.limbs);
      }
    #endif
  #else // StdC | BASE
    for (size_t i = 0; i < NUM_INPUTS; i++)
      {
        radix_t x = rnd_nums[i];
        sink = x*P0;
      }
  #endif
  }
    break;
  case 5: { // V_SQR
  #if SSE
    #if R25
      for (size_t i = 0; i < NUM_VECTORS; i++)
      {
        __m256i* converted_x = rnd_nums_converted[i];
        _v_r25_sqr(converted_x, converted_z);
      }
    #else
      for (size_t i = 0; i < NUM_INPUTS; i++)
      {
        limb_t* x = rnd_nums[i].limbs;
        // limb_t* y = P0.limbs;
        _r51_sqr(x, sink.limbs);
      }
    #endif
  #else // STDC | BASE
    for (size_t i = 0; i < NUM_INPUTS; i++)
      {
        radix_t x = rnd_nums[i];
        sink = x.sqr();
      }
  #endif
  }
    break;
  case 6: // V_INV
  {
  #if SSE
    #if R25
      for (size_t i = 0; i < NUM_VECTORS; i++)
      {
        __m256i* converted_x = rnd_nums_converted[i];
        _v_r25_inv(converted_x, converted_z);
      }
    #else
      for (size_t i = 0; i < NUM_INPUTS; i++)
      {
        limb_t* x = rnd_nums[i].limbs;
        _r51_inv(x, sink.limbs);
      }
    #endif
  #else // STDC | BASE
    for (size_t i = 0; i < NUM_INPUTS; i++)
      {
        radix_t x = rnd_nums[i];
        sink = x.inv();
      }
  #endif
  }
    break;
  case 7: // V_SCALAR_MUL
  {
  #if SSE & R25
    for (size_t i = 0; i < NUM_VECTORS; i++)
    {
      radix_t e0 = rnd_nums[i];
      radix_t e1 = rnd_nums[i+1];
      radix_t e2 = rnd_nums[i+2];
      radix_t e3 = rnd_nums[i+3];
      radix_t E[4]; E[0] = e0; E[1] = e1; E[2] = e2; E[3] = e3;
      v_scalar_mult(vp_P, E, P_sink, P_sink, P_sink, P_sink);
    }
  #else
    for (size_t i = 0; i < NUM_INPUTS; i++)
    {
      radix_t e0 = rnd_nums[i];
//      scalar_mult(P, e0, P_sink);
#if BASE
      scalar_mult(P, e0, P_sink);
#else
      scalar_mult_branchless(P, e0, P_sink);
#endif
    }
  #endif
  }
    break;
  case 8: { // PUBLIC KEYS
    publickeys(N, (char (*)[32])sk_unhexed, (char (*)[64])pk_computed);
  }
    break;

  default:
    break;
  }
}



double median(std::vector<double> &v)
{
  size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + n, v.end());
  return v[n];
}

void benchmark(int index)
{
  long num_runs = 1;
  myInt64 start, end;
  double multiplier = 1;  
  double cycles = 0;
  long required = 1e8;

  do
  {
    num_runs = num_runs * multiplier;
    start = start_tsc();
    for (size_t i = 0; i < num_runs; i++)
    {
      run_func(index);
    }
    end = stop_tsc(start);

    cycles = (double)end;
    multiplier = (required) / (cycles);

  } while (multiplier > 2);

  // std::cout << "num runs: " << num_runs << std::endl;

  std::vector<double> cyclesList;
  for (int j = 0; j < 100; j++)
  {
    start = start_tsc();
    for (int i = 0; i < num_runs; i++)
    {
      run_func(index);
    }
    end = stop_tsc(start);
    cyclesList.push_back((double)end / num_runs);
    // std::cout << ".";
  }
  double median_cyc = median(cyclesList);
  double ops_per_cyc = (double)ops_count[index] / median_cyc;
  // std::cout << std::endl << std::fixed  << std::setprecision(2)
  //           << func_names[index] << " took " << median_cyc << " cycles; "
  //           << ops_per_cyc << " ops/cycle" << std::endl;
    std::cout << std::fixed  << std::setprecision(2)
            << func_names[index] << ";" << median_cyc << ";"
            << ops_per_cyc << std::endl;
}

int main(int argc, char **argv) {
  // read NUM_INPUTS from command line
  if (argc==2) {
    NUM_INPUTS = atoi(argv[1]);
  }

  // it's still a mix between variable size and single inputs
  for (int i = 0; i < NUM_BENCHMARKS; i++)
  {
    if (vectorized[i]) ops_count[i] *= NUM_INPUTS;
  }

  // how many vectors are needed
  NUM_VECTORS = (3+NUM_INPUTS)/4;

  std::cout.setf(std::ios::unitbuf);
  // std::cout.imbue(std::locale("")); // uncomment to display thousand separator
  
  // Parsing the keys
  std::string line;
  std::ifstream infile(KEY_INPUT);
  if (!infile.is_open()) {
    std::cerr << "Error opening keys file" << std::endl;
    return 1;
  }
  std::string delimiter = ":";
  std::vector<std::string> sk;
  std::vector<std::string> pk;
  std::string token;
  std::string tokens[2];
  size_t pos;
  while(std::getline(infile, line)) {
    int i = 0;
    while ((pos = line.find(delimiter)) != std::string::npos) {
      token = line.substr(0, pos);
      line.erase(0, pos + delimiter.length());
      tokens[i++] = token;
    }
    if(i != 2) {
      std::cout<<"Error in one of the lines"<<std::endl;
      break;
    }
    sk.push_back(tokens[0]);
    pk.push_back(tokens[1]);
  }
  
  int n = sk.size();
  sk_unhexed = (char *)calloc(n, 32);
  pk_computed = (char *)calloc(n, 64);
  for (int i = 0 ; i < n ; i++) {
    memcpy(sk_unhexed+(i*32), unhexlify(sk[i].c_str(),0), 32);
  }

  // prepare numbers
  for (size_t i = 0; i < NRND_NUM; i++)
  {
    ulli bar[5];
    bar[0] = rnd_numbers[i*5+0];
    bar[1] = rnd_numbers[i*5+1];
    bar[2] = rnd_numbers[i*5+2];
    bar[3] = rnd_numbers[i*5+3];
    bar[4] = rnd_numbers[i*5+4];
    radix_t foo(bar);
    rnd_nums.push_back(foo);
  }

  #if SSE & R25
  // prepare x values
  for (size_t j = 0; j < NRND_NUM; j+=4)
  {
    __m256i converted_x [10];
    // prepare vectorized limbs
    for(int i=0;i<10;i++)
    {
      converted_x[i] = _mm256_set_epi64x(rnd_nums[j+3].limbs[i], rnd_nums[j+2].limbs[i], rnd_nums[j+1].limbs[i], rnd_nums[j].limbs[i]);
    }
    rnd_nums_converted.push_back(converted_x);
  }
  // prepare a single y
  // TODO: maybe change to have different y each run
  for (int i = 0; i < 10; i++)
  {
    converted_y[i] = _mm256_set_epi64x(_y4.limbs[i], _y3.limbs[i], _y2.limbs[i], _y1.limbs[i]);
  }

  // prepare point for v_scalar_mult
  vp_P[0][0] = _x1; vp_P[0][1] = _y1;
  vp_P[1][0] = _x2; vp_P[1][1] = _y2;
  vp_P[2][0] = _x3; vp_P[2][1] = _y3;
  vp_P[3][0] = _x4; vp_P[3][1] = _y4;

  std::cerr << "Numbers prepared:\trnd_nums len: " << rnd_nums.size() << " rnd_nums_converted len: " << rnd_nums_converted.size() << std::endl;
  #endif
  std::cout<< "NUM_INPUTS: " << NUM_INPUTS << std::endl;
  std::cout<< "NUM_VECTORS: " << NUM_VECTORS << std::endl;
  
  std::cout << "name;cycles;ops/cycle" << std::endl;

  // Benchmarking
  for (int i = 0; i < NUM_BENCHMARKS; i++)
  {
    benchmark(i);
  }
  
  
  free(sk_unhexed);
  free(pk_computed);
}
