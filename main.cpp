#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include "ed25519.h"
// #include "StdC-optimize/radix/radix51.h"
#include "util.h"
#include "tsc_x86.h"

#include <vector>

typedef u_int64_t ulli;

int main(int argc, char const *argv[])
{
  std::string line;

  if(argc !=2)
  {
    std::cout<<"Usage: ./proj path_to_file"<<std::endl;
    return 0;
  }
  
  std::ifstream infile(argv[1]);
  if (!infile.is_open())
  {
    std::cout<<"Error opening file"<<std::endl;
      return 1;
  }
  
  std::cerr << "Parsing input... ";
  // Reading all the keys in the file
  std::string delimiter = ":";
  std::vector<std::string> sk;
  std::vector<std::string> pk;
  std::string token;
  std::string tokens[2];
  size_t pos;
  while(std::getline(infile, line))
  {
  //  std::cout<<line<<std::endl;
    int i=0;
    while ((pos = line.find(delimiter)) != std::string::npos)
    {
      token = line.substr(0, pos);
      //std::cout << token << std::endl;
      line.erase(0, pos + delimiter.length());
      tokens[i++] = token;
    }
    
    if(i != 2)
    {
      std::cout<<"Error in one of the lines"<<std::endl;
      break;
    }
    
    sk.push_back(tokens[0]);
//    std::cerr << sk[sk.size() - 1] << std::endl;
    pk.push_back(tokens[1]);
  }
  std::cerr << "done" << std::endl;
    
  // Unhexlifying all the private keys
  std::cerr << "Converting from text to hex... ";
  int n = sk.size();
  char sk_unhexed[n][32];
  char pk_computed[n][64];
  for (int i = 0 ; i < n ; i++) {
    memcpy(sk_unhexed[i], unhexlify(sk[i].c_str(),0), 32);
  }
  std::cerr << "done" << std::endl;
  
//  for (int i = 0 ; i < n ; i++) tmp_print_bytes(sk_unhexed[i], 32); //std::cerr << hexlify(sk_unhexed[i], 32) << std::endl;
  
  // Computing all the public keys
  std::cerr << "Computing the public keys... ";
  publickeys(n, sk_unhexed, pk_computed);
  std::cerr << "done" << std::endl;
  
  std::cerr << "Verifying... ";
  int wrong = 0;
//  std::cerr << std::endl;
  for (int i = 0 ; i < n ; i++) {
//    std::cerr << i << ": " << (strncmp(pk_computed[i], pk[i].c_str(), 64) ? "WRONG" : "ok") << ". " << std::endl;
    if (strncmp(pk_computed[i], pk[i].c_str(), 64)) {
//      std::cerr << "Expected: " << pk[i].c_str() << std::endl;
//      std::cerr << "Obtained: "; for (int j = 0 ; j < 64 ; j++) std::cerr << pk_computed[i][j]; std::cerr << std::endl;
      wrong++;
    }
  }
  std::cerr << wrong << " public keys were wrong." << std::endl;
  
  return 0;
}
    
/*
  std::string delimiter = ":";
  std::vector<std::string> sk;
  std::vector<std::string> pk;
  std::string m;
  std::string sig;
  std::string token;
  std::string tokens[4];
  size_t pos;
  int ok;
  int count = 0;
  int count_wrong_pk = 0;
  while(std::getline(infile, line))
  {
  //  std::cout<<line<<std::endl;
    int i=0;
    while ((pos = line.find(delimiter)) != std::string::npos)
    {
      token = line.substr(0, pos);
      //std::cout << token << std::endl;
      line.erase(0, pos + delimiter.length());
      tokens[i++] = token;
    }
    
    if(i != 4)
    {
      std::cout<<"Error in one of the lines"<<std::endl;
      break;
    }
    sk = tokens[0].substr(0, 64);
    pk = tokens[1];
    m = tokens[2];
    sig = tokens[3];
    // Unhexlify returns ulli array with 4 elements
    char *sk_unhexed = unhexlify(sk.c_str(),0);
    char *pk_unhexed = (unhexlify(pk.c_str(),16)); // Previously 16 - why ?
    
/// PARAMETERS
//    std::cerr << "PARAMETERS:" << std::endl;
//    std::cerr << "sk = " << sk << " (" << strlen(sk_unhexed) << ")" << std::endl;
//    std::cerr << "pk = " << pk << " (" << strlen(pk_unhexed) << ")" << std::endl;
//
//    char *hash = (char *)malloc(SHA512_DIGEST_LENGTH);
//    H((const unsigned char *)sk_unhexed, 32, hash);
//    char *hash_hex = hexlify(hash, 64);
//    std::cerr << "sk = "; tmp_print_bytes(sk_unhexed, 32);
//    std::cerr << "SHA512(sk_unhexed)  = " << hash_hex << std::endl;
//    free(hash);
//    free(hash_hex);
//
//    std::cerr << "-----------------------" << std::endl;
    
    
    std::cerr << "Input " << count++ << ": ";
    
/// PUBLIC KEY
    std::cerr << "Testing publickey()... ";
    char *pk_computed = (char *)calloc(65, sizeof(char));
    pk_computed[64] = '\0';
    publickey(sk_unhexed, pk_computed); // Benchmark inside publickey()
    std::cerr << (strncmp(pk_computed, pk.c_str(), 64) ? "WRONG" : "ok") << ". ";
    if (strncmp(pk_computed, pk.c_str(), 64)) count_wrong_pk++;
    free(pk_computed);

 
/// SIGNATURE
//    std::cerr << "Testing signature()... ";
//    char *sig_computed = (char *)calloc(129, sizeof(char));
//    sig_computed[128] = '\0';
//    signature(&m[0], sk_unhexed, pk_unhexed, sig_computed);
//    std::cerr << (strncmp(sig_computed, sig.c_str(), 128) ? "WRONG" : "ok") << ". ";
//    if (strncmp(sig_computed, sig.c_str(), 128)) count_wrong_sig++;
//    free(sig_computed);

    
/// VERIFICATION
//    std::cerr << "Verifying signature... ";
//    int res = verify(&sig[0], &m[0], &pk[0]);
//    std::cerr << (res ? "FAIL" : "ok") << ".";
//    if (res) count_wrong_ver++;
  

/// FREEING MEMORY
    free(sk_unhexed);
    free(pk_unhexed);
    std::cerr << std::endl;
//    return 0; // Uncomment to debug individual functions
  }
  
  
  std::cerr << "-----------------------" << std::endl;
  std::cerr << "Finished processing " << count << " inputs." << std::endl;
  std::cerr << count_wrong_pk << " public key(s) were wrong (success rate of " << (100 - 100*(count_wrong_pk/(float)count)) << " %)." << std::endl;
//  std::cerr << count_wrong_sig << " signature(s) were wrong (success rate of " << (100 - 100*(count_wrong_sig/(float)count)) << " %)." << std::endl;
//  std::cerr << count_wrong_ver << " verification(s) were wrong (success rate of " << (100 - 100*(count_wrong_ver/(float)count)) << " %)." << std::endl;
  
  return 0;
}
*/
