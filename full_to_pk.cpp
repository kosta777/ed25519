// only goal of this file is to convert a full (skpk:pk:m:sig: format) input file into a key-only input file (sk:pk: format)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
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
  
  // Reading all the keys in the file
  std::string delimiter = ":";
  std::vector<std::string> sk;
  std::vector<std::string> pk;
  std::string token;
  std::string tokens[4];
  size_t pos;
  
  
  std::ofstream out("output.txt");
  
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
    
    out << tokens[0].substr(0, 64) << ":" << tokens[1] << ":" << std::endl;
  }
  out.close();

  return 0;
}
