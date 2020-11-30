#ifndef FITS_INTERFACE_HPP
#define FITS_INTERFACE_HPP

#include <string>
#include <vector>

class FitsInterface{
public:
  static void readFits(int Ni,int Nj,double* z,const std::string filepath);
  static void writeFits(int Ni,int Nj,double* z,const std::string filepath);
  static void writeFits(int Ni,int Nj,double* z,std::vector<std::string> key,std::vector<std::string> value,std::vector<std::string> description,const std::string filepath);
};

#endif /* FITS_INTERFACE_HPP */
