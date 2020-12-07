#ifndef IMAGE_PLANE_HPP
#define IMAGE_PLANE_HPP

#include <string>
#include <vector>
#include <map>

#include "tableDefinition.hpp"

class RectGrid;

class Cross {
public:
  int i0;
  int j0;
  double* coeff_x;
  double* coeff_y;

  Cross(int a,int b){
    this->i0 = a;
    this->j0 = b;
    int n = 8;
    this->coeff_x = (double*) malloc(n*sizeof(double));
    this->coeff_y = (double*) malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
      this->coeff_x[i] = 0.0;
      this->coeff_y[i] = 0.0;
    }
  }
  ~Cross(){
    free(coeff_x);
    free(coeff_y);
  }
};


class ImagePlane {
public:
  int Nx;                    //pixels in x direction
  int Ny;                    //pixels in y direction
  int Nm;                    //total pixels in the image data
  RectGrid* grid;
  int Nmask;                 //pixels in the mask
  mytable B;
  mytable C;
  mytable S;
  std::string noise_flag;

  ImagePlane(const std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);  // used to read images
  ImagePlane(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);                             // used to create images
  ImagePlane(const ImagePlane& image);
  ~ImagePlane();

  void writeImage(const std::string filename);
  void readB(const std::string filepath,int i,int j,int ci,int cj);
  void readC(const std::string flag,const std::string filepath);
  void readS(const std::string filepath);
  
private:
  void setCroppedLimitsEven(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset);
  void setCroppedLimitsOdd(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset);
};

#endif /* IMAGE_PLANE_HPP */
