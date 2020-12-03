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
  int Ni;                    //pixels in x direction
  int Nj;                    //pixels in y direction
  int Nm;                    //total pixels in the image data
  RectGrid* grid;
  int Nmask;                 //pixels in the mask
  double* defl_x;            //deflected x coordinates
  double* defl_y;            //deflected y coordinates
  int* active;               //active image pixels used in the construction of the adaptive grid
  InterpolationCell** cells = 0;    //for each image pixel, the corresponding source pixels and interpolation weights
  Cross** crosses = 0;       //array of cross structures for creating the D_s(s_p)*D_psi for reconstructing perturbations
  InterpolationCell** dpsi_cells = 0; //for each image pixel, the corresponding dpsi pixels and interpolation weights, if applicable
  mytable B;
  mytable C;
  mytable S;
  std::string noise_flag;

  ImagePlane(const std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);  // used to read images
  ImagePlane(int Ni,int Nj,double xmin,double xmax,double ymin,double ymax);                             // used to create images
  ImagePlane(const ImagePlane& image);
  ~ImagePlane();

  void writeImage(const std::string filename);
  void readB(const std::string filepath,int i,int j,int ci,int cj);
  void readC(const std::string flag,const std::string filepath);
  void readS(const std::string filepath);
  //  void maskData(std::map<int,int> lookup,ImagePlane* masked);
  void printCross(int k);
  void lowerResRebinAdditive(ImagePlane* newImage);
  void lowerResRebinIntegrate(ImagePlane* newImage);

  
private:
  void setCroppedLimitsEven(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset);
  void setCroppedLimitsOdd(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset);
};

#endif /* IMAGE_PLANE_HPP */
