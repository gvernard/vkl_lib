#include "rectGrid.hpp"
#include "fitsInterface.hpp"

#include <cmath>
#include <iostream>

int main(){

  // create a RectGrid object
  int Nx = 85;
  int Ny = 50;
  double b1_min = 4.7035368;
  double b1_max = 8.2338171;
  double b2_min = 4.0344134;
  double b2_max = 7.1372285;
  RectGrid mygrid(Nx,Ny,b1_min,b1_max,b2_min,b2_max);


  // calculate some function z at its pixels
  for(int i=0;i<mygrid.Ny;i++){
    double y = mygrid.center_y[i];  
    for(int j=0;j<mygrid.Nx;j++){
      double x = mygrid.center_x[j];
      mygrid.z[i*mygrid.Nx+j] = 2*cos(x*y/3.0) + cos(5.0*x);
    }
  }
  FitsInterface::writeFits(Nx,Ny,mygrid.z,"vkl_z.fits");
  

  // calculate the first derivative along x
  mygrid.calculate_zx("1");
  FitsInterface::writeFits(Nx,Ny,mygrid.zx,"vkl_zx.fits");

  // calculate the first derivative along y
  mygrid.calculate_zy("1");
  FitsInterface::writeFits(Nx,Ny,mygrid.zx,"vkl_zy.fits");

  // calculate the mixed derivative along x and y
  mygrid.calculate_zxy("1");
  FitsInterface::writeFits(Nx,Ny,mygrid.zxy,"vkl_zxy.fits");

  // calculate the second derivative along x
  mygrid.calculate_zxx("2");
  FitsInterface::writeFits(Nx,Ny,mygrid.zxx,"vkl_zxx.fits");

  // calculate the second derivative along y
  mygrid.calculate_zyy("2");
  FitsInterface::writeFits(Nx,Ny,mygrid.zyy,"vkl_zyy.fits");

    
}
