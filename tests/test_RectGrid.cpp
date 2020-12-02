#include "rectGrid.hpp"
#include "fitsInterface.hpp"

#include <cmath>
#include <iostream>
#include <cstring>

int main(){

  // create a RectGrid object
  int Nx = 85;
  int Ny = 50;
  double b1_min = 4.7035368;
  double b1_max = 8.2338171;
  double b2_min = 4.0344134;
  double b2_max = 7.1372285;
  //  RectGrid mygrid(Nx,Ny,b1_min,b1_max,b2_min,b2_max);
  std::map<std::string,std::string> options{{"dev1_accu","1"},{"dev2_accu","2"},{"interp","bilinear"}};
  RectGrid mygrid(Nx,Ny,b1_min,b1_max,b2_min,b2_max,options);

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
  mygrid.calculate_zx();
  FitsInterface::writeFits(Nx,Ny,mygrid.zx,"vkl_zx.fits");

  // calculate the first derivative along y
  mygrid.calculate_zy();
  FitsInterface::writeFits(Nx,Ny,mygrid.zx,"vkl_zy.fits");

  // calculate the mixed derivative along x and y
  mygrid.calculate_zxy();
  FitsInterface::writeFits(Nx,Ny,mygrid.zxy,"vkl_zxy.fits");

  // calculate the second derivative along x
  mygrid.calculate_zxx();
  FitsInterface::writeFits(Nx,Ny,mygrid.zxx,"vkl_zxx.fits");

  // calculate the second derivative along y
  mygrid.calculate_zyy();
  FitsInterface::writeFits(Nx,Ny,mygrid.zyy,"vkl_zyy.fits");









  

  // create a RectGrid object
  Nx = 5;
  Ny = 5;
  b1_min = -4.;
  b1_max = 0.2;
  b2_min = -2.0;
  b2_max = 2.1;
  options = {{"dev1_accu","1"},{"dev2_accu","2"},{"interp","bicubic"}};
  RectGrid mygrid2(Nx,Ny,b1_min,b1_max,b2_min,b2_max,options);
  /*
  for(int j=0;j<mygrid2.Nx;j++){
    std::cout << mygrid2.center_x[j] << " ";
  }
  std::cout << std::endl;
  std::cout << "stepx: " << mygrid2.step_x << std::endl;
  for(int i=0;i<mygrid2.Ny;i++){
    std::cout << mygrid2.center_y[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "stepy: " << mygrid2.step_y << std::endl;
  */

  // calculate some function z at its pixels
  //  double dum[25] = {0.8,1.0,0.3,0.7,0.9,0.7,0.4,0.6,0.1,0.2,0.4,0.7,0.2,0.9,0.0,0.1,0.2,0.3,0.4,0.5,0.4,0.7,0.0,0.3,0.1};
  double dum[25] = {0.4,0.7,0.,0.3,0.1,0.1,0.2,0.3,0.4,0.5,0.4,0.7,0.2,0.9,0.,0.7,0.4,0.6,0.1,0.2,0.8,1.,0.3,0.7,0.9};
  std::memcpy(mygrid2.z,dum,Nx*Ny*sizeof(mygrid2.z));
  

  
  //RectGrid mygrid_out(10*Nx,10*Ny,mygrid2.center_x[0],mygrid2.center_x[mygrid2.Nx-1],mygrid2.center_y[0],mygrid2.center_y[mygrid2.Ny-1],options);
  /*
  RectGrid mygrid_out(30*Nx,30*Ny,b1_min+1,b1_max-1,b2_min+1,b2_max-1,options);

  for(int j=0;j<mygrid_out.Nx;j++){
    std::cout << mygrid_out.center_x[j] << " ";
  }
  std::cout << std::endl;
  std::cout << "stepx: " << mygrid_out.step_x << std::endl;
  for(int i=0;i<mygrid_out.Ny;i++){
    std::cout << mygrid_out.center_y[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "stepy: " << mygrid_out.step_y << std::endl;

  
  for(int i=0;i<mygrid_out.Ny;i++){
    for(int j=0;j<mygrid_out.Nx;j++){
      double x = mygrid_out.center_x[j];
      double y = mygrid_out.center_y[i];
      if( mygrid2.point_between_pixel_centers(x,y,0) ){
	mygrid_out.z[i*mygrid_out.Nx+j] = mygrid2.interp2d(x,y);
      } else {
	mygrid_out.z[i*mygrid_out.Nx+j] = 0;
      }
    }
  }
  */
  RectGrid* mygrid_out = mygrid2.embeddedNewGrid(30*Nx,30*Ny);

  FitsInterface::writeFits(mygrid_out->Nx,mygrid_out->Ny,mygrid_out->z,"vkl_interp.fits");
  
  
}

