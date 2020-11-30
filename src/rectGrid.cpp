#include "rectGrid.hpp"

#include <CCfits/CCfits>


RectGrid::RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  this->common_constructor(Nx,Ny,xmin,xmax,ymin,ymax);
}

RectGrid::RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,const std::string filepath){
  this->common_constructor(Nx,Ny,xmin,xmax,ymin,ymax);

  // Now read the input fits file into z
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  image.readAllKeys();
  // Here on can check if sizes agree
  std::valarray<float> contents(image.axis(0)*image.axis(1));
  this->readFits(filepath,contents);

  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->z[i*Nx+j] = contents[i*Nx+j];
    }
  }
}

void RectGrid::common_constructor(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  this->Nz     = this->Nx*this->Ny;
  this->width  = xmax - xmin;
  this->height = ymax - ymin;
  this->step_x = this->width/this->Nx;
  this->step_y = this->height/this->Ny;
  
  for(int i=0;i<this->Ny;i++){
    this->center_y[i] = ymin + (this->step_y/2.0) + i*this->step_y;
    this->bound_y[i]  = ymin + i*this->step_y;
  }  
  this->bound_y[this->Ny] = ymax;
  
  for(int j=0;j<this->Nx;j++){
    this->center_x[j] = xmin + (this->step_x/2.0) + j*this->step_x;
    this->bound_x[j]  = xmin + j*this->step_x;
  }
  this->bound_x[this->Nx] = xmax;
}

void RectGrid::readFits(const std::string filepath,std::valarray<float>& contents){
  std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filepath,CCfits::Read,true));
  CCfits::PHDU& image = pInfile->pHDU();
  
  std::valarray<float> tmp;
  image.readAllKeys();
  image.read(tmp);
  
  int Nj = image.axis(0);
  int Ni = image.axis(1);

  //convert FITS standard (bottom to top) to the one used in this code (top to bottom)
  for(int j=0;j<Nj;j++){
    for(int i=0;i<Ni;i++){
      contents[j*Ni+i] = tmp[(Nj-j-1)*Ni+i];
    }
  }
}

RectGrid::~RectGrid(){
  free(center_x);
  free(center_y);
  free(z);
  free(bound_x);
  free(bound_y);
}

bool RectGrid::point_in_grid(double x,double y){
  if( x < this->bound_x[0] || this->bound_x[this->Nx] < x || y < this->bound_y[0] || this->bound_y[this->Ny] < y ){
    return false;
  } else {
    return true;
  }
}

bool RectGrid::point_between_pixel_centers(double x,double y,int boundary_size){
  if( x < this->center_x[0] || this->center_x[this->Nx-1-boundary_size] < x || y < this->center_y[0] || this->center_y[this->Ny-1-boundary_size] < y ){
    return false;
  } else {
    return true;
  }
}

void RectGrid::match_point_to_pixel(double x,double y,int& i0,int& j0){
  i0 = (int) floor( (y-this->bound_y[0])/this->step_y );
  j0 = (int) floor( (x-this->bound_x[0])/this->step_x );
}


