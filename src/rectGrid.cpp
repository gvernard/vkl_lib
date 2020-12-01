#include "rectGrid.hpp"
#include "fitsInterface.hpp"

#include <cmath>
#include <algorithm>

RectGrid::RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  this->common_constructor(Nx,Ny,xmin,xmax,ymin,ymax);
}

RectGrid::RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,const std::string filepath){
  this->common_constructor(Nx,Ny,xmin,xmax,ymin,ymax);
  this->readFits(this->z,filepath);
}

void RectGrid::common_constructor(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  this->Nx = Nx;
  this->Ny = Ny;
  this->width  = xmax - xmin;
  this->height = ymax - ymin;
  this->step_x = this->width/this->Nx;
  this->step_y = this->height/this->Ny;

  this->Nz = this->Nx*this->Ny;
  this->z  = (double*) calloc(this->Nx*this->Ny,sizeof(double));

  this->center_y = (double*) calloc(this->Ny,sizeof(double));
  this->bound_y  = (double*) calloc(this->Ny+1,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    this->center_y[i] = ymin + (this->step_y/2.0) + i*this->step_y;
    this->bound_y[i]  = ymin + i*this->step_y;
  }  
  this->bound_y[this->Ny] = ymax;
  
  this->center_x = (double*) calloc(this->Nx,sizeof(double));
  this->bound_x  = (double*) calloc(this->Nx+1,sizeof(double));
  for(int j=0;j<this->Nx;j++){
    this->center_x[j] = xmin + (this->step_x/2.0) + j*this->step_x;
    this->bound_x[j]  = xmin + j*this->step_x;
  }
  this->bound_x[this->Nx] = xmax;
}

void RectGrid::readFits(double* field,const std::string filepath){
  FitsInterface::readFits(this->Ny,this->Nx,field,filepath);
}

void multiplyVectorByScalar(std::vector<double> &v,double k){
  std::transform(v.begin(),v.end(),v.begin(), [k](double &c){ return c*k; });
}


RectGrid::~RectGrid(){
  free(center_x);
  free(center_y);
  free(z);
  free(zx);
  free(zy);
  free(zxy);
  free(zxx);
  free(zyy);
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

bool RectGrid::match_point_to_closest_4(double x,double y,int* i,int* j){
  if( !this->point_between_pixel_centers(x,y,0) ){
    return false;
  }
  int i0,j0;
  match_point_to_pixel(x,y,i0,j0);
  double dy = y - this->center_y[i0];
  double dx = x - this->center_x[j0];
  int f1,f2;
  if( dx<0 ){
    f1 = -1;
    f2 = 0;
  } else {
    f1 = 0;
    f2 = 1;
  }
  j[0] = j0 + f1;
  j[1] = j0 + f2;
  j[2] = j0 + f1;
  j[3] = j0 + f2;
  if( dy<0 ){
    f1 = 0;
    f2 = -1;
  } else {
    f1 = 1;
    f2 = 0;
  }
  i[0] = i0 + f1;
  i[1] = i0 + f1;
  i[2] = i0 + f2;
  i[3] = i0 + f2;  
  return true;
}

bool RectGrid::match_point_to_closest_16(double x,double y,int* i,int* j){
  if( !this->point_between_pixel_centers(x,y,1) ){
    return false;
  }
  int i0,j0;
  match_point_to_pixel(x,y,i0,j0);
  double dy = y - this->center_y[i0];
  double dx = x - this->center_x[j0];
  int f1,f2,f3,f4;
  if( dx<0 ){
    f1 = -2;
    f2 = -1;
    f3 = 0;
    f4 = 1;
  } else {
    f1 = -1;
    f2 = 0;
    f3 = 1;
    f4 = 2;
  }
  j[0]  = j0 + f1;
  j[1]  = j0 + f2;
  j[2]  = j0 + f3;
  j[3]  = j0 + f4;
  j[4]  = j0 + f1;
  j[5]  = j0 + f2;
  j[6]  = j0 + f3;
  j[7]  = j0 + f4;
  j[8]  = j0 + f1;
  j[9]  = j0 + f2;
  j[10] = j0 + f3;
  j[11] = j0 + f4;
  j[12] = j0 + f1;
  j[13] = j0 + f2;
  j[14] = j0 + f3;
  j[15] = j0 + f4;
  if( dy<0 ){
    f1 = -1;
    f2 = 0;
    f3 = 1;
    f4 = 2;
  } else {
    f1 = 2;
    f2 = 1;
    f3 = 0;
    f4 = -1;
  }
  i[0]  = i0 + f1;
  i[1]  = i0 + f1;
  i[2]  = i0 + f1;
  i[3]  = i0 + f1;  
  i[4]  = i0 + f2;
  i[5]  = i0 + f2;
  i[6]  = i0 + f2;
  i[7]  = i0 + f2;  
  i[8]  = i0 + f3;
  i[9]  = i0 + f3;
  i[10] = i0 + f3;
  i[11] = i0 + f3;  
  i[12] = i0 + f4;
  i[13] = i0 + f4;
  i[14] = i0 + f4;
  i[15] = i0 + f4;  
  return true;
}


void RectGrid::calculate_zx(std::string accuracy){
  this->zx = (double*) calloc(this->Nz,sizeof(double));
  this->calculate_derivative_1(this->Nx,this->Ny,this->center_x,this->z,this->zx,accuracy);
}

void RectGrid::calculate_zy(std::string accuracy){
  this->zy = (double*) calloc(this->Nz,sizeof(double));
  double* tmp_1 = (double*) calloc(this->Nz,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zy[j*this->Ny+i] = this->z[i*this->Nx+j];
    }
  }
  this->calculate_derivative_1(this->Ny,this->Nx,this->center_y,this->zy,tmp_1,accuracy);
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zy[i*this->Nx+j] = tmp_1[j*this->Ny+i];
    }
  }
  free(tmp_1);
}

void RectGrid::calculate_zxy(std::string accuracy){
  if( this->zy == NULL ){
    this->calculate_zy(accuracy);
  }
  this->zxy = (double*) calloc(this->Nz,sizeof(double));
  this->calculate_derivative_1(this->Nx,this->Ny,this->center_x,this->zy,this->zxy,accuracy);  
}

void RectGrid::calculate_zxx(std::string accuracy){
  this->zxx = (double*) calloc(this->Nz,sizeof(double));
  this->calculate_derivative_2(this->Nx,this->Ny,this->center_x,this->z,this->zxx,accuracy);
}

void RectGrid::calculate_zyy(std::string accuracy){
  this->zyy = (double*) calloc(this->Nz,sizeof(double));
  double* tmp_1 = (double*) calloc(this->Nz,sizeof(double));
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zyy[j*this->Ny+i] = this->z[i*this->Nx+j];
    }
  }
  this->calculate_derivative_2(this->Ny,this->Nx,this->center_y,this->zyy,tmp_1,accuracy);
  for(int i=0;i<this->Ny;i++){
    for(int j=0;j<this->Nx;j++){
      this->zyy[i*this->Nx+j] = tmp_1[j*this->Ny+i];
    }
  }
  free(tmp_1);
}


void RectGrid::calculate_derivative_1(int Nh,int Nv,double* h,double* zz,double* zout,std::string accuracy){
  // calculate the FIRST derivative along the horizontal axis
  int h0,v0;
  std::vector<int> rel_index_v;
  std::vector<int> rel_index_h;
  std::vector<double> coeff;
  double dh = h[1]-h[0];
  
  // 1st column
  if( accuracy == "1" ){
    rel_index_v = {0,0};
    rel_index_h = this->derivative_1_forward_1_index;
    coeff       = this->derivative_1_forward_1_coeff;
  } else {
    rel_index_v = {0,0,0};
    rel_index_h = this->derivative_1_forward_2_index;
    coeff       = this->derivative_1_forward_2_coeff;
  }
  h0 = 0;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh;
  }

  // last column
  if( accuracy == "1" ){
    rel_index_v = {0,0};
    rel_index_h = this->derivative_1_backward_1_index;
    coeff       = this->derivative_1_backward_1_coeff;
  } else {
    rel_index_v = {0,0,0};
    rel_index_h = this->derivative_1_backward_2_index;
    coeff       = this->derivative_1_backward_2_coeff;
  }
  h0 = Nh-1;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh;
  }

  // middle chunk
  rel_index_v = {0,0,0};
  rel_index_h = this->derivative_1_central_2_index;
  coeff       = this->derivative_1_central_2_coeff;
  for(int v=0;v<Nv;v++){
    for(int h=1;h<Nh-1;h++){
      zout[v*Nh+h] = this->weighted_sum(v,h,rel_index_v,rel_index_h,coeff,Nh,zz)/dh;
    }
  }  
}


void RectGrid::calculate_derivative_2(int Nh,int Nv,double* h,double* zz,double* zout,std::string accuracy){
  // calculate the SECOND derivative along the horizontal axis
  int h0,v0;
  std::vector<int> rel_index_v;
  std::vector<int> rel_index_h;
  std::vector<double> coeff;
  double dh2 = pow(h[1]-h[0],2);
  
  // 1st column
  if( accuracy == "1" ){
    rel_index_v = {0,0,0};
    rel_index_h = this->derivative_2_forward_1_index;
    coeff       = this->derivative_2_forward_1_coeff;
  } else {
    rel_index_v = {0,0,0,0};
    rel_index_h = this->derivative_2_forward_2_index;
    coeff       = this->derivative_2_forward_2_coeff;
  }
  h0 = 0;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh2;
  }

  // last column
  if( accuracy == "1" ){
    rel_index_v = {0,0,0};
    rel_index_h = this->derivative_2_backward_1_index;
    coeff       = this->derivative_2_backward_1_coeff;
  } else {
    rel_index_v = {0,0,0,0};
    rel_index_h = this->derivative_2_backward_2_index;
    coeff       = this->derivative_2_backward_2_coeff;
  }
  h0 = Nh-1;
  for(int v=0;v<Nv;v++){
    zout[v*Nh+h0] = this->weighted_sum(v,h0,rel_index_v,rel_index_h,coeff,Nh,zz)/dh2;
  }

  // middle chunk
  rel_index_v = {0,0,0};
  rel_index_h = this->derivative_2_central_2_index;
  coeff       = this->derivative_2_central_2_coeff;
  for(int v=0;v<Nv;v++){
    for(int h=1;h<Nh-1;h++){
      zout[v*Nh+h] = this->weighted_sum(v,h,rel_index_v,rel_index_h,coeff,Nh,zz)/dh2;
    }
  }  
}

double RectGrid::weighted_sum(int i0,int j0,std::vector<int> rel_i,std::vector<int> rel_j,std::vector<double> coeff,int z_Nx,double* zz){
  double sum = 0.0;
  for(int k=0;k<rel_i.size();k++){
    int index_z  = (i0+rel_i[k])*z_Nx + (j0+rel_j[k]);
    sum += coeff[k]*zz[index_z];
  }
  return sum;
}
