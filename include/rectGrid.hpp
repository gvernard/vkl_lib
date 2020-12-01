#ifndef RECT_GRID_HPP
#define RECT_GRID_HPP

#include <string>
#include <vector>

struct triplet{
  int i;    // index in the y direction
  int j;    // index in the x direction
  double v; // value
};

class RectGrid{
public:
  int Nx;
  int Ny;
  int Nz;              // Total number of pixels in the grid
  double* center_x;    // 1D, size of Nx: pixel center x coordinates
  double* center_y;    // 1D, size of Ny: pixel center y coordinates
  double* z   = NULL;  // 1D, size of Nx*Ny: Values of the 2D surface in row-major format, i.e. z[i*Nx+j]
  double* zx  = NULL;  // same as z: first derivative along the x axis
  double* zy  = NULL;  // same as z: first derivative along the y axis
  double* zxy = NULL;  // same as z: mixed derivative along the x and y axes
  double* zxx = NULL;  // same as z: second derivative along the x axis
  double* zyy = NULL;  // same as z: second derivative along the y axis
  
  double* bound_x; // 1D, size of Nx+1: pixel x boundaries, can be accessed by x[j],x[j+1], where 0<j<Nx
  double* bound_y; // 1D, size of Ny+1: pixel y boundaries, can be accessed by y[i],y[i+1], where 0<i<Ny
  double width;
  double height;
  double step_x;
  double step_y;
  
  RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);
  RectGrid(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,const std::string filepath);
  ~RectGrid();

  bool point_in_grid(double x,double y);
  bool point_between_pixel_centers(double x,double y,int boundary_size);
  void match_point_to_pixel(double x,double y,int& i0,int& j0);

  bool match_point_to_closest_4(double x,double y,int* i,int* j);
  bool match_point_to_closest_16(double x,double y,int* i,int* j);

  void calculate_zx(std::string accuracy);
  void calculate_zy(std::string accuracy);
  void calculate_zxy(std::string accuracy);
  void calculate_zxx(std::string accuracy);
  void calculate_zyy(std::string accuracy);

  
private:
  void common_constructor(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);
  void readFits(double* field,const std::string filepath);
  void multiplyVectorByScalar(std::vector<double> &v,double k);
  
  void calculate_derivative_1(int Nh,int Nv,double* h,double* zz,double* zout,std::string accuracy);
  void calculate_derivative_2(int Nh,int Nv,double* h,double* zz,double* zout,std::string accuracy);

  // Finite difference coefficients and relative indices.
  // 1st derivatives
  static const std::vector<int>    derivative_1_forward_1_index;
  static const std::vector<double> derivative_1_forward_1_coeff;
  static const std::vector<int>    derivative_1_forward_2_index;
  static const std::vector<double> derivative_1_forward_2_coeff;
  static const std::vector<int>    derivative_1_backward_1_index;
  static const std::vector<double> derivative_1_backward_1_coeff;
  static const std::vector<int>    derivative_1_backward_2_index;
  static const std::vector<double> derivative_1_backward_2_coeff;
  static const std::vector<int>    derivative_1_central_2_index;
  static const std::vector<double> derivative_1_central_2_coeff;
  // 2nd derivatives
  static const std::vector<int>    derivative_2_forward_1_index;
  static const std::vector<double> derivative_2_forward_1_coeff;
  static const std::vector<int>    derivative_2_forward_2_index;
  static const std::vector<double> derivative_2_forward_2_coeff;
  static const std::vector<int>    derivative_2_backward_1_index;
  static const std::vector<double> derivative_2_backward_1_coeff;
  static const std::vector<int>    derivative_2_backward_2_index;
  static const std::vector<double> derivative_2_backward_2_coeff;
  static const std::vector<int>    derivative_2_central_2_index;
  static const std::vector<double> derivative_2_central_2_coeff;

  double weighted_sum(int i0,int j0,std::vector<int> rel_i,std::vector<int> rel_j,std::vector<double> coeff,int z_Nx,double* zz);
};

// 2nd derivatives
const std::vector<int>    RectGrid::derivative_1_forward_1_index({0,1});
const std::vector<double> RectGrid::derivative_1_forward_1_coeff({-1,1});
const std::vector<int>    RectGrid::derivative_1_forward_2_index({0,1,2});
const std::vector<double> RectGrid::derivative_1_forward_2_coeff({-1.5,2.0,-0.5});
const std::vector<int>    RectGrid::derivative_1_backward_1_index({-1,0});
const std::vector<double> RectGrid::derivative_1_backward_1_coeff({-1,1});
const std::vector<int>    RectGrid::derivative_1_backward_2_index({-2,-1,0});
const std::vector<double> RectGrid::derivative_1_backward_2_coeff({0.5,-2,1.5});
const std::vector<int>    RectGrid::derivative_1_central_2_index({-1,0,1});
const std::vector<double> RectGrid::derivative_1_central_2_coeff({-0.5,0.0,0.5});
// 2nd derivatives
const std::vector<int>    RectGrid::derivative_2_forward_1_index({0,1,2});
const std::vector<double> RectGrid::derivative_2_forward_1_coeff({1,-2,1});  
const std::vector<int>    RectGrid::derivative_2_forward_2_index({0,1,2,3});
const std::vector<double> RectGrid::derivative_2_forward_2_coeff({2,-5,4,-1});
const std::vector<int>    RectGrid::derivative_2_backward_1_index({-2,-1,0});
const std::vector<double> RectGrid::derivative_2_backward_1_coeff({1,-2,1});  
const std::vector<int>    RectGrid::derivative_2_backward_2_index({-3,-2,-1,0});
const std::vector<double> RectGrid::derivative_2_backward_2_coeff({-1,4,-5,2});
const std::vector<int>    RectGrid::derivative_2_central_2_index({-1,0,1});
const std::vector<double> RectGrid::derivative_2_central_2_coeff({1,-2,1});



#endif /* RECT_GRID_HPP */
