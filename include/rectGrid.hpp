#ifndef RECT_GRID_HPP
#define RECT_GRID_HPP

#include <cstdlib>
#include <string>
#include <valarray>

class RectGrid{
public:
  int Nx;
  int Ny;
  int Nz;           // Total number of pixels in the grid
  double* center_x; // 1D, size of Nx: pixel center x coordinates
  double* center_y; // 1D, size of Ny: pixel center y coordinates
  double* z;        // 1D, size of Nx*Ny: Values of the 2D surface in row-major format, i.e. z[i*Nx+j]

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

private:
  void common_constructor(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);
  void readFits(const std::string filepath,std::valarray<float>& contents);
};

#endif /* RECT_GRID_HPP */
