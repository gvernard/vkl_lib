#ifndef IRRE_GRID_HPP
#define IRRE_GRID_HPP

#include <string>
#include <vector>
#include <map>

namespace vkl {
  
  class Triangle{
  public:
    double x[3];
    double y[3];
    double z[3];

    Triangle(double* x,double* y);
    Triangle(double* x,double* y,double* z);
    bool interp_bary(double x,double y,double& value);
    double triangle_area();
  };

  class Cell{
  public:
    int Nc;
    double* x = NULL;
    double* y = NULL;
    double xmin;
    double xmax;
    double ymin;
    double ymax;

    bool point_in_cell(double x,double y);
    double area();

  private:
    void find_bounding_box();
    bool point_in_bounding_box(double x,double y);
    bool sides_are_intersecting(double v1x1, double v1y1, double v1x2, double v1y2,double v2x1, double v2y1, double v2x2, double v2y2);
  }


    class IrreGrid{
    public:
      int Nz;
      double* center_x;
      double* center_y;
      double* z  = NULL;
      double* zx = NULL;
      double* zy = NULL;
      Cell boundary; 
    };

}



#endif /* IRRE_GRID_HPP */
