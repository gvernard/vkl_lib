#ifndef SOURCE_PLANE_HPP
#define SOURCE_PLANE_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "tableDefinition.hpp"
#include "rectGrid.hpp"

class BaseCovKernel;

class BaseSourcePlane {
public:
  std::string source_type;         // adaptive or regular
  int Sm;                          // total number of pixels
  BaseCovKernel* kernel = NULL;    // pointer to kernel class
  std::map<std::string,mytable> H; // list of regularization matrices H for different schemes
  int eigenSparseMemoryAllocForH;  // estimate of the non-zero elements per row of the regularization matrix H

  BaseSourcePlane(){};
  BaseSourcePlane(BaseCovKernel* kernel);
  BaseSourcePlane(const BaseSourcePlane& other);
  ~BaseSourcePlane(){};
  
  virtual BaseSourcePlane* clone() = 0;
  virtual void constructH(const std::string reg_scheme) = 0;
  virtual void outputSource(const std::string path) = 0;
  virtual void outputSourceErrors(double* errors,const std::string path) = 0;

  void set_kernel(BaseCovKernel* kernel);
  void clear_H(std::string reg_scheme="all");
  bool find_reg(std::string reg_scheme);
  void print_reg();
};


class FixedSource: public BaseSourcePlane,public RectGrid {
public:
  FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,BaseCovKernel* kernel=NULL);
  FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath,BaseCovKernel* kernel=NULL);
  FixedSource(const FixedSource& other);
  ~FixedSource(){};
  
  virtual FixedSource* clone();
  virtual void constructH(std::string reg_scheme);
  virtual void outputSource(const std::string path);
  virtual void outputSourceErrors(double* errors,const std::string path){};

private:
  std::vector<mytriplet> getChunk(int i0,int j0,double hx,double hy,std::vector<int> ind_x,std::vector<double> coeff_x,std::vector<int> ind_y,std::vector<double> coeff_y);
  void appendToH(mytable& H,std::vector<mytriplet>& chunk);
};

/*
class AdaptiveSource: public BaseSourcePlane {
public:
  AdaptiveSource(int Sm,std::string reg_scheme);
  AdaptiveSource(std::string mode,int Sm,int spacing,std::string reg_scheme);
  AdaptiveSource(const AdaptiveSource& source);
  ~AdaptiveSource();

  virtual AdaptiveSource* clone();
  virtual void constructH();
  virtual void outputSource(const std::string path);
  virtual void outputSourceErrors(double* errors,const std::string path){};

  void constructDerivatives();
  void createAdaGrid(ImagePlane* image,CollectionMassModels* mycollection);
  void createDelaunay();
  void writeTriangles();
  void setGrid(std::map<std::string,std::string> pars){};
  bool pointInPolygon(double x,double y){};
  void writeVertexFaces(int index,std::string output);
  //  void boundPolygon();

private:

  struct xypoint {
    double x;
    double y;
  };
  
  struct a_triangle {
    int a;
    int b;
    int c;
  };

  xypoint intersection_point_x(xypoint p0,xypoint p1,xypoint p2);
  xypoint intersection_point_y(xypoint p0,xypoint p1,xypoint p2);
  double triangleArea(a_triangle triangle);

  std::string mode;
  int spacing;
  int n_triangles;
  std::vector<a_triangle> triangles;
  std::vector< std::vector<int> > opposite_edges_per_vertex;

  std::vector<int> mask_triangles;
};
*/


class FactorySourcePlane{//This is a singleton class.
public:
  FactorySourcePlane(FactorySourcePlane const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactorySourcePlane const&) = delete;

  static FactorySourcePlane* getInstance(){
    static FactorySourcePlane dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseSourcePlane* createSourcePlane(std::map<std::string,std::string> source){
    if( source["source_type"] == "fixed" ){
      return new FixedSource(stoi(source["Nx"]),stoi(source["Ny"]),stof(source["xmin"]),stof(source["xmax"]),stof(source["ymin"]),stof(source["xmax"]));
    } else if( source["source_type"] == "adaptive" ){
      //      if( source["mode"] == "random" ){
      //	return new AdaptiveSource(stoi(source["sm"]),source["reg_s"]);
      //      } else if( source["mode"] == "image" || source["mode"] == "grid" ){
      //	return new AdaptiveSource(source["mode"],stoi(source["sm"]),stoi(source["spacing"]),source["reg_s"]);
      //      }
    } else {
      return NULL;
    }
  }

private:
  FactorySourcePlane(){};
};


#endif /* SOURCE_PLANE_HPP */
