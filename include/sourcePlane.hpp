#ifndef SOURCE_PLANE_HPP
#define SOURCE_PLANE_HPP

#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <iostream>

#include "covKernels.hpp"
#include "tableDefinition.hpp"


class BaseSourcePlane {
public:
  std::string type;                // adaptive or regular
  int Sm;                          // total number of pixels
  std::string reg;                 // name of the regularization scheme
  int eigenSparseMemoryAllocForH;  // estimate of the non-zero elements per row of the regularization matrix H
  BaseCovKernel* kernel;           // pointer to kernel class
  mytable H;

  BaseSourcePlane(){};
  BaseSourcePlane(const BaseSourcePlane& other){
    type = other.type;
    Sm   = other.Sm;
    reg  = other.reg;
    H.Ti = other.H.Ti;
    H.Tj = other.H.Tj;
    eigenSparseMemoryAllocForH = other.eigenSparseMemoryAllocForH;
  };
  ~BaseSourcePlane(){
    if( this->reg == "covariance_kernel" ){
      delete this->kernel;
    }
  };
  
  virtual BaseSourcePlane* clone() = 0;
  virtual void constructH() = 0;
  virtual void outputSource(const std::string path) = 0;
  virtual void outputSourceErrors(double* errors,const std::string path) = 0;
};


class FixedSource: public BaseSourcePlane {
public:
  RectGrid* grid;

  FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string reg_scheme);
  FixedSource(const FixedSource& source) : BaseSourcePlane(source) {};
  virtual FixedSource* clone(){
    return new FixedSource(*this);
  };
  ~FixedSource(){
    delete(grid);
  }
  
  virtual void constructH();
  virtual void outputSource(const std::string path);
  virtual void outputSourceErrors(double* errors,const std::string path){};
};


class AdaptiveSource: public BaseSourcePlane {
public:
  AdaptiveSource(int Sm,std::string reg_scheme);
  AdaptiveSource(std::string mode,int Sm,int spacing,std::string reg_scheme);
  AdaptiveSource(const AdaptiveSource& source);
  virtual AdaptiveSource* clone(){
    return new AdaptiveSource(*this);
  };
  ~AdaptiveSource();

  virtual void constructH();
  virtual void outputSource(const std::string path);
  virtual void outputSourceErrors(double* errors,const std::string path);

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



class FactorySourcePlane{//This is a singleton class.
public:
  FactorySourcePlane(FactorySourcePlane const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactorySourcePlane const&) = delete;

  static FactorySourcePlane* getInstance(){
    static FactorySourcePlane dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseSourcePlane* createSourcePlane(std::map<std::string,std::string> source){
    if( source["type"] == "fixed" ){
      return new FixedSource(stoi(source["sx"]),stoi(source["sy"]),stof(source["size"]),source["reg_s"]);
      //    } else if( source["type"] == "floating" ){
      //      return new FloatingSource(stoi(source["sx"]),stoi(source["sy"]),stof(source["size"]),stof(source["x0"]),stof(source["y0"]),source["reg"]);
    } else if( source["type"] == "adaptive" ){
      if( source["mode"] == "random" ){
	return new AdaptiveSource(stoi(source["sm"]),source["reg_s"]);
      } else if( source["mode"] == "image" || source["mode"] == "grid" ){
	return new AdaptiveSource(source["mode"],stoi(source["sm"]),stoi(source["spacing"]),source["reg_s"]);
      }
    } else {
      return NULL;
    }
  }

private:
  FactorySourcePlane(){};
};


#endif /* SOURCE_PLANE_HPP */
