#ifndef SOURCE_PLANE_HPP
#define SOURCE_PLANE_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "covKernels.hpp"
#include "tableDefinition.hpp"
#include "rectGrid.hpp"

class BaseSourcePlane {
public:
  std::string source_type;         // adaptive or regular
  int Sm;                          // total number of pixels
  int eigenSparseMemoryAllocForH;  // estimate of the non-zero elements per row of the regularization matrix H
  BaseCovKernel* kernel;           // pointer to kernel class
  std::map<std::string,mytable> H;

  BaseSourcePlane(){};
  BaseSourcePlane(const BaseSourcePlane& other){
    source_type = other.source_type;
    Sm   = other.Sm;
    H    = other.H;
    eigenSparseMemoryAllocForH = other.eigenSparseMemoryAllocForH;
  };
  ~BaseSourcePlane(){
    if( this->find_reg("covariance_kernel") ){
      delete this->kernel;
    }
  };
  
  virtual BaseSourcePlane* clone() = 0;
  virtual void constructH(const std::string reg_scheme) = 0;
  virtual void outputSource(const std::string path) = 0;
  virtual void outputSourceErrors(double* errors,const std::string path) = 0;

  void clear_H(std::string reg_scheme="all"){
    if( reg_scheme == "all" ){
      this->H.clear();
    } else {
      if( !this->find_reg(reg_scheme) ){
	this->H.erase(reg_scheme);
      }
    }
  }

  bool find_reg(std::string reg_scheme){
    if( this->H.find(reg_scheme) == this->H.end() ){
      return false;
    } else {
      return true;
    }
  }

  void print_reg(){
    for(std::map<std::string,mytable>::iterator it=this->H.begin();it!=this->H.end();it++){
      std::cout << it->first << std::endl;
    }
  }
};


class FixedSource: public BaseSourcePlane,public RectGrid {
public:
  FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax): RectGrid(Nx,Ny,xmin,xmax,ymin,ymax){
    source_type = "fixed";
    Sm   = this->Nz;
  }
  FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath): RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
    source_type = "fixed";
    Sm   = this->Nz;
  }
  FixedSource(const FixedSource& other): BaseSourcePlane(other), RectGrid(other) {
    source_type = "fixed";
    Sm   = other.Sm;
  };
  virtual FixedSource* clone(){
    return new FixedSource(*this);
  };
  ~FixedSource(){};
  
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
