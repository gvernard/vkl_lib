#ifndef SOURCE_PROFILE_HPP
#define SOURCE_PROFILE_HPP

#include <cmath>
#include <string>
#include <map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "nonLinearPars.hpp"
#include "imagePlane.hpp"
#include "rectGrid.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

class BaseProfile {
public:
  std::string type;
  int output_res;
  
  BaseProfile(){};
  ~BaseProfile(){};

  virtual double value(double x,double y) = 0;
  virtual void outputProfile(std::string filename) = 0;
  virtual void scaleProfile() = 0;
  virtual void inverseScaleProfile() = 0;
  
  void profile(int Sj,int Si,double* sx,double* sy,double* s);
  void writeProfile(std::string filename,double half_range);
};




class BaseAnalyticFunction {
public:
  std::string type;
  std::map<std::string,double> pars;

  BaseAnalyticFunction(){};
  ~BaseAnalyticFunction(){};

  virtual double function_value(double x,double y) = 0;
  virtual std::vector<double> extent() = 0;
};


class Sersic: public BaseAnalyticFunction {
public:
  Sersic(std::map<std::string,double> pars);
  double function_value(double x,double y);
  std::vector<double> extent();
  void scaleProfile();
  void inverseScaleProfile();
private:
  const double fac = 0.01745329251; // conversion from degrees to rad
};


class proGauss: public BaseAnalyticFunction { // name Gauss is taken in nonLinearPars.hpp so I use proGauss (pro for profile)
public:
  proGauss(std::map<std::string,double> pars);
  double function_value(double x,double y);
  std::vector<double> extent();
  void scaleProfile();
  void inverseScaleProfile(){};
private:
  const double fac = 0.01745329251; // conversion from degrees to rad
};


class FactoryAnalyticFunction {
public:
  FactoryAnalyticFunction(FactoryAnalyticFunction const&) = delete;
  void operator=(FactoryAnalyticFunction const&) = delete;

  static FactoryAnalyticFunction* getInstance(){
    static FactoryAnalyticFunction dum;
    return &dum;
  }

  BaseAnalyticFunction* createAnalyticFunction(const std::string &name,std::map<std::string,double> pars){
    if( name == "sersic" ){
      return new Sersic(pars);
    } else if ( name == "gauss" ){
      return new proGauss(pars);
    } else {
      return NULL;
    }
  }

private:
  FactoryAnalyticFunction(){};
};


class Analytic: public BaseProfile {
public:
  std::vector<BaseAnalyticFunction*> components;

  Analytic(std::vector<std::string> names,std::vector<std::map<std::string,double> > all_pars);
  ~Analytic(){
    for(int i=0;i<this->components.size();i++){
      delete(components[i]);
    }
  }

  double value(double x,double y);
  void outputProfile(std::string filename);
  void scaleProfile(){};
  void inverseScaleProfile(){};
};





class myDelaunay: public BaseProfile {
public:
  int N;
  double* x;
  double* y;
  double* src;

  myDelaunay(std::string filename);
  ~myDelaunay(){
    free(x);
    free(y);
    free(src);
    free(convex_hull);
  }
  double value(double x,double y);
  void outputProfile(std::string filename);
  void scaleProfile(){};
  void inverseScaleProfile(){};
  double sourceExtent();
  
private:
  struct atriangle {
    int a;
    int b;
    int c;
  };
  std::vector<atriangle> triangles;
  Point* convex_hull;
  int ch_size;
};



class Custom: public BaseProfile,public RectGrid {
public:
  Custom(std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,double Mtot,std::string interp): BaseProfile(),RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
    this->type = "custom";
    this->Mtot = Mtot;
    scaleProfile();
    this->set_interp(interp);
  }
  ~Custom(){}

  double Mtot;

  double value(double x,double y);
  void outputProfile(std::string filename);
  void scaleProfile();
  void inverseScaleProfile(){};
};


#endif /* SOURCE_PROFILE_HPP */
