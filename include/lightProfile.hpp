#ifndef SOURCE_PROFILE_HPP
#define SOURCE_PROFILE_HPP

#include <string>
#include <map>
#include <vector>

#include "rectGrid.hpp"

namespace vkl {

  //#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
  //#include <CGAL/Polygon_2_algorithms.h>
  //typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  //typedef K::Point_2 Point;

  class BaseProfile {
  public:
    int Npars;
    std::string profile_type;
    std::map<std::string,double> ppars;

    BaseProfile(){};
    BaseProfile(const BaseProfile& other);
    ~BaseProfile(){};

    virtual void updateProfilePars(std::map<std::string,double> pars) = 0;
    virtual double value(double x,double y) = 0;
    virtual bool is_in_range(double xin,double yin) = 0;
    virtual void get_extent(double& xmin,double& xmax,double& ymin,double& ymax) = 0;
  
  protected:
    BaseProfile(int Npars,std::string profile_type): Npars(Npars), profile_type(profile_type){};
  };


  class CollectionProfiles {
  public:
    std::vector<BaseProfile*> profiles;
  
    CollectionProfiles(){};
    CollectionProfiles(const CollectionProfiles& other);
    ~CollectionProfiles();
    double all_values(double xin,double yin);
    void getExtent(double& xmin,double& xmax,double& ymin,double& ymax);
    void write_all_profiles(const std::string filepath);
  };



  class Sersic: public BaseProfile {
  public:
    Sersic(std::map<std::string,double> pars);
    Sersic(const Sersic& other);
    ~Sersic(){};
    void updateProfilePars(std::map<std::string,double> pars);
    double value(double x,double y);
    bool is_in_range(double xin,double yin);
    void get_extent(double& xmin,double& xmax,double& ymin,double& ymax);
  private:
    void set_extent();
    double p_xmin;
    double p_xmax;
    double p_ymin;
    double p_ymax;  
    double bn;
    double cospa;
    double sinpa;
    double rr;
  };


  class Gauss: public BaseProfile {
  public:
    Gauss(std::map<std::string,double> pars);
    Gauss(const Gauss& other);
    ~Gauss(){};
    void updateProfilePars(std::map<std::string,double> pars);
    double value(double x,double y);
    bool is_in_range(double xin,double yin);
    void get_extent(double& xmin,double& xmax,double& ymin,double& ymax);
  private:
    void set_extent();
    double p_xmin;
    double p_xmax;
    double p_ymin;
    double p_ymax;  
    double sdev;
    double cospa;
    double sinpa;
  };


  class Custom: public BaseProfile,public RectGrid {
  public:
    Custom(std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,double Mtot,std::string interp);
    Custom(const Custom& other);  
    ~Custom(){};
    void updateProfilePars(std::map<std::string,double> pars){};
    double value(double x,double y);
    bool is_in_range(double xin,double yin);
    void get_extent(double& xmin,double& xmax,double& ymin,double& ymax);
  private:
    double Mtot;
    void scaleProfile();
  };

  /*
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
  */

  class FactoryProfile {
  public:
    FactoryProfile(FactoryProfile const&) = delete;
    void operator=(FactoryProfile const&) = delete;

    static FactoryProfile* getInstance(){
      static FactoryProfile dum;
      return &dum;
    }

    BaseProfile* createProfile(const std::string profile_name,std::map<std::string,std::string> pars){
      if( profile_name == "sersic" ){
	std::map<std::string,double> tmp_pars;
	for(std::map<std::string,std::string>::iterator it=pars.begin();it!=pars.end();it++){
	  tmp_pars.insert( std::pair<std::string,double>(it->first,std::stod(it->second)) );
	}
	return new Sersic(tmp_pars);
      } else if( profile_name == "gauss" ){
	std::map<std::string,double> tmp_pars;
	for(std::map<std::string,std::string>::iterator it=pars.begin();it!=pars.end();it++){
	  tmp_pars.insert( std::pair<std::string,double>(it->first,std::stod(it->second)) );
	}
	return new Gauss(tmp_pars);
      } else if( profile_name == "irregular" ){
	return NULL;
      } else if( profile_name == "custom" ){
	std::string filepath = pars["filepath"];
	std::string interp   = pars["interp"];
	int Nx      = std::stoi(pars["Nx"]);
	int Ny      = std::stoi(pars["Ny"]);
	double xmin = std::stod(pars["xmin"]);
	double xmax = std::stod(pars["xmax"]);
	double ymin = std::stod(pars["ymin"]);
	double ymax = std::stod(pars["ymax"]);
	if( pars.find("M_tot") != pars.end() ){
	  double Mtot = std::stod(pars["M_tot"]);
	  return new Custom(filepath,Nx,Ny,xmin,xmax,ymin,ymax,Mtot,interp);
	} else {
	  return new Custom(filepath,Nx,Ny,xmin,xmax,ymin,ymax,0.0,interp);
	}
      } else {
	return NULL;
      }
    }

  private:
    FactoryProfile(){};
  };

}

#endif /* SOURCE_PROFILE_HPP */
