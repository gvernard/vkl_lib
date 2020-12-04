#ifndef MASS_MODELS_HPP
#define MASS_MODELS_HPP

#include <cmath>
#include <vector>
#include <string>
#include <map>

#include "nonLinearPars.hpp"
#include "imagePlane.hpp"
#include "sourcePlane.hpp"

extern "C"{
  void fastelldefl_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* defl);
  void fastellmag_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* defl,double* jacob);
}


class BaseMassModel{
public:
  int n;
  std::string mass_type;
  std::map<std::string,double> mpars;

  BaseMassModel(){};
  ~BaseMassModel(){
    mpars.clear();
  };
  
  void setMassPars(std::vector<Nlpar*> nlpars);
  void printMassPars();
  virtual void defl(double xin,double yin,double& xout,double& yout) = 0;
  virtual double kappa(double xin,double yin)   = 0;
  virtual void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi) = 0;
  virtual double psi(double xin,double yin) = 0;
};


class Sie: public BaseMassModel{
public:
  Sie(std::vector<Nlpar*> nlpars);
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin);
};

class Spemd: public BaseMassModel{
public:
  Spemd(std::vector<Nlpar*> nlpars);
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin){};
};

class Pert: public BaseMassModel,public FixedSource {
public:
  Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax): FixedSource(Nx,Ny,xmin,xmax,ymin,ymax){
    this->mass_type = "pert";
  }
  Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath): FixedSource(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
    this->mass_type = "pert";
    this->updatePert();
  }
  Pert(int Nx,int Ny,ImagePlane* image): FixedSource(Nx,Ny,image->grid->xmin,image->grid->xmax,image->grid->ymin,image->grid->ymax){
    this->mass_type = "pert";
  }
  ~Pert(){};
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin);
  void replaceDpsi(double* new_dpsi);
  void addDpsi(double* corrections);
  void updatePert();
};


class FactoryMassModel{//This is a singleton class.
public:
  FactoryMassModel(FactoryMassModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryMassModel const&) = delete;

  static FactoryMassModel* getInstance(){
    static FactoryMassModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  BaseMassModel* createMassModel(const std::string &modelname,std::vector<Nlpar*> nlpars){
    if( modelname == "sie" ){
      return new Sie(nlpars);
    } else if ( modelname == "spemd" ){
      return new Spemd(nlpars);
    } else {
      return NULL;
    }
  }

  BaseMassModel* createMassModel(const std::string &modelname,std::vector<Nlpar*> nlpars,double Dls,double Dos){
    if( modelname == "sie" ){
      for(int i=0;i<nlpars.size();i++){
	if( nlpars[i]->nam == "sigma" ){
	  nlpars[i]->nam = "b";
	  nlpars[i]->val = 2.592*pow(nlpars[i]->val/299.792458,2)*Dls/Dos;
	}
      }
      return new Sie(nlpars);
    } else if ( modelname == "spemd" ){
      return new Spemd(nlpars);
    } else {
      return NULL;
    }
  }

  BaseMassModel* createMassModel(const std::string &modelname,std::map<std::string,std::string> pars){
    if( modelname == "pert" ){
      return new Pert(std::stoi(pars["Ni"]),std::stoi(pars["Nj"]),std::stof(pars["xmin"]),std::stof(pars["xmax"]),std::stof(pars["ymin"]),std::stof(pars["ymax"]),pars["filename"]);
    } else {
      return NULL;
    }
  }

private:
  FactoryMassModel(){};
};



class CollectionMassModels {
public:
  int n;
  std::vector<BaseMassModel*> models;
  std::map<std::string,double> mpars;
  
  CollectionMassModels();
  CollectionMassModels(std::vector<Nlpar*> nlpars);
  ~CollectionMassModels();
  void printPhysPars();
  void setPhysicalPars(std::vector<Nlpar*> nlpars);
  void all_defl(double xin,double yin,double& xout,double& yout);
  void all_defl(ImagePlane* image);
  double all_kappa(double xin,double yin);
  void all_kappa(ImagePlane* image,ImagePlane* kappa_tot);
  void all_gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  void all_gamma(ImagePlane* image,ImagePlane* gamma_mag,ImagePlane* gamma_phi);
  double all_psi(double xin,double yin);
  void all_psi(ImagePlane* image,ImagePlane* psi_tot);
  double detJacobian(double xin,double yin);
  void detJacobian(ImagePlane* image,ImagePlane* det);
};

#endif /* MASS_MODELS_HPP */
