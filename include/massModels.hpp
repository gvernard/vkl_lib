#ifndef MASS_MODELS_HPP
#define MASS_MODELS_HPP

#include <cmath>
#include <vector>
#include <string>
#include <map>

#include "imagePlane.hpp"
#include "sourcePlane.hpp"

class BaseMassModel{
public:
  int Npars;
  std::string mass_type;
  std::map<std::string,double> mpars;

  BaseMassModel(){};
  ~BaseMassModel(){};
  
  virtual void updateMassPars(std::map<std::string,double> pars) = 0;   // for parametric mass models
  virtual void updateMassPars(std::string way,double* new_dpsi) = 0;      // for Pert mass model
  virtual void defl(double xin,double yin,double& xout,double& yout) = 0;
  virtual double kappa(double xin,double yin)   = 0;
  virtual void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi) = 0;
  virtual double psi(double xin,double yin) = 0;

  void printMassPars();

protected:
  BaseMassModel(int Npars,std::string mass_type): Npars(Npars), mass_type(mass_type){};
};

class CollectionMassModels {
public:
  std::vector<BaseMassModel*> models;
  
  CollectionMassModels(){};
  ~CollectionMassModels(){};
  void all_defl(double xin,double yin,double& xout,double& yout);
  double all_kappa(double xin,double yin);
  void all_gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double all_psi(double xin,double yin);
  double detJacobian(double xin,double yin);
};

class ExternalShear: public BaseMassModel{
public:
  ExternalShear(std::map<std::string,double> pars);
  void updateMassPars(std::map<std::string,double> pars);
  void updateMassPars(std::string way,double* new_dpsi){};
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin);
};

class Sie: public BaseMassModel{
public:
  Sie(std::map<std::string,double> pars);
  void updateMassPars(std::map<std::string,double> pars);
  void updateMassPars(std::string way,double* new_dpsi){};
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin);
private:
  void check_close_to_origin(double& x,double& y);
};

class Spemd: public BaseMassModel{
public:
  Spemd(std::map<std::string,double> pars);
  void updateMassPars(std::map<std::string,double> pars);
  void updateMassPars(std::string way,double* new_dpsi){};
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin);
};

class Pert: public BaseMassModel,public FixedSource {
public:
  Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);
  Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath);
  Pert(int Nx,int Ny,ImagePlane* image);
  ~Pert(){};
  void updateMassPars(std::map<std::string,double> pars){};
  void updateMassPars(std::string way,double* new_dpsi);
  void defl(double xin,double yin,double& xout,double& yout);
  double kappa(double xin,double yin);
  void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
  double psi(double xin,double yin);

  void updateDerivatives();
};


class FactoryParametricMassModel{//This is a singleton class.
public:
  FactoryParametricMassModel(FactoryParametricMassModel const&) = delete;//Stop the compiler generating methods of copy the object.
  void operator=(FactoryParametricMassModel const&) = delete;

  static FactoryParametricMassModel* getInstance(){
    static FactoryParametricMassModel dum;//Guaranteed to be destroyed. Instantiated on first call.
    return &dum;
  }

  // This creates only parametric models. Pert models have to be created manually beacuse they take different constructor options
  BaseMassModel* createParametricMassModel(const std::string modelname,std::map<std::string,double> pars){
    if( modelname == "sie" ){
      return new Sie(pars);
    } else if( modelname == "spemd" ){
      return new Spemd(pars);
    } else if( modelname == "external_shear" ){
      return new ExternalShear(pars);
    } else if( modelname == "pert" ){
      // throw exception
    } else {
      return NULL;
    }
  }

private:
  FactoryParametricMassModel(){};
};

#endif /* MASS_MODELS_HPP */
