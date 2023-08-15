#ifndef MASS_MODELS_HPP
#define MASS_MODELS_HPP

#include <cmath>
#include <vector>
#include <string>
#include <map>

#include "imagePlane.hpp"
#include "sourcePlane.hpp"

namespace vkl {

  class BaseMassModel{
  public:
    int Npars;
    std::string mass_type;
    std::map<std::string,double> mpars;

    BaseMassModel(){};
    BaseMassModel(const BaseMassModel& other);
    ~BaseMassModel(){};
  
    virtual void updateMassPars(std::map<std::string,double> pars) = 0;   // for parametric mass models
    virtual void updateMassPars(std::string way,double* new_dpsi) = 0;      // for Pert mass model
    virtual void defl(double xin,double yin,double& xout,double& yout) = 0;
    virtual double kappa(double xin,double yin)   = 0;
    virtual void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi) = 0;
    virtual double psi(double xin,double yin) = 0;
    virtual void getExtent(double& xmin,double& xmax,double& ymin,double& ymax) = 0;
  
    void printMassPars();

  protected:
    BaseMassModel(int Npars,std::string mass_type): Npars(Npars), mass_type(mass_type){};
  };

  class CollectionMassModels {
  public:
    std::vector<BaseMassModel*> models;
  
    CollectionMassModels(){};
    CollectionMassModels(const CollectionMassModels& other);
    ~CollectionMassModels();
    void all_defl(double xin,double yin,double& xout,double& yout);
    double all_kappa(double xin,double yin);
    void all_gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
    double all_psi(double xin,double yin);
    double detJacobian(double xin,double yin);
    void getExtent(double& xmin,double& xmax,double& ymin,double& ymax);
  };

  class ExternalShear: public BaseMassModel{
  public:
    ExternalShear(std::map<std::string,double> pars);
    ExternalShear(const ExternalShear& other) : BaseMassModel(other){};
    void updateMassPars(std::map<std::string,double> pars);
    void updateMassPars(std::string way,double* new_dpsi){};
    void defl(double xin,double yin,double& xout,double& yout);
    double kappa(double xin,double yin);
    void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
    double psi(double xin,double yin);
    void getExtent(double& xmin,double& xmax,double& ymin,double& ymax);
  };

  class Sie: public BaseMassModel{
  public:
    Sie(std::map<std::string,double> pars);
    Sie(const Sie& other) : BaseMassModel(other){}
    void updateMassPars(std::map<std::string,double> pars);
    void updateMassPars(std::string way,double* new_dpsi){};
    void defl(double xin,double yin,double& xout,double& yout);
    double kappa(double xin,double yin);
    void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
    double psi(double xin,double yin);
    void getExtent(double& xmin,double& xmax,double& ymin,double& ymax);
  private:
    void check_close_to_origin(double& x,double& y);
  };

  class Spemd: public BaseMassModel{
  public:
    Spemd(std::map<std::string,double> pars);
    Spemd(const Spemd& other) : BaseMassModel(other){}
    void updateMassPars(std::map<std::string,double> pars);
    void updateMassPars(std::string way,double* new_dpsi){};
    void defl(double xin,double yin,double& xout,double& yout);
    double kappa(double xin,double yin);
    void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
    double psi(double xin,double yin);
    void getExtent(double& xmin,double& xmax,double& ymin,double& ymax);
  };

  class Pert: public BaseMassModel,public FixedSource {
  public:
    Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax);
    Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath);
    Pert(int Nx,int Ny,ImagePlane* image);
    Pert(const Pert& other) : BaseMassModel(other),FixedSource(other){}
    ~Pert(){};
    void updateMassPars(std::map<std::string,double> pars){};
    void updateMassPars(std::string way,double* new_dpsi);
    void defl(double xin,double yin,double& xout,double& yout);
    double kappa(double xin,double yin);
    void gamma(double xin,double yin,double& gamma_mag,double& gamma_phi);
    double psi(double xin,double yin);
    void getExtent(double& xmin,double& xmax,double& ymin,double& ymax);

    void updateDerivatives();
  };


  class FactoryMassModel{
  public:
    FactoryMassModel(FactoryMassModel const&) = delete;
    void operator=(FactoryMassModel const&) = delete;

    static FactoryMassModel* getInstance(){
      static FactoryMassModel dum;
      return &dum;
    }

    // This creates only parametric models. Pert models have to be created manually beacuse they take different constructor options
    BaseMassModel* createMassModel(const std::string model_name,std::map<std::string,std::string> pars){
      if( model_name == "sie" ){
	std::map<std::string,double> tmp_pars;
	for(std::map<std::string,std::string>::iterator it=pars.begin();it!=pars.end();it++){
	  tmp_pars.insert( std::pair<std::string,double>(it->first,std::stod(it->second)) );
	}
	return new Sie(tmp_pars);
      } else if( model_name == "spemd" ){
	std::map<std::string,double> tmp_pars;
	for(std::map<std::string,std::string>::iterator it=pars.begin();it!=pars.end();it++){
	  tmp_pars.insert( std::pair<std::string,double>(it->first,std::stod(it->second)) );
	}
	return new Spemd(tmp_pars);
      } else if( model_name == "external_shear" ){
	std::map<std::string,double> tmp_pars;
	for(std::map<std::string,std::string>::iterator it=pars.begin();it!=pars.end();it++){
	  tmp_pars.insert( std::pair<std::string,double>(it->first,std::stod(it->second)) );
	}
	return new ExternalShear(tmp_pars);
      } else if( model_name == "pert" ){
	std::string filepath = pars["filepath"];
	int Nx      = std::stoi(pars["Nx"]);
	int Ny      = std::stoi(pars["Ny"]);
	double xmin = std::stod(pars["xmin"]);
	double xmax = std::stod(pars["xmax"]);
	double ymin = std::stod(pars["ymin"]);
	double ymax = std::stod(pars["ymax"]);
	return new Pert(Nx,Ny,xmin,xmax,ymin,ymax,filepath);
      } else {
	return NULL;
      }
    }

    BaseMassModel* copyMassModel(BaseMassModel* other){
      if( other->mass_type == "sie" ){
	Sie* ptr = static_cast<Sie*> (other);
	return new Sie(*ptr);
      } else if( other->mass_type == "spemd" ){
	Spemd* ptr = static_cast<Spemd*> (other);
	return new Spemd(*ptr);
      } else if( other->mass_type == "external_shear" ){
	ExternalShear* ptr = static_cast<ExternalShear*> (other);
	return new ExternalShear(*ptr);
      } else if( other->mass_type == "pert" ){
	Pert* ptr = static_cast<Pert*> (other);
        return new Pert(*ptr);
      } else {
	return NULL;
      }
    }
    
  private:
    FactoryMassModel(){};
  };

}
  
#endif /* MASS_MODELS_HPP */
