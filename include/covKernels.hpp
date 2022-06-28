#ifndef COVARIANCE_KERNELS_HPP
#define COVARIANCE_KERNELS_HPP

#include <string>
#include <map> 

namespace vkl {

  class BaseCovKernel {
  public:
    int Npars;
    std::string cov_type;
    std::map<std::string,double> cpars;
  
    BaseCovKernel(){};
    BaseCovKernel(const BaseCovKernel& other);
    ~BaseCovKernel(){};
  
    virtual void updateCovPars(std::map<std::string,double> pars) = 0;
    virtual double getCovariance(double r) = 0;
    virtual double getCovarianceSelf() = 0;

    void printCovPars();
    bool larger_than_cmax(double cov);
  protected:
    BaseCovKernel(int Npars,std::string cov_type): Npars(Npars), cov_type(cov_type){};
  };


  class GaussKernel: public BaseCovKernel {
  public:
    GaussKernel(std::map<std::string,double> pars);
    double getCovariance(double r);
    double getCovarianceSelf();
    void updateCovPars(std::map<std::string,double> pars);
  };

  class ExpKernel: public BaseCovKernel {
  public:
    ExpKernel(std::map<std::string,double> pars);
    double getCovariance(double r);
    double getCovarianceSelf();
    void updateCovPars(std::map<std::string,double> pars);
  };

  class ExpGaussKernel: public BaseCovKernel {
  public:
    ExpGaussKernel(std::map<std::string,double> pars);
    double getCovariance(double r);
    double getCovarianceSelf();
    void updateCovPars(std::map<std::string,double> pars);
  };

  class FactoryCovKernel {//This is a singleton class.
  public:
    FactoryCovKernel(FactoryCovKernel const&) = delete;//Stop the compiler generating methods of copy the object.
    void operator=(FactoryCovKernel const&) = delete;

    static FactoryCovKernel* getInstance(){
      static FactoryCovKernel dum;//Guaranteed to be destroyed. Instantiated on first call.
      return &dum;
    }

    BaseCovKernel* createCovKernel(const std::string kernel_type,std::map<std::string,double> pars){
      if( kernel_type == "exp" ){
	return new ExpKernel(pars);
      } else if( kernel_type == "gauss" ){
	return new GaussKernel(pars);
      } else if( kernel_type == "exp_gauss" ){
	return new ExpGaussKernel(pars);
      } else {
	// throw exception
	return NULL;
      }
    }

  private:
    FactoryCovKernel(){};
  };

}

#endif /* COVARIANCE_KERNELS_HPP */
