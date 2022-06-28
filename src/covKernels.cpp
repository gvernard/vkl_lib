#include "covKernels.hpp"

#include <cmath>

using namespace vkl;

//Abstract class: BaseCovKernel
//===============================================================================================================
BaseCovKernel::BaseCovKernel(const BaseCovKernel& other){
  this->Npars = other.Npars;
  this->cov_type = other.cov_type;
  for(std::map<std::string,double>::iterator it=cpars.begin();it!=cpars.end();it++){
    this->cpars[it->first] = it->second;
  }
};
void BaseCovKernel::printCovPars(){
  printf("      %10s\n",this->cov_type.c_str());
  for(std::map<std::string,double>::iterator it=cpars.begin();it!=cpars.end();it++){
    printf("%s: %10.5f\n",it->first.c_str(),it->second);
  }
}
bool BaseCovKernel::larger_than_cmax(double cov){
  if( cov > cpars["cmax"] ){
    return true;
  } else {
    return false;
  }
}


//Derived class: ExpKernel
//===============================================================================================================
ExpKernel::ExpKernel(std::map<std::string,double> pars): BaseCovKernel(2,"exp"){
  updateCovPars(pars);
}
void ExpKernel::updateCovPars(std::map<std::string,double> pars){
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    cpars[it->first] = it->second;
  }
}
double ExpKernel::getCovariance(double r){
  return exp(-r/cpars["sdev"]);
}
double ExpKernel::getCovarianceSelf(){
  return 1.0;
}


//Derived class: GaussKernel
//===============================================================================================================
GaussKernel::GaussKernel(std::map<std::string,double> pars): BaseCovKernel(2,"gauss"){
  updateCovPars(pars);
}
void GaussKernel::updateCovPars(std::map<std::string,double> pars){
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    cpars[it->first] = it->second;
  }
}
double GaussKernel::getCovariance(double r){
  return exp(-r*r/(2*pow(cpars["sdev"],2)));
}
double GaussKernel::getCovarianceSelf(){
  return 1.1;
}


//Derived class: ExpGaussKernel
//===============================================================================================================
ExpGaussKernel::ExpGaussKernel(std::map<std::string,double> pars): BaseCovKernel(2,"exp_gauss"){
  updateCovPars(pars);
}
void ExpGaussKernel::updateCovPars(std::map<std::string,double> pars){
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    cpars[it->first] = it->second;
  }
  cpars["fac"] = 1.0/(cpars["sdev"]*sqrt(2*M_PI));
}
double ExpGaussKernel::getCovariance(double r){
  double cov = cpars["fac"]*exp(-pow(r,cpars["expo"])/(2*pow(cpars["sdev"],2)));
  return cov;
}
double ExpGaussKernel::getCovarianceSelf(){
  return cpars["fac"];
}


//Derived class: INSERT YOUR CUSTOM COVARIANCE CLASS IMPLEMENTATION HERE
//===============================================================================================================
