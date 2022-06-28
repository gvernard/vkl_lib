#include "sourcePlane.hpp"

using namespace vkl;


BaseSourcePlane:: BaseSourcePlane(BaseCovKernel* kernel){
  this->set_kernel(kernel);
}


BaseSourcePlane::BaseSourcePlane(const BaseSourcePlane& other){
  source_type = other.source_type;
  Sm     = other.Sm;
  kernel = other.kernel;
  H      = other.H;
  eigenSparseMemoryAllocForH = other.eigenSparseMemoryAllocForH;
};

void BaseSourcePlane::set_kernel(BaseCovKernel* kernel){
  this->kernel = kernel;
}
  
void BaseSourcePlane::clear_H(std::string reg_scheme){
  if( reg_scheme == "all" ){
    this->H.clear();
  } else {
    if( !this->find_reg(reg_scheme) ){
      this->H.erase(reg_scheme);
    }
  }
}

bool BaseSourcePlane::find_reg(std::string reg_scheme){
  if( this->H.find(reg_scheme) == this->H.end() ){
    return false;
  } else {
    return true;
  }
}

void BaseSourcePlane::print_reg(){
  for(std::map<std::string,mytable>::iterator it=this->H.begin();it!=this->H.end();it++){
    std::cout << it->first << std::endl;
  }
}
