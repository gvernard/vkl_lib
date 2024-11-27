#include <cmath>
#include <iostream>
#include <cstring>

#include "json/json.h"

#include "vkllib.hpp"


int main(int argc, char* argv[]){
  Json::Value mass_models;

  double theta_E = 1.57;
  double x0 = 0.0;
  double y0 = 0.0;
  double pa = 23.4;
  double q = 0.999;
  double gam = 2.0;

  

  // Define a basic SIE mass model
  Json::Value pars;
  pars["theta_E"] = theta_E;
  pars["q"]  = q;
  pars["pa"] = pa;
  pars["x0"] = x0;
  pars["y0"] = y0;
  
  Json::Value mass_model;
  mass_model["type"] = "sie";
  mass_model["pars"] = pars;
  mass_models.append(mass_model);  



  // Define a Multipole mass model
  Json::Value pars2;
  pars2["theta_E"] = theta_E;
  pars2["x0"] = x0;
  pars2["y0"] = y0;  
  pars2["gam"]  = gam;
  pars2["phi_m"] = 0;
  pars2["m"] = 4;
  pars2["k_m"] = 1.0;

  Json::Value mass_model2;
  mass_model2["type"] = "multipole";
  mass_model2["pars"] = pars2;
  mass_models.append(mass_model2);  


  
  // Define a SPEMD mass model
  Json::Value pars3;
  pars3["theta_E"] = theta_E;
  pars3["q"]  = q;
  pars3["pa"] = pa;
  pars3["x0"] = x0;
  pars3["y0"] = y0;
  pars3["gam"]  = gam;
  pars3["s"] = 0.0001;

  Json::Value mass_model3;
  mass_model3["type"] = "spemd";
  mass_model3["pars"] = pars3;
  mass_models.append(mass_model3);




  



  // Create the mass model collection
  vkl::CollectionMassModels mmodel = vkl::JsonParsers::parse_mass_model(mass_models);


  


  // create a RectGrid object
  int Nx = 200;
  int Ny = 200;
  double b1_min = -2.0;
  double b1_max = 2.0;
  double b2_min = -2.0;
  double b2_max = 2.0;
  //  RectGrid mygrid(Nx,Ny,b1_min,b1_max,b2_min,b2_max);
  std::map<std::string,std::string> options{{"dev1_accu","1"},{"dev2_accu","2"},{"interp","bilinear"}};
  vkl::RectGrid mygrid_A(Nx,Ny,b1_min,b1_max,b2_min,b2_max,options);
  vkl::RectGrid mygrid_B(Nx,Ny,b1_min,b1_max,b2_min,b2_max,options);
  vkl::RectGrid mygrid_C(Nx,Ny,b1_min,b1_max,b2_min,b2_max,options);

  // calculate some function z at its pixels
  double xdefl,ydefl;
  for(int i=0;i<mygrid_A.Ny;i++){
    double y = mygrid_A.center_y[i];  
    for(int j=0;j<mygrid_A.Nx;j++){
      double x = mygrid_A.center_x[j];

      mmodel.models[0]->defl(x,y,xdefl,ydefl);
      mygrid_A.z[i*mygrid_A.Nx+j] = xdefl;

      mmodel.models[1]->defl(x,y,xdefl,ydefl);
      mygrid_B.z[i*mygrid_A.Nx+j] = xdefl;

      mmodel.models[2]->defl(x,y,xdefl,ydefl);
      mygrid_C.z[i*mygrid_A.Nx+j] = xdefl;
    }
  }
  vkl::FitsInterface::writeFits(Nx,Ny,mygrid_A.z,"sie_x.fits");
  vkl::FitsInterface::writeFits(Nx,Ny,mygrid_B.z,"multi_x.fits");
  vkl::FitsInterface::writeFits(Nx,Ny,mygrid_C.z,"spemd_x.fits");





  


  
  /*
  double xin = 5.9398;
  double yin = 6.14047;
  double xdefl,ydefl;
  for(int i=0;i<mmodel.models.size();i++){
    std::cout << mmodel.models[i]->mass_type << std::endl;
    mmodel.models[i]->defl(xin,yin,xdefl,ydefl);
    mmodel.models[i]->printMassPars();
    printf("(%8.6f,%8.6f)  <->  (%8.6f,%8.6f)\n",xin,yin,xin-xdefl,yin-ydefl);    
  }
  */

  
  return 0;
}
