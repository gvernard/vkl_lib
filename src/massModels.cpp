#include "massModels.hpp"

#include <algorithm>

extern "C"{
  void fastelldefl_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* defl);
  void fastellmag_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* defl,double* jacob);
  void ellipphi_(double* x1,double* x2,double* b,double* g,double* q,double* s2,double* phi);
}



//Abstract class: BaseMassModel
//===============================================================================================================
void BaseMassModel::printMassPars(){
  for(std::map<std::string,double>::iterator it=mpars.begin();it!=mpars.end();it++){
    std::cout << it->first << " " << it->second << std::endl;
    //printf("%8s: %7.3f\n",iterator->second->nam,iterator->second->val);
  }
}

//Class: CollectionMassModels
//===============================================================================================================
CollectionMassModels::CollectionMassModels(const CollectionMassModels& other){
  this->models.resize(other.models.size());
  for(int i=0;i<other.models.size();i++){
    this->models[i] = other.models[i];
  }
}

CollectionMassModels::~CollectionMassModels(){
  for(int i=0;i<this->models.size();i++){
    delete(models[i]);
  }
}

void CollectionMassModels::all_defl(double xin,double yin,double& xout,double& yout){
  double ax   = 0.0;
  double ay   = 0.0;
  double dumx = 0.0;
  double dumy = 0.0;
  for(int i=0;i<this->models.size();i++){
    this->models[i]->defl(xin,yin,dumx,dumy);
    ax += dumx;
    ay += dumy;
  }
  xout = xin - ax;
  yout = yin - ay;
}

double CollectionMassModels::all_kappa(double xin,double yin){
  double kappa = 0.0;
  for(int i=0;i<this->models.size();i++){
    kappa += this->models[i]->kappa(xin,yin);
  }
  return kappa;
}

void CollectionMassModels::all_gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  double gamma_x = 0.0;
  double gamma_y = 0.0;
  double mag = 0.0;
  double phi = 0.0;
  for(int i=0;i<this->models.size();i++){
    this->models[i]->gamma(xin,yin,mag,phi);
    gamma_x += mag*cos(2*phi);
    gamma_y += mag*sin(2*phi);
  }
  gamma_mag = hypot(gamma_x,gamma_y);
  gamma_phi = 0.5*atan2(gamma_y,gamma_x);
}

double CollectionMassModels::detJacobian(double xin,double yin){
  double gamma_mag,gamma_phi;
  this->all_gamma(xin,yin,gamma_mag,gamma_phi);
  return pow(1.0 - this->all_kappa(xin,yin),2) - pow(gamma_mag,2);
}

double CollectionMassModels::all_psi(double xin,double yin){
  double psi = 0.0;
  for(int i=0;i<this->models.size();i++){
    psi += this->models[i]->psi(xin,yin);
  }
  return psi;
}

void CollectionMassModels::getExtent(double& xmin,double& xmax,double& ymin,double& ymax){
  std::vector<double> vx;
  std::vector<double> vy;
  double tmp_xmin,tmp_xmax,tmp_ymin,tmp_ymax;
  for(int i=0;i<this->models.size();i++){
    if( this->models[i]->mass_type == "pert" ){
      Pert* pert = static_cast<Pert*> (this->models[i]);
      vx.push_back(pert->xmin);
      vx.push_back(pert->xmax);
      vy.push_back(pert->ymin);
      vy.push_back(pert->ymax);
    } else if( this->models[i]->mass_type == "external_shear" ){
      vx.push_back(-2); // just assuming a range of 2x2 arcsec
      vx.push_back(2);
      vy.push_back(-2);
      vy.push_back(2);
    } else {
      vx.push_back(this->models[i]->mpars["x0"] - 2.5*this->models[i]->mpars["theta_E"]);
      vx.push_back(this->models[i]->mpars["x0"] + 2.5*this->models[i]->mpars["theta_E"]);
      vy.push_back(this->models[i]->mpars["y0"] - 2.5*this->models[i]->mpars["theta_E"]);
      vy.push_back(this->models[i]->mpars["y0"] + 2.5*this->models[i]->mpars["theta_E"]);
    }
  }

  auto resultx = std::minmax_element(std::begin(vx),std::end(vx));
  xmin = *resultx.first;
  xmax = *resultx.second;
  auto resulty = std::minmax_element(std::begin(vy),std::end(vy));
  ymin = *resulty.first;
  ymax = *resulty.second;
}


//Derived class from BaseMassModel: ExternalShear
//===============================================================================================================
ExternalShear::ExternalShear(std::map<std::string,double> pars): BaseMassModel(2,"external_shear"){
  updateMassPars(pars);
}

void ExternalShear::updateMassPars(std::map<std::string,double> pars){
  // First modify all parameter values that come from the outside
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    mpars[it->first] = it->second;
  }
  // Then convert any of the following accordingly, if they exist in the incoming pars
  if( pars.find("phi") != pars.end() ){
    mpars["phi"] = (pars["phi"] + 90.0) * 0.01745329251; // in rad
  }
  mpars["gx"] = mpars["g"]*cos(2*mpars["phi"]);
  mpars["gy"] = mpars["g"]*sin(2*mpars["phi"]);
}

void ExternalShear::defl(double xin,double yin,double& xout,double& yout){
  xout = mpars["gx"]*xin + mpars["gy"]*yin;
  yout = mpars["gy"]*xin - mpars["gx"]*yin;
}

double ExternalShear::kappa(double xin,double yin){
  return 0.0;
}

void ExternalShear::gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  gamma_mag = mpars["g"];
  gamma_phi = mpars["phi"];
}

double ExternalShear::psi(double xin,double yin){
  return 0.5*mpars["gx"]*(xin*xin-yin*yin) + mpars["gy"]*xin*yin;
}

//Derived class from BaseMassModel: Sie (Singular Isothermal Ellipsoid)
//===============================================================================================================
Sie::Sie(std::map<std::string,double> pars): BaseMassModel(5,"sie"){
  updateMassPars(pars);
}

void Sie::updateMassPars(std::map<std::string,double> pars){
  // First modify all parameter values that come from the outside
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    mpars[it->first] = it->second;
  }
  // Then convert any of the following accordingly, if they exist in the incoming pars
  if( pars.find("pa") != pars.end() ){
    mpars["pa"] = (mpars["pa"] + 90.0) * 0.01745329251; // in rad
  }
  if( pars.find("theta_E") != pars.end() || pars.find("q") != pars.end() ){
    mpars["b"] = mpars["theta_E"]/sqrt(mpars["q"]);
  }
}

void Sie::defl(double xin,double yin,double& xout,double& yout){
  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);
  
  this->check_close_to_origin(x_t,y_t);
  
  double fac   = sqrt(1.0-mpars["q"]*mpars["q"]);
  double omega = sqrt(mpars["q"]*mpars["q"]*x_t*x_t + y_t*y_t); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.

  double ax_t = (mpars["b"]*mpars["q"]/fac)*atan(x_t*fac/omega);
  double ay_t = (mpars["b"]*mpars["q"]/fac)*atanh(y_t*fac/omega);
  
  //rotate back according to position angle, no need to translate (this is equivalent to rotating by -pa using the same equations as above)
  xout = ax_t*cos(mpars["pa"]) - ay_t*sin(mpars["pa"]);
  yout = ax_t*sin(mpars["pa"]) + ay_t*cos(mpars["pa"]);
}

double Sie::kappa(double xin,double yin){
  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);

  this->check_close_to_origin(x_t,y_t);

  double omega = x_t*x_t + y_t*y_t/(mpars["q"]*mpars["q"]); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  return 0.5*mpars["b"]/sqrt(omega);
}

void Sie::gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);

  this->check_close_to_origin(x_t,y_t);
  
  double omega     = x_t*x_t + y_t*y_t/(mpars["q"]*mpars["q"]); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double z         = atan2(y_t,x_t);        // this is because cos2z = (y^2-x^2)/(x^2+y^2), and sin2z = 2xy/(x^2+y^2)
  double gamma_x_t = -0.5*mpars["b"]*cos(2*z)/sqrt(omega);
  double gamma_y_t = -0.5*mpars["b"]*sin(2*z)/sqrt(omega);
  
  //rotate back according to position angle
  gamma_mag = hypot(gamma_x_t,gamma_y_t);
  gamma_phi = 0.5*atan2(gamma_y_t,gamma_x_t) + mpars["pa"]; // in rad
}

double Sie::psi(double xin,double yin){
  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);
  
  this->check_close_to_origin(x_t,y_t);
  
  double fac   = sqrt(1.0-mpars["q"]*mpars["q"]);
  double omega = sqrt(mpars["q"]*mpars["q"]*x_t*x_t + y_t*y_t); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.

  double ax_t = (mpars["b"]*mpars["q"]/fac)*atan(x_t*fac/omega);
  double ay_t = (mpars["b"]*mpars["q"]/fac)*atanh(y_t*fac/omega);
  
  double psi = x_t*ax_t + y_t*ay_t;
  return psi;
}

void Sie::check_close_to_origin(double& x_t,double& y_t){
  if( fabs(x_t) < 0.0001 && fabs(y_t) < 0.0001 ){
    if( std::signbit(x_t) ){
      x_t = -0.0001;
    } else {
      x_t =  0.0001;
    }
    if( std::signbit(y_t) ){
      y_t = -0.0001;
    } else {
      y_t =  0.0001;
    }
  }
}

//Derived class from BaseMassModel: Spemd (Softened Power-law Elliptical Mass Density)
//===============================================================================================================
Spemd::Spemd(std::map<std::string,double> pars): BaseMassModel(7,"spemd"){
  updateMassPars(pars);
}

void Spemd::updateMassPars(std::map<std::string,double> pars){
  // First modify all parameter values that come from the outside
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    mpars[it->first] = it->second;
  }
  // Then convert any of the following accordingly, if they exist in the incoming pars
  if( pars.find("pa") != pars.end() ){
    mpars["pa"] = (mpars["pa"] + 90.0) * 0.01745329251; // in rad
  }
  if( pars.find("gam") != pars.end() ){
    mpars["e"] = (mpars["gam"]-1.0)/2.0;
  }
  if( pars.find("theta_E") != pars.end() || pars.find("q") != pars.end() ){
    mpars["b"] = pow(mpars["theta_E"]/sqrt(mpars["q"]),2*mpars["e"])*(1.0-mpars["e"]);
    //mpars["b"] = pow(pars["b"],2.0*mpars["e"])*(2.0-2.0*mpars["e"])/(2.0*mpars["q"]); // this is the old implementation
    //mpars["b"] = pow(pars["b"],2.0*mpars["e"])*(2.0-2.0*mpars["e"])/(2.0*pow(mpars["q"],mpars["e"])); // this is for GLEE
  }
  if( pars.find("s") != pars.end() ){
    mpars["s2"] = pow(pars["s"],2);
  }
}

void Spemd::defl(double xin,double yin,double& xout,double& yout){
  double defl[2] = {0.0,0.0};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);

  fastelldefl_(&x_t,&y_t,&mpars["b"],&mpars["e"],&mpars["q"],&mpars["s2"],defl);

  double ax_t = defl[0];
  double ay_t = defl[1];

  //rotate back according to position angle, no need to translate
  xout = ax_t*cos(mpars["pa"]) - ay_t*sin(mpars["pa"]);
  yout = ax_t*sin(mpars["pa"]) + ay_t*cos(mpars["pa"]);
}

double Spemd::kappa(double xin,double yin){
  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);

  double omega2 = x_t*x_t + y_t*y_t/(mpars["q"]*mpars["q"]); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  return mpars["b"]/pow(omega2+mpars["s2"],mpars["e"]);
  //return (2.0*pow(mpars["q"],mpars["e"])/(1+mpars["q"]))*mpars["b"]/pow(omega+mpars["s2"],mpars["e"]);
}

void Spemd::gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  double defl[2] = {0.0,0.0};
  double jacob[2*2] = {0.0,0.0,0.0,0.0};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);

  fastellmag_(&x_t,&y_t,&mpars["b"],&mpars["e"],&mpars["q"],&mpars["s2"],defl,jacob);

  double gamma_x_t = (jacob[0] - jacob[3])/2.0;
  double gamma_y_t = jacob[1];
  
  //rotate back according to position angle
  gamma_mag = hypot(gamma_x_t,gamma_y_t);
  gamma_phi = 0.5*atan2(gamma_y_t,gamma_x_t) + mpars["pa"]; // in rad
}

double Spemd::psi(double xin,double yin){
  double psi;
  
  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-mpars["x0"])*cos(mpars["pa"]) + (yin-mpars["y0"])*sin(mpars["pa"]);
  double y_t = -(xin-mpars["x0"])*sin(mpars["pa"]) + (yin-mpars["y0"])*cos(mpars["pa"]);

  ellipphi_(&x_t,&y_t,&mpars["b"],&mpars["e"],&mpars["q"],&mpars["s2"],&psi);
  
  return psi;
}

//Derived class from BaseMassModel: Pert (perturbations on a grid)
//===============================================================================================================
Pert::Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax): BaseMassModel(Nx*Ny,"pert"),FixedSource(Nx,Ny,xmin,xmax,ymin,ymax){};
Pert::Pert(int Nx,int Ny,ImagePlane* image): BaseMassModel(Nx*Ny,"pert"),FixedSource(Nx,Ny,image->grid->xmin,image->grid->xmax,image->grid->ymin,image->grid->ymax){};
Pert::Pert(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath): BaseMassModel(Nx*Ny,"pert"),FixedSource(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
  this->updateDerivatives();
}

void Pert::updateMassPars(std::string way,double* new_dpsi){
  if( way == "replace" ){
    for(int k=0;k<this->Sm;k++){
      this->z[k] = new_dpsi[k];
    }
  } else if( way == "add" ){
    for(int k=0;k<this->Sm;k++){
      this->z[k] += new_dpsi[k];
    }
  }
  this->updateDerivatives();
}

void Pert::updateDerivatives(){
  this->calculate_zx();
  this->calculate_zy();
  this->calculate_zxx();
  this->calculate_zyy();
  this->calculate_zxy();
}

void Pert::defl(double xin,double yin,double& xout,double& yout){
  xout = (this->*interp2d)(xin,yin,this->zx);
  yout = (this->*interp2d)(xin,yin,this->zy);
}

double Pert::psi(double xin,double yin){
  return (this->*interp2d)(xin,yin,this->z);
}

double Pert::kappa(double xin,double yin){
  double zxx = (this->*interp2d)(xin,yin,this->zxx);
  double zyy = (this->*interp2d)(xin,yin,this->zyy);
  return (zxx+zyy)/2.0;
}

void Pert::gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  double zxx = (this->*interp2d)(xin,yin,this->zxx);
  double zyy = (this->*interp2d)(xin,yin,this->zyy);
  double gx  = (zxx-zyy)/2.0;
  double gy  = (this->*interp2d)(xin,yin,this->zxy);
  gamma_mag  = hypot(gx,gy);
  gamma_phi  = 0.5*atan2(gy,gx); // in rad
}


//Derived class from BaseMassModel: ADD YOUR DESIRED MASS MODEL IMPLEMENTATION HERE
//===============================================================================================================
