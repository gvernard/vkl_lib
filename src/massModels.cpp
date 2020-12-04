#include "massModels.hpp"

//Abstract class: BaseMassModel
//===============================================================================================================
void BaseMassModel::setMassPars(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    if( nlpars[i]->nam == "pa" ){
      this->mpars[nlpars[i]->nam] = (nlpars[i]->val + 90.0) * 0.01745329251; // in rad
      //this->mpars[nlpars[i]->nam] = (nlpars[i]->val) * 0.01745329251; // in rad
    } else {
      this->mpars[nlpars[i]->nam] = nlpars[i]->val;
    }
  }
}

void BaseMassModel::printMassPars(){
  typedef std::map<std::string,double>::iterator some_type;
  for(some_type iterator=this->mpars.begin();iterator!=this->mpars.end();iterator++){
    std::cout << iterator->first << " " << iterator->second << std::endl;
    //    printf("%8s: %7.3f\n",iterator->second->nam,iterator->second->val);
  }
}


//Derived class from BaseMassModel: Sie (Singular Isothermal Ellipsoid)
//===============================================================================================================
Sie::Sie(std::vector<Nlpar*> nlpars){
  this->n = 5;
  this->mass_type = "sie";
  setMassPars(nlpars);
}

void Sie::defl(double xin,double yin,double& xout,double& yout){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"];
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);
  
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
  
  double fac   = sqrt(1.0-q*q);
  double omega = sqrt(q*q*x_t*x_t + y_t*y_t); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.

  double ax_t = (b/fac)*atan(x_t*fac/omega);
  double ay_t = (b/fac)*atanh(y_t*fac/omega);
  
  //rotate back according to position angle, no need to translate (this is equivalent to rotating by -pa using the same equations as above)
  double ax =  ax_t*cos(pa) - ay_t*sin(pa);
  double ay =  ax_t*sin(pa) + ay_t*cos(pa);
  
  xout = ax;
  yout = ay;
}

double Sie::kappa(double xin,double yin){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"];
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

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

  double omega = q*q*x_t*x_t + y_t*y_t; // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double k = 0.5*b/sqrt(omega);
  return k;
}


void Sie::gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  double b  = this->mpars["b"];
  double q  = this->mpars["q"];
  double pa = this->mpars["pa"];
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

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

  double omega     = q*q*x_t*x_t + y_t*y_t; // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double z         = atan2(y_t,x_t);        // this is because cos2z = (y^2-x^2)/(x^2+y^2), and sin2z = 2xy/(x^2+y^2)
  double gamma_x_t = -0.5*b*cos(2*z)/sqrt(omega);
  double gamma_y_t = -0.5*b*sin(2*z)/sqrt(omega);
  
  //rotate back according to position angle
  gamma_mag = hypot(gamma_x_t,gamma_y_t);
  gamma_phi = 0.5*atan2(gamma_y_t,gamma_x_t) + pa; // in rad
}

double Sie::psi(double xin,double yin){
  double ax,ay;
  this->defl(xin,yin,ax,ay);
  double psi = xin*ax + yin*ay;
  return psi;
}

//Derived class from BaseMassModel: Spemd (Softened Power-law Elliptical Mass Density)
//===============================================================================================================
Spemd::Spemd(std::vector<Nlpar*> nlpars){
  this->n = 7;
  this->mass_type = "spemd";
  setMassPars(nlpars);
}

void Spemd::defl(double xin,double yin,double& xout,double& yout){
  double q  = this->mpars["q"];
  double e  = this->mpars["e"];
  double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e)/(q*2.0); // this is correct!
  //double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e);
  //double b = 2*(1-e)*pow(this->mpars["b"],2.0)/((1+q)*(pow(this->mpars["b"],2-2*e)-1));
  double pa = this->mpars["pa"];
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];
  double s2 = this->mpars["s"] * this->mpars["s"];
  double defl[2] = {0.0,0.0};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  fastelldefl_(&x_t,&y_t,&b,&e,&q,&s2,defl);

  double ax_t = defl[0];
  double ay_t = defl[1];

  //rotate back according to position angle, no need to translate
  double ax =  ax_t*cos(pa) - ay_t*sin(pa);
  double ay =  ax_t*sin(pa) + ay_t*cos(pa);
  
  xout = ax;
  yout = ay;
}

double Spemd::kappa(double xin,double yin){
  double q  = this->mpars["q"];
  double e  = this->mpars["e"];
  double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e)/(q*2.0); // this is correct!
  //double b = 2*(1-e)*pow(this->mpars["b"],2.0)/((1+q)*(pow(this->mpars["b"],2-2*e)-1));
  double pa = this->mpars["pa"];
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];
  double s2 = this->mpars["s"] * this->mpars["s"];

  //rotate the coordinate system according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  double omega = x_t*x_t + y_t*y_t/(q*q); // this does not depend on using omega, omega', or zeta, as long as I change correctly to these elliptical radii.
  double k = b/pow(omega+s2,e);
  return k;
}

void Spemd::gamma(double xin,double yin,double& gamma_mag,double& gamma_phi){
  double q  = this->mpars["q"];
  double e  = this->mpars["e"];
  double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e)/(q*2.0); // this is correct!
  //double b  = pow(this->mpars["b"],2.0*e)*(2.0-2.0*e);
  //double b = 2*(1-e)*pow(this->mpars["b"],2.0)/((1+q)*(pow(this->mpars["b"],2-2*e)-1));
  double pa = this->mpars["pa"];
  double x0 = this->mpars["x0"];
  double y0 = this->mpars["y0"];
  double s2 = this->mpars["s"] * this->mpars["s"];
  double defl[2] = {0.0,0.0};
  double jacob[2*2] = {0.0,0.0,0.0,0.0};

  //rotate according to position angle and translate to the lens center
  double x_t =  (xin-x0)*cos(pa) + (yin-y0)*sin(pa);
  double y_t = -(xin-x0)*sin(pa) + (yin-y0)*cos(pa);

  fastellmag_(&x_t,&y_t,&b,&e,&q,&s2,defl,jacob);

  double gamma_x_t = (jacob[0] - jacob[3])/2.0;
  double gamma_y_t = jacob[1];
  
  //rotate back according to position angle
  gamma_mag = hypot(gamma_x_t,gamma_y_t);
  gamma_phi = 0.5*atan2(gamma_y_t,gamma_x_t) + pa; // in rad
}

//Derived class from BaseMassModel: Pert (perturbations on a grid)
//===============================================================================================================
void Pert::replaceDpsi(double* new_dpsi){
  for(int k=0;k<this->Sm;k++){
    this->z[k] = new_dpsi[k];
  }
  updatePert();
}

void Pert::addDpsi(double* corrections){
  for(int k=0;k<this->Sm;k++){
    this->z[k] += corrections[k];
  }
  updatePert();
}

void Pert::updatePert(){
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





//Class: CollectionMassModels
//===============================================================================================================
CollectionMassModels::CollectionMassModels(){
  this->mpars["gx"] = 0.0;
  this->mpars["gy"] = 0.0;
}
CollectionMassModels::CollectionMassModels(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    if( nlpars[i]->nam == "phi" ){
      nlpars[i]->val = (nlpars[i]->val + 90.0)*0.01745329251; // convert from EoN to normal cartesian
    }
  }
  this->setPhysicalPars(nlpars);
};
CollectionMassModels::~CollectionMassModels(){
  for(int i=0;i<this->models.size();i++){
    delete(this->models[i]);
  }
  mpars.clear();
};

void CollectionMassModels::printPhysPars(){
  typedef std::map<std::string,double>::iterator some_type;
  for(some_type iterator=this->mpars.begin();iterator!=this->mpars.end();iterator++){
    std::cout << iterator->first << " " << iterator->second << std::endl;
    //    printf("%8s: %7.3f\n",iterator->second->nam,iterator->second->val);
  }
}

void CollectionMassModels::setPhysicalPars(std::vector<Nlpar*> nlpars){
  for(int i=0;i<nlpars.size();i++){
    this->mpars[nlpars[i]->nam] = nlpars[i]->val;
  }
  this->mpars["gx"] = this->mpars["g"]*cos(2*this->mpars["phi"]);
  this->mpars["gy"] = this->mpars["g"]*sin(2*this->mpars["phi"]);
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
  xout = (1.0-this->mpars["gx"])*xin - this->mpars["gy"]*yin - ax;
  yout = (1.0+this->mpars["gx"])*yin - this->mpars["gy"]*xin - ay;
}

void CollectionMassModels::all_defl(ImagePlane* image){
  double dumx = 0.0;
  double dumy = 0.0;
  double xin,yin;
  
  for(int j=0;j<image->Nm;j++){
    xin = image->grid->center_x[j];
    yin = image->grid->center_y[j];
    double ax   = 0.0;
    double ay   = 0.0;
    for(int i=0;i<this->models.size();i++){
      this->models[i]->defl(xin,yin,dumx,dumy);
      ax += dumx;
      ay += dumy;
    }
    image->defl_x[j] = (1.0-this->mpars["gx"])*xin - this->mpars["gy"]*yin - ax;
    image->defl_y[j] = (1.0+this->mpars["gx"])*yin - this->mpars["gy"]*xin - ay;
  }
}

double CollectionMassModels::all_kappa(double xin,double yin){
  double k = 0.0;
  for(int i=0;i<this->models.size();i++){
    k += this->models[i]->kappa(xin,yin);
  }
  return k;
}

void CollectionMassModels::all_kappa(ImagePlane* image,ImagePlane* kappa_tot){
  double xin,yin;
  for(int j=0;j<image->Nm;j++){
    xin = image->grid->center_x[j];
    yin = image->grid->center_y[j];
    double k = 0.0;
    for(int i=0;i<this->models.size();i++){
      k += this->models[i]->kappa(xin,yin);
    }
    kappa_tot->grid->z[j] = k;
  }
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
  gamma_x += this->mpars["gx"];
  gamma_y += this->mpars["gy"];
  gamma_mag = hypot(gamma_x,gamma_y);
  gamma_phi = 0.5*atan2(gamma_y,gamma_x);
}

void CollectionMassModels::all_gamma(ImagePlane* image,ImagePlane* gamma_mag,ImagePlane* gamma_phi){
  double xin,yin;
  double mag,phi;
  for(int j=0;j<image->Nm;j++){
    xin = image->grid->center_x[j];
    yin = image->grid->center_y[j];
    double gamma_x = 0.0;
    double gamma_y = 0.0;
    for(int i=0;i<this->models.size();i++){
      this->models[i]->gamma(xin,yin,mag,phi);
      gamma_x += mag*cos(2*phi);
      gamma_y += mag*sin(2*phi);
    }
    gamma_x += this->mpars["gx"];
    gamma_y += this->mpars["gy"];

    gamma_mag->grid->z[j] = hypot(gamma_x,gamma_y);
    gamma_phi->grid->z[j] = 0.5*atan2(gamma_y,gamma_x);
  }
}

double CollectionMassModels::detJacobian(double xin,double yin){
  double gamma_mag,gamma_phi;
  this->all_gamma(xin,yin,gamma_mag,gamma_phi);
  return pow(1.0 - this->all_kappa(xin,yin),2) - pow(gamma_mag,2);
}

void CollectionMassModels::detJacobian(ImagePlane* image,ImagePlane* detA){
  double xin,yin;
  for(int j=0;j<image->Nm;j++){
    xin = image->grid->center_x[j];
    yin = image->grid->center_y[j];
    detA->grid->z[j] = this->detJacobian(xin,yin);
  }
}

double CollectionMassModels::all_psi(double xin,double yin){
  double psi = 0.0;
  for(int i=0;i<this->models.size();i++){
    psi += this->models[i]->psi(xin,yin);
  }
  psi += 0.5*this->mpars["gx"]*(xin*xin-yin*yin) + this->mpars["gy"]*xin*yin;
  return psi;
}

void CollectionMassModels::all_psi(ImagePlane* image,ImagePlane* psi_tot){
  double xin,yin;
  for(int j=0;j<image->Nm;j++){
    xin = image->grid->center_x[j];
    yin = image->grid->center_y[j];
    double psi = 0.0;
    for(int i=0;i<this->models.size();i++){
      psi += this->models[i]->psi(xin,yin);
    }
    psi += this->mpars["gx"]*(xin*xin-yin*yin)/2.0 + this->mpars["gy"]*xin*yin;
    psi_tot->grid->z[j] = psi;
  }
}
