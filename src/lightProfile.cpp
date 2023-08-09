#define _USE_MATH_DEFINES

#include "lightProfile.hpp"
#include "fitsInterface.hpp"
#include "rectGrid.hpp"

#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm> // for minmax_element
#include <iostream>
#include <stdexcept>
#include <assert.h>

/*
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2_algorithms.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int,K>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2<unsigned int,K>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                          Delaunay;
typedef K::Point_2                                                     Point;
*/


using namespace vkl;

//BaseClass: BaseProfile
//===============================================================================================================
BaseProfile::BaseProfile(const BaseProfile& other){
  this->Npars = other.Npars;
  this->profile_type = other.profile_type;
  this->upsilon = other.upsilon;
  this->M_tot = other.M_tot;
  this->ZP = other.ZP;
}



//Class: CollectionProfiles
//===============================================================================================================
CollectionProfiles::CollectionProfiles(const CollectionProfiles& other){
  this->profiles.resize(other.profiles.size());
  for(int i=0;i<other.profiles.size();i++){
    this->profiles[i] = other.profiles[i];
  }
}

CollectionProfiles::~CollectionProfiles(){
  for(int i=0;i<this->profiles.size();i++){
    delete(this->profiles[i]);
  }
}

double CollectionProfiles::all_values(double xin,double yin){
  double value = 0.0;
  for(int i=0;i<this->profiles.size();i++){
    value += this->profiles[i]->value(xin,yin);
  }
  return value;
}

double CollectionProfiles::all_values_to_mass(double xin,double yin){
  double value = 0.0;
  for(int i=0;i<this->profiles.size();i++){
    value += this->profiles[i]->value_to_mass(xin,yin);
  }
  return value;
}

void CollectionProfiles::getExtent(double& xmin,double& xmax,double& ymin,double& ymax){
  std::vector<double> vx;
  std::vector<double> vy;
  for(int i=0;i<this->profiles.size();i++){
    this->profiles[i]->get_extent(xmin,xmax,ymin,ymax);
    vx.push_back(xmin);
    vx.push_back(xmax);
    vy.push_back(ymin);
    vy.push_back(ymax);
  }

  auto resultx = std::minmax_element(std::begin(vx),std::end(vx));
  xmin = *resultx.first;
  xmax = *resultx.second;
  auto resulty = std::minmax_element(std::begin(vy),std::end(vy));
  ymin = *resulty.first;
  ymax = *resulty.second;
}

void CollectionProfiles::write_all_profiles(const std::string filepath){
  double xmin,xmax,ymin,ymax;
  this->getExtent(xmin,xmax,ymin,ymax);
  double w = xmax - xmin;
  double h = ymax - ymin;
  if( h>w ){
    w = h;
  }
  
  RectGrid grid(300,300,xmin,xmin+w,ymin,ymin+w);
  for(int i=0;i<grid.Ny;i++){
    for(int j=0;j<grid.Nx;j++){
      grid.z[i*grid.Nx+j] = this->all_values(grid.center_x[j],grid.center_y[i]);
    }
  }
  std::vector<std::string> keys{"xmin","xmax","ymin","ymax"};
  std::vector<std::string> values{std::to_string(grid.xmin),std::to_string(grid.xmax),std::to_string(grid.ymin),std::to_string(grid.ymax)};
  std::vector<std::string> descriptions{"left limit of the frame","right limit of the frame","bottom limit of the frame","top limit of the frame"};
  FitsInterface::writeFits(grid.Nx,grid.Ny,grid.z,keys,values,descriptions,filepath);
}



// Class Sersic
//===============================================================================================================
Sersic::Sersic(std::map<std::string,double> pars): BaseProfile(9,"sersic"){
  std::map<std::string,double>::iterator end = pars.end();

  // Check that either M_tot or i_eff is present
  std::map<std::string,double>::iterator it1,it2;
  it1 = pars.find("M_tot");
  it2 = pars.find("i_eff");
  assert( !(it1 == end && it2 == end) && "Error: At least one of 'M_tot' or 'i_eff' must be present in Sersic profile constructor!" );
  assert( !(it1 != end && it2 != end) && "Error: Both 'M_tot' and 'i_eff' cannot be present in Sersic profile constructor!" );

  // Check that all mandatory parameters by name are present
  std::vector<std::string> missing;  
  std::vector<std::string> mandatory{"ZP","r_eff","x0","y0","n","pa","q"};
  for(int i=0;i<mandatory.size();i++){
    if( pars.find(mandatory[i]) == end ){
      missing.push_back(mandatory[i]);
    }
  }
  if( missing.size() != 0 ){
    std::stringstream message("Error: The following parameters are missing in Sersic profile constructor: ");
    for(int i=0;i<missing.size();i++){
      message << " " << missing[i];
    }
    std::string tmp = message.str();
    throw std::runtime_error(tmp.c_str());
  }
  
  // ZP will never be updated so it is not touched by updateProfilePars
  this->ZP = pars["ZP"];
  
  updateProfilePars(pars);
}

Sersic::Sersic(const Sersic& other): BaseProfile(other){
  this->Reff = other.Reff;
  this->Ieff = other.Ieff;
  this->x0   = other.x0;
  this->y0   = other.y0;
  this->n    = other.n;
  this->pa   = other.pa;
  this->q    = other.q;
  this->upsilon_exp = other.upsilon_exp;

  this->p_xmin = other.p_xmin;
  this->p_xmax = other.p_xmax;
  this->p_ymin = other.p_ymin;
  this->p_ymax = other.p_ymax;
  this->bn = other.bn;
  this->cospa = other.cospa;
  this->sinpa = other.sinpa;
}

void Sersic::updateProfilePars(std::map<std::string,double> pars){
  std::map<std::string,double>::iterator end = pars.end();

  // Update BaseProfile parameter first
  if( pars.find("upsilon") != end ){
    this->upsilon = pars["upsilon"]; // upsilon is a base class variable
  }

  // Single parameter updates
  if( pars.find("x0") != end ){
    this->x0 = pars["x0"];
  }
  if( pars.find("y0") != end ){
    this->y0 = pars["y0"];
  }
  if( pars.find("n") != end ){
    this->n = pars["n"];
    this->bn = 1.9992*this->n - 0.3271;//From Capaccioli 1989
  }
  if( pars.find("pa") != end ){
    this->pa = (pars["pa"] + 90.0) * 0.01745329251; // in rad
    this->cospa = cos(this->pa);
    this->sinpa = sin(this->pa);
  }
  if( pars.find("upsilon_exp") != end ){
    this->upsilon_exp = pars["upsilon_exp"];
  }

  // Update Reff and q that depend on each other
  std::map<std::string,double>::iterator it_q,it_r;
  it_q = pars.find("q");
  it_r = pars.find("r_eff"); // r_eff must be the semi-major axis
  if( it_q != end && it_r != end ){ // both 'q' and 'r_eff' are found
    this->q = pars["q"];
    this->Reff = this->q*pars["r_eff"];
  } else if( it_q != end && it_r == end ){ // 'q' is found but not 'r_eff'
    double old_q = this->q;
    double old_Reff = this->Reff;
    this->q = pars["q"];
    this->Reff = this->q*(old_Reff/old_q);
  } else if( it_q == end && it_r != end ){ // 'r_eff' is found but not 'q'
    this->Reff = this->q*pars["r_eff"];
  }

  // Update M_tot (I_eff) and then I_eff (M_tot).
  if( pars.find("M_tot") != end ){
    // See eq. 4 in Peng et al. 2010
    this->M_tot = pars["M_tot"];
    double fac = pow(this->Reff/this->q,2)*2*M_PI*this->n*exp(this->bn)*tgamma(2*this->n)*this->q/pow(this->bn,2*this->n);
    this->Ieff = pow(10.0,-0.4*(this->M_tot-this->ZP))/fac;
  } else if( pars.find("M_tot") != end ){
    this->Ieff = pars["i_eff"];
    double fac = pow(this->Reff/this->q,2)*2*M_PI*this->n*exp(this->bn)*tgamma(2*this->n)*this->q/pow(this->bn,2*this->n);
    this->M_tot = -2.5*log10(fac*this->Ieff) + this->ZP;
  }

  set_extent();
}

double Sersic::value(double x,double y){
  if( is_in_range(x,y) ){
    double u,v,r,fac2;
    u =  (x - this->x0)*this->cospa + (y - this->y0)*this->sinpa;
    v = -(x - this->x0)*this->sinpa + (y - this->y0)*this->cospa;
    r = hypot(this->q*u,v);
    fac2 = pow(r/this->Reff,1.0/this->n) - 1.0;
    return this->Ieff*exp(-this->bn*fac2);
  } else {
    return 0.0;
  }
}

double Sersic::value_to_mass(double x,double y){
  if( is_in_range(x,y) ){
    double u,v,r,fac,light;
    light = this->value(x,y);
    u =  (x - this->x0)*this->cospa + (y - this->y0)*this->sinpa;
    v = -(x - this->x0)*this->sinpa + (y - this->y0)*this->cospa;
    r = hypot(this->q*u,v);
    fac = this->upsilon*pow(r/this->Reff,this->upsilon_exp);
    return fac*light;
  } else {
    return 0.0;
  }
}

bool Sersic::is_in_range(double xin,double yin){
  if( xin < this->p_xmin || this->p_xmax < xin || yin < this->p_ymin || this->p_ymax < yin ){
    return false;
  } else {
    return true;
  }
}

void Sersic::get_extent(double& xmin,double& xmax,double& ymin,double& ymax){
  xmin = this->p_xmin;
  xmax = this->p_xmax;
  ymin = this->p_ymin;
  ymax = this->p_ymax;
}

void Sersic::set_extent(){
  double dx = fabs(3*this->Reff);
  this->p_xmin = this->x0 - dx;
  this->p_xmax = this->x0 + dx;
  double dy = fabs(3*this->Reff);
  this->p_ymin = this->y0 - dy;
  this->p_ymax = this->y0 + dy;
}

// Class Gauss
//===============================================================================================================
Gauss::Gauss(std::map<std::string,double> pars): BaseProfile(8,"gauss"){
  std::map<std::string,double>::iterator end = pars.end();

  // Check that either M_tot or i_eff is present
  std::map<std::string,double>::iterator it1,it2;
  it1 = pars.find("M_tot");
  it2 = pars.find("i_eff");
  assert( !(it1 == end && it2 == end) && "Error: At least one of 'M_tot' or 'i_eff' must be present in Gauss profile constructor!" );
  assert( !(it1 != end && it2 != end) && "Error: Both 'M_tot' and 'i_eff' cannot be present in Gauss profile constructor!" );

  // Check that all mandatory parameters by name are present
  std::vector<std::string> missing;  
  std::vector<std::string> mandatory{"ZP","r_eff","x0","y0","pa","q"};
  for(int i=0;i<mandatory.size();i++){
    if( pars.find(mandatory[i]) == end ){
      missing.push_back(mandatory[i]);
    }
  }
  if( missing.size() != 0 ){
    std::stringstream message("Error: The following parameters are missing in Gauss profile constructor: ");
    for(int i=0;i<missing.size();i++){
      message << " " << missing[i];
    }
    std::string tmp = message.str();
    throw std::runtime_error(tmp.c_str());
  }
  
  // ZP will never be updated so it is not touched by updateProfilePars
  this->ZP = pars["ZP"];

  updateProfilePars(pars);
}

Gauss::Gauss(const Gauss& other): BaseProfile(other){
  this->Reff = other.Reff;
  this->Ieff = other.Ieff;
  this->x0   = other.x0;
  this->y0   = other.y0;
  this->pa   = other.pa;
  this->q    = other.q;
  this->upsilon_exp = other.upsilon_exp;

  this->p_xmin = other.p_xmin;
  this->p_xmax = other.p_xmax;
  this->p_ymin = other.p_ymin;
  this->p_ymax = other.p_ymax;
  this->sdev_fac = other.sdev_fac;
  this->cospa = other.cospa;
  this->sinpa = other.sinpa;
}

void Gauss::updateProfilePars(std::map<std::string,double> pars){
  std::map<std::string,double>::iterator end = pars.end();

  // Update BaseProfile parameter first
  if( pars.find("ZP") != end ){
    this->ZP = pars["ZP"];
  }  
  if( pars.find("upsilon") != end ){
    this->upsilon = pars["upsilon"]; // upsilon is a base class variable
  }

  // Single parameter updates
  if( pars.find("x0") != end ){
    this->x0 = pars["x0"];
  }
  if( pars.find("y0") != end ){
    this->y0 = pars["y0"];
  }
  if( pars.find("pa") != end ){
    this->pa = (pars["pa"] + 90.0) * 0.01745329251; // in rad
    this->cospa = cos(this->pa);
    this->sinpa = sin(this->pa);
  }
  if( pars.find("upsilon_exp") != end ){
    this->upsilon_exp = pars["upsilon_exp"];
  }
  if( pars.find("q") != end ){
    this->q = pars["q"];
  }
  if( pars.find("r_eff") != end ){
    this->Reff = pars["r_eff"];
    this->sdev_fac = 2.0*this->Reff*this->Reff;
  }

  // Update M_tot (I_eff) and then I_eff (M_tot).
  if( pars.find("M_tot") != end ){
    this->M_tot = pars["M_tot"];
    double fac = M_PI*this->sdev_fac;
    this->Ieff = this->q*pow(10.0,-0.4*(this->M_tot-this->ZP))/fac;
  } else if( pars.find("i_eff") != end ){
    this->Ieff = pars["i_eff"];
    double fac = M_PI*this->sdev_fac;
    this->M_tot = -2.5*log10(fac*this->Ieff/this->q) + this->ZP;
  }
  
  set_extent();
}

double Gauss::value(double x,double y){
  if( is_in_range(x,y) ){
    double u,v,r2;
    u =   (x - this->x0)*this->cospa + (y - this->y0)*this->sinpa;
    v = - (x - this->x0)*this->sinpa + (y - this->y0)*this->cospa;
    r2 = (this->q*this->q*u*u + v*v)/this->sdev_fac;
    return this->Ieff*exp(-r2);
  } else {
    return 0.0;
  }
}

double Gauss::value_to_mass(double x,double y){
  if( is_in_range(x,y) ){
    double u,v,r2,fac,light;
    light = this->value(x,y);
    u =  (x - this->x0)*this->cospa + (y - this->y0)*this->sinpa;
    v = -(x - this->x0)*this->sinpa + (y - this->y0)*this->cospa;
    r2 = (this->q*this->q*u*u + v*v);
    fac = this->upsilon*pow(r2/this->sdev_fac,this->upsilon_exp);
    return fac*light;
  } else {
    return 0.0;
  }
}

bool Gauss::is_in_range(double xin,double yin){
  if( xin < this->p_xmin || this->p_xmax < xin || yin < this->p_ymin || this->p_ymax < yin ){
    return false;
  } else {
    return true;
  }
}

void Gauss::get_extent(double& xmin,double& xmax,double& ymin,double& ymax){
  xmin = this->p_xmin;
  xmax = this->p_xmax;
  ymin = this->p_ymin;
  ymax = this->p_ymax;
}

void Gauss::set_extent(){
  double dimg = 5.0*this->Reff;
  this->p_xmin = this->x0 - dimg;
  this->p_xmax = this->x0 + dimg;
  this->p_ymin = this->y0 - dimg;
  this->p_ymax = this->y0 + dimg;
}



//Derived class from BaseProfile: Custom
//===============================================================================================================
Custom::Custom(std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,double ZP,double M_tot,std::string interp,double upsilon): BaseProfile(1,"custom",upsilon,ZP,M_tot),RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
  this->set_interp(interp);

  // The given image must always be in units of electrons/s.
  if( M_tot == 0.0 ){
    double total_flux = this->sum(total_flux);
    this->M_tot = -2.5*log(total_flux) + ZP;
  } else {
    // If M_tot is given, then we recalibrate the given image by the total flux.
    double new_total_flux = pow(10.0,-0.4*(this->M_tot - ZP));
    double old_total_flux = this->sum();
    double factor = new_total_flux/old_total_flux;
    for(int i=0;i<this->Nz;i++){
      this->z[i] *= factor;
    }
  }

  double total_flux = this->sum();
  double total_flux_mag = -2.5*log(total_flux) + ZP;
  std::cout << "Custom source total flux is: " << total_flux << " " << total_flux_mag << " (should be " << this->M_tot << ")" << std::endl;
}

double Custom::value(double x,double y){
  if( is_in_range(x,y) ){
    double val = (this->*interp2d)(x,y,this->z);
    return val;
  } else {
    return 0.0;
  }
}

double Custom::value_to_mass(double x,double y){
  if( is_in_range(x,y) ){
    double val = (this->*interp2d)(x,y,this->z);
    return this->upsilon*val; // No radial dependence, I need some characteristic size of the light profile for that, like the R_1/2, which I don't have yet. 
  } else {
    return 0.0;
  }
}

bool Custom::is_in_range(double xin,double yin){
  if( xin < this->xmin || this->xmax < xin || yin < this->ymin || this->ymax < yin ){
    return false;
  } else {
    return true;
  }
}

void Custom::get_extent(double& x_min,double& x_max,double& y_min,double& y_max){
  x_min = this->xmin;
  x_max = this->xmax;
  y_min = this->ymin;
  y_max = this->ymax;
}





/*
//Derived class from BaseProfile: Delaunay
//===============================================================================================================
myDelaunay::myDelaunay(std::string filename){
  this->type = "delaunay";
  this->output_res = 500;


  // Read v,x,y from file to the class variables
  std::ifstream file(filename);
  std::string line;
  double xx,yy,vv;
  std::vector<double> xvec;
  std::vector<double> yvec;
  std::vector<double> vvec;

  while( std::getline(file,line) ){
    std::istringstream ss(line);
    ss >> vv >> xx >> yy;
    xvec.push_back(xx);
    yvec.push_back(yy);
    vvec.push_back(vv);
  }

  int N = xvec.size();
  this->N = N;
  this->x   = (double*) calloc(N,sizeof(double));
  this->y   = (double*) calloc(N,sizeof(double));
  this->src = (double*) calloc(N,sizeof(double));
  for(int i=0;i<N;i++){
    this->x[i]   = xvec[i];
    this->y[i]   = yvec[i];
    this->src[i] = vvec[i];
    //    std::cout << this->src[i] << " " << this->x[i] << " " << this->y[i] << std::endl;
  }

  // Create the Dealaunay triangulation
  std::vector< std::pair<Point,int> > points;
  for(int i=0;i<N;i++){
    points.push_back( std::make_pair(Point(this->x[i],this->y[i]),i) );
  }
  Delaunay triangulation;
  triangulation.insert(points.begin(),points.end());

  //Get each Delaunay triangle in my own struct, and number them
  //[constructing this->triangles]
  Delaunay::Finite_faces_iterator fit;
  Delaunay::Face_handle face;
  atriangle triangle;
  int index = 0;
  this->triangles.resize( triangulation.number_of_faces() );
  for(fit=triangulation.finite_faces_begin();fit!=triangulation.finite_faces_end();fit++){
    face = fit;
    face->info() = index;

    triangle.a = (int) face->vertex(0)->info();
    triangle.b = (int) face->vertex(1)->info();
    triangle.c = (int) face->vertex(2)->info();
    
    this->triangles[index] = triangle;
    index++;
  }


  // FILE* fh = fopen("triangles.dat","w");
  // for(int q=0;q<this->triangles.size();q++){
  //   fprintf(fh,"%10.5f%10.5f%10.5f",this->x[this->triangles[q].a],this->x[this->triangles[q].b],this->x[this->triangles[q].c]);
  //   fprintf(fh,"%10.5f%10.5f%10.5f",this->y[this->triangles[q].a],this->y[this->triangles[q].b],this->y[this->triangles[q].c]);
  //   fprintf(fh,"\n");
  // }
  // fclose(fh);


  // Get the convex hull of the triangulation
  Delaunay::Vertex_circulator vc = triangulation.incident_vertices(triangulation.infinite_vertex());
  Delaunay::Vertex_circulator done(vc);
  std::vector<Point> dum;
  do{
    dum.push_back(vc->point());
  }while( ++vc != done );

  this->ch_size = dum.size();
  this->convex_hull = (Point*) malloc(this->ch_size*sizeof(Point));
  for(int i=0;i<this->ch_size;i++){
    this->convex_hull[i] = dum[i];
    //    std::cout << dum[i].x() << " " << dum[i].y() << std::endl;
  }
  //  std::cout << std::endl << std::endl;
}

double myDelaunay::value(double x,double y){
  double val;
  // If the point is outside the convex hull of the triangulation set its value to zero, else interpolate within the triangle it in.
  if( CGAL::bounded_side_2(this->convex_hull,this->convex_hull+this->ch_size,Point(x,y),K()) == CGAL::ON_UNBOUNDED_SIDE ){
    val = 0.0;
    std::cout << "out " << x << " " << y << std::endl;
  } else {
    double wa,wb,wc;
    double ybc,xac,xcb,yac,xxc,yyc,den;
    atriangle triangle;

    for(int j=0;j<this->triangles.size();j++){
      triangle = this->triangles[j];
      
      ybc = this->y[triangle.b] - this->y[triangle.c];//(yb-yc)
      xac = this->x[triangle.a] - this->x[triangle.c];//(xa-xc)
      xcb = this->x[triangle.c] - this->x[triangle.b];//(xc-xb)
      yac = this->y[triangle.a] - this->y[triangle.c];//(ya-yc)
      xxc = x                   - this->x[triangle.c];//(x -xc)
      yyc = y                   - this->y[triangle.c];//(y -yc)
      den = ybc*xac + xcb*yac;
      
      wa = ( ybc*xxc+xcb*yyc)/den;
      wb = (-yac*xxc+xac*yyc)/den;
      wc = 1.0 - wa - wb;
      
      if( 0.0 <= wa && wa <= 1.0 && 0.0 <= wb && wb <= 1.0 && 0.0 <= wc && wc <= 1.0 ){
	val = wa*this->src[triangle.a] + wb*this->src[triangle.b] + wc*this->src[triangle.c];
	break;
      }
    }
  }

  return val;
}

void myDelaunay::outputProfile(std::string filename){
  // find the extent of the grid
  double half_range = 0.1;
  double drange = 0.1;

  double sum = 0.0;
  for(int i=0;i<this->N;i++){
    sum += this->src[i];
  }

  double part = 0.0;
  while( part/sum < 0.95 ){
    part = 0.0;
    for(int i=0;i<this->N;i++){
      if( -half_range < this->x[i] && this->x[i] < half_range && -half_range < this->y[i] && this->y[i] < half_range ){
	part += this->src[i];
      }
    }
    half_range += 0.1;
  }

  this->writeProfile(filename,half_range);  
}

double myDelaunay::sourceExtent(){
  // the extent of the source (in arcsec) must include 99% of the source flux
  double sum = 0.0;
  for(int i=0;i<this->N;i++){
    sum += this->src[i];
  }
  double limit = 0.99*sum;
  
  double psum,xmin,xmax,ymin,ymax;
  double hsize = 0.0;
  double dsize = 0.25;
  int jmax = 20;
  for(int j=1;j<jmax;j++){
    hsize = j*dsize;

    xmin = -hsize;
    xmax =  hsize;
    ymin = -hsize;
    ymax =  hsize;

    psum = 0.0;
    for(int i=0;i<this->N;i++){
      if( xmin < this->x[i] && this->x[i] < xmax && ymin < this->y[i] && this->y[i] < ymax ){
	psum += this->src[i];
      }
    }
    if( psum > limit ){
      break;
    }
  }

  if( hsize == jmax*dsize ){
    // Source is larger than jmax*dsize arcsec!
    return 0;
  }
  double size = 2*hsize; // arcsec

  return size;
}
*/
