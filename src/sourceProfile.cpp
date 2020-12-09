#define _USE_MATH_DEFINES

#include "sourceProfile.hpp"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>

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



//Derived class from BaseAnalyticFunction: Sersic
//===============================================================================================================
Sersic::Sersic(std::map<std::string,double> pars): BaseProfile(7,"sersic"){
  updateProfilePars(pars);
}

void Sersic::updateProfilePars(std::map<std::string,double> pars){
  // First modify all parameter values that come from the outside
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    ppars[it->first] = it->second;
  }
  // Then convert any of the following accordingly, if they exist in the incoming pars
  if( pars.find("n") != pars.end() ){
    this->bn = 1.9992*this->ppars["n"] - 0.3271;//From Capaccioli 1989    
  }
  if( pars.find("pa") != pars.end() ){
    ppars["pa"] = (pars["pa"] + 90.0) * 0.01745329251; // in rad
    this->cospa = cos(ppars["pa"]);
    this->sinpa = sin(ppars["pa"]);
  }
  if( pars.find("M_tot") != pars.end() ){
    double fac = pow(ppars["r_eff"],2)*2*M_PI*ppars["n"]*exp(this->bn)*tgamma(2*ppars["n"])/pow(this->bn,2*ppars["n"]);
    ppars["i_eff"] = ppars["q"]*pow(10.0,-0.4*ppars["M_tot"])/fac;
  } else {
    double fac = pow(ppars["r_eff"],2)*2*M_PI*ppars["n"]*exp(this->bn)*tgamma(2*ppars["n"])/pow(this->bn,2*ppars["n"]);
    ppars["M_tot"] = -2.5*log10(fac*ppars["i_eff"]/ppars["q"]);
  }
  set_extent();
}

double Sersic::value(double x,double y){
  if( in_range(x,y) ){
    double u,v,r,fac2;
    u =  (x - ppars["x0"])*this->cospa + (y - ppars["y0"])*this->sinpa;
    v = -(x - ppars["x0"])*this->sinpa + (y - ppars["y0"])*this->cospa;
    r = hypot(ppars["q"]*u,v);
    fac2 = pow(r/ppars["r_eff"],1.0/ppars["n"]) - 1.0;
    return ppars["i_eff"]*exp(-this->bn*fac2);
  } else {
    return 0.0;
  }
}

bool Sersic::in_range(double xin,double yin){
  if( xin < this->p_xmin || this->p_xmax < xin || yin < this->p_ymin || this->p_ymax < yin ){
    return false;
  } else {
    return true;
  }
}

void Sersic::set_extent(){
  double dx = 3*ppars["_reff"]*this->cospa;
  this->p_xmin = ppars["x0"] - dx;
  this->p_xmax = ppars["x0"] + dx;
  double dy = 3*ppars["r_eff"]*this->sinpa;
  this->p_ymin = ppars["y0"] - dy;
  this->p_ymax = ppars["y0"] + dy;
}

//Derived class from BaseAnalyticFunction: Gauss
//===============================================================================================================
Gauss::Gauss(std::map<std::string,double> pars): BaseProfile(6,"gauss"){
  updateProfilePars(pars);
}

void Gauss::updateProfilePars(std::map<std::string,double> pars){
  // First modify all parameter values that come from the outside
  for(std::map<std::string,double>::iterator it=pars.begin();it!=pars.end();it++){
    ppars[it->first] = it->second;
  }
  // Then convert any of the following accordingly, if they exist in the incoming pars
  if( pars.find("pa") != pars.end() ){
    ppars["pa"] = (pars["pa"] + 90.0) * 0.01745329251; // in rad
    this->cospa = cos(ppars["pa"]);
    this->sinpa = sin(ppars["pa"]);
  }
  if( pars.find("r_eff") != pars.end() ){
    this->sdev = 2*ppars["r_eff"]*ppars["r_eff"];
  }
  if( pars.find("M_tot") != pars.end() ){
    double fac = 2.0*M_PI*pow(ppars["r_eff"],2);
    ppars["i_eff"] = ppars["q"]*pow(10.0,-0.4*ppars["M_tot"])/fac;
  } else {
    double fac = 2.0*M_PI*pow(ppars["r_eff"],2);
    ppars["M_tot"] = -2.5*log10(fac*ppars["i_eff"]/ppars["q"]);
  }
  set_extent();
}

double Gauss::value(double x,double y){
  if( in_range(x,y) ){
    double u,v,r2;
    u =   (x - ppars["x0"])*this->cospa + (y - ppars["y0"])*sinpa;
    v = - (x - ppars["x0"])*this->sinpa + (y - ppars["y0"])*cospa;
    r2 = (ppars["q"]*ppars["q"]*u*u + v*v)/this->sdev;
    return ppars["i_eff"]*exp(-r2);
  } else {
    return 0.0;
  }
}

bool Gauss::in_range(double xin,double yin){
  if( xin < this->p_xmin || this->p_xmax < xin || yin < this->p_ymin || this->p_ymax < yin ){
    return false;
  } else {
    return true;
  }
}

void Gauss::set_extent(){
  double dimg = 3.0*ppars["r_eff"];
  this->p_xmin = ppars["x0"] - dimg;
  this->p_xmax = ppars["x0"] + dimg;
  this->p_ymin = ppars["y0"] - dimg;
  this->p_ymax = ppars["y0"] + dimg;
}



//Derived class from BaseProfile: Custom
//===============================================================================================================
Custom::Custom(std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,double Mtot,std::string interp): BaseProfile(0,"custom"),RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
  this->Mtot = Mtot;
  this->set_interp(interp);
  scaleProfile();
}

double Custom::value(double x,double y){
  if( in_range(x,y) ){
    double val = (this->*interp2d)(x,y,this->z);
    return val;
  } else {
    return 0.0;
  }
}

bool Custom::in_range(double xin,double yin){
  if( xin < this->xmin || this->xmax < xin || yin < this->ymin || this->ymax < yin ){
    return false;
  } else {
    return true;
  }
}

// private
void Custom::scaleProfile(){
  double sum = 0.0;
  double dS = this->step_x*this->step_y;
  for(int i=0;i<this->Nz;i++){
    sum += this->z[i]*dS;
  }
  double factor = pow(10.0,-0.4*this->Mtot)/sum;
  for(int i=0;i<this->Nz;i++){
    this->z[i] *= factor;
  }
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
