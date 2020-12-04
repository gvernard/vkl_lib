#define _USE_MATH_DEFINES

#include "sourceProfile.hpp"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include "fitsInterface.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int,K>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2<unsigned int,K>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                          Delaunay;
typedef K::Point_2                                                     Point;

//Abstract class: BaseProfile
//===============================================================================================================
void BaseProfile::profile(int Sj,int Si,double* sx,double* sy,double* s){
  for(int i=0;i<Si;i++){
    for(int j=0;j<Sj;j++){
      s[i*Sj+j] += this->value(sx[i*Sj+j],sy[i*Sj+j]);
    }
  }
}

void BaseProfile::writeProfile(std::string filename,double half_range){
  // produces a square image centered at (0,0) on the source plane

  // create grid of source brightness profile
  double* array = (double*) malloc(output_res*output_res*sizeof(double));
  double dpix = 2.0*half_range/output_res;
  for(int i=0;i<output_res;i++){
    for(int j=0;j<output_res;j++){
      double x = j*dpix - half_range;
      double y = i*dpix - half_range;
      //array[(output_res-1-i)*output_res+j] = this->value(x,y);
      array[i*output_res+j] = this->value(x,y);
    }
  }

  // Total brightness
  double sum = 0.0;
  for(long i=0;i<output_res*output_res;i++){
    sum += array[i]*pow(2*half_range/output_res,2);
  }
  //  std::cout << sum << std::endl;
  
  //Write FITS:
  std::vector<std::string> key{"WIDTH","HEIGHT"};
  std::vector<std::string> val{std::to_string(2*half_range),std::to_string(2*half_range)};
  std::vector<std::string> txt{"width of the image in arcsec","height of the image in arcsec"};
  FitsInterface::writeFits(output_res,output_res,array,key,val,txt,filename);
  free(array);
}





//Derived class from BaseAnalyticFunction: Sersic
//===============================================================================================================
Sersic::Sersic(std::map<std::string,double> pars){
  this->type          = "sersic";
  this->pars["n"]     = pars["n"];
  this->pars["r_eff"] = pars["r_eff"];
  this->pars["q"]     = pars["q"];
  this->pars["x0"]    = pars["x0"];
  this->pars["y0"]    = pars["y0"];
  this->pars["pa"]    = pars["pa"];

  if( pars.find("M_tot") == pars.end() ){
    // not found
    this->pars["i_eff"] = pars["i_eff"];
    this->inverseScaleProfile();    
  } else {
    // found
    this->pars["M_tot"] = pars["M_tot"];
    this->scaleProfile();
  }
}

double Sersic::function_value(double x,double y){
  double bn = 1.9992*this->pars["n"] - 0.3271;//From Capaccioli 1989
  double u,v,r,fac2;
  double cosphi = cos(this->pars["pa"]*this->fac);
  double sinphi = sin(this->pars["pa"]*this->fac);
  
  u =   (x - this->pars["x0"])*cosphi + (y - this->pars["y0"])*sinphi;
  v = - (x - this->pars["x0"])*sinphi + (y - this->pars["y0"])*cosphi;
  r = sqrt(this->pars["q"]*this->pars["q"]*u*u + v*v);
  fac2 = pow(r/this->pars["r_eff"],1.0/this->pars["n"]) - 1.0;
  return this->pars["i_eff"]*exp(-bn*fac2);
}

void Sersic::scaleProfile(){
  double bn = 1.9992*this->pars["n"] - 0.3271;//From Capaccioli 1989
  double den = pow(this->pars["r_eff"],2)*2*M_PI*this->pars["n"]*exp(bn)*tgamma(2*this->pars["n"])/pow(bn,2*this->pars["n"]);
  this->pars["i_eff"] = this->pars["q"]*pow(10.0,-0.4*this->pars["M_tot"])/den;
}

void Sersic::inverseScaleProfile(){
  double bn = 1.9992*this->pars["n"] - 0.3271;//From Capaccioli 1989
  double fac = this->pars["i_eff"]*pow(this->pars["r_eff"],2)*2*M_PI*this->pars["n"]*exp(bn)*tgamma(2*this->pars["n"])/(this->pars["q"]*pow(bn,2*this->pars["n"]));
  this->pars["M_tot"] = -2.5*log10(fac);
}

std::vector<double> Sersic::extent(){
  double dx = 3*this->pars["_reff"]*cos(this->pars["pa"]*this->fac);
  double xmin = this->pars["x0"] - dx;
  double xmax = this->pars["x0"] + dx;
  double dy = 3*this->pars["r_eff"]*sin(this->pars["pa"]*this->fac);
  double ymin = this->pars["y0"] - dy;
  double ymax = this->pars["y0"] + dy;
  std::vector<double> ranges = {xmin,xmax,ymin,ymax};
  return ranges;
}

//Derived class from BaseAnalyticFunction: Gauss
//===============================================================================================================
proGauss::proGauss(std::map<std::string,double> pars){
  this->type          = "gauss";
  this->pars["r_eff"] = pars["r_eff"];
  this->pars["M_tot"] = pars["M_tot"];
  this->pars["q"]     = pars["q"];
  this->pars["x0"]    = pars["x0"];
  this->pars["y0"]    = pars["y0"];
  this->pars["pa"]    = pars["pa"];
  this->scaleProfile();
}

double proGauss::function_value(double x,double y){
  double u,v,r2;
  double cosphi = cos(this->pars["pa"]*this->fac);
  double sinphi = sin(this->pars["pa"]*this->fac);  
  double sdev   = 2*this->pars["r_eff"]*this->pars["r_eff"];
  
  u =   (x - this->pars["x0"])*cosphi + (y - this->pars["y0"])*sinphi;
  v = - (x - this->pars["x0"])*sinphi + (y - this->pars["y0"])*cosphi;
  //    u =   x*cosphi + y*sinphi;
  //    v = - x*sinphi + y*cosphi;
  r2 = (this->pars["q"]*this->pars["q"]*u*u + v*v)/sdev;
  //    return (this->ieff*exp(-r2)/(sqrt(sdev*3.14159)));
  return this->pars["i_eff"]*exp(-r2);
}

void proGauss::scaleProfile(){
  double den = 2.0*M_PI*pow(this->pars["r_eff"],2);
  this->pars["i_eff"] = pow(10.0,-0.4*this->pars["M_tot"])*this->pars["q"]/den;
}

std::vector<double> proGauss::extent(){
  double dimg = 3.0*this->pars["r_eff"];
  double xmin = this->pars["x0"] - dimg;
  double xmax = this->pars["x0"] + dimg;
  double ymin = this->pars["y0"] - dimg;
  double ymax = this->pars["y0"] + dimg;
  std::vector<double> ranges = {xmin,xmax,ymin,ymax};
  return ranges;
}







//Derived class from BaseProfile: Analytic
//===============================================================================================================
Analytic::Analytic(std::vector<std::string> names,std::vector<std::map<std::string,double> > par_maps){
  this->type = "analytic";
  this->output_res = 500;
  for(int i=0;i<names.size();i++){
    BaseAnalyticFunction* function = FactoryAnalyticFunction::getInstance()->createAnalyticFunction(names[i],par_maps[i]);
    this->components.push_back( function );
  }
}

double Analytic::value(double x,double y){
  double value = 0.0;
  for(int i=0;i<this->components.size();i++){
    value += this->components[i]->function_value(x,y);
  }
  return value;
}

void Analytic::outputProfile(std::string filename){
  double half_range = 0.0;
  for(int i=0;i<this->components.size();i++){
    std::vector<double> ranges = this->components[i]->extent();
    for(int j=0;j<ranges.size();j++){
      if( fabs(ranges[j]) > half_range ){
	half_range = fabs(ranges[j]);
      }
    }
  }
  this->writeProfile(filename,half_range);
}






//Derived class from BaseProfile: Custom
//===============================================================================================================
double Custom::value(double x,double y){
  double val = (this->*interp2d)(x,y,this->z);
  return val;
}

void Custom::outputProfile(std::string filename){
  std::vector<double> ranges = {this->xmin,this->xmax,this->ymin,this->ymax};
  double half_range = 0.0;
  for(int i=0;i<ranges.size();i++){
    if( fabs(ranges[i]) > half_range ){
      half_range = fabs(ranges[i]);
    }
  }
  this->output_res = 2.0*half_range*this->Nx/this->width;
  this->writeProfile(filename,half_range);
}

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

  /*
  FILE* fh = fopen("triangles.dat","w");
  for(int q=0;q<this->triangles.size();q++){
    fprintf(fh,"%10.5f%10.5f%10.5f",this->x[this->triangles[q].a],this->x[this->triangles[q].b],this->x[this->triangles[q].c]);
    fprintf(fh,"%10.5f%10.5f%10.5f",this->y[this->triangles[q].a],this->y[this->triangles[q].b],this->y[this->triangles[q].c]);
    fprintf(fh,"\n");
  }
  fclose(fh);
  */

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
