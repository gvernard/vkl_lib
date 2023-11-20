#include "imagePlane.hpp"

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iostream>

#include "tableDefinition.hpp"
#include "fitsInterface.hpp"
#include "rectGrid.hpp"

using namespace vkl;

//ImagePlane class implementation
//============================================================================================
ImagePlane::ImagePlane(const std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  this->Nx   = Ny;
  this->Ny   = Nx;
  grid = new RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath);

  Nm   = Nx*Ny;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;
}

ImagePlane::ImagePlane(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  this->Nx   = Nx;
  this->Ny   = Ny;
  grid = new RectGrid(Nx,Ny,xmin,xmax,ymin,ymax);

  Nm   = Nx*Ny;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;
}

ImagePlane::ImagePlane(const ImagePlane& image){
  Nx = image.Nx;
  Ny = image.Ny;
  Nm = image.Nm;
  grid = image.grid;
}

ImagePlane::~ImagePlane(){
  delete(grid);
}

void ImagePlane::writeImage(const std::string filename){  
  std::vector<std::string> key{"WIDTH","HEIGHT"};
  std::vector<std::string> val{std::to_string(this->grid->width),std::to_string(this->grid->height)};
  std::vector<std::string> txt{"width of the image in arcsec","height of the image in arcsec"};
  FitsInterface::writeFits(this->Nx,this->Ny,this->grid->z,key,val,txt,filename);
}

void ImagePlane::readS(const std::string filepath){
  if( filepath == "0" ){
    
    for(int i=0;i<this->Nm;i++){
      this->S.tri.push_back({i,i,1.});
    }
    this->Nmask = this->Nm;

  } else {

    double* fits = (double*) malloc(this->Nm*sizeof(double));
    FitsInterface::readFits(this->Nx,this->Ny,fits,filepath);
    this->Nmask = 0.0;
    for(int i=0;i<this->Nm;i++){
      if( fits[i] == 1 ){
	this->S.tri.push_back({i,i,1.0});
	this->Nmask++;
      } else {
	this->S.tri.push_back({i,i,0.0});
      }
    }
    free(fits);
    
  }
}

void ImagePlane::readC(const std::string flag,const std::string filepath){
  double value;
  this->noise_flag = flag;

  if( flag == "uniform" ){

    //There is only one element in the file, repeated along the diagonal
    std::ifstream infile(filepath);
    infile >> value;
    for(int i=0;i<this->Nm;i++){
      this->C.tri.push_back({i,i,value});
    }
    infile.close();

  } else if( flag == "map" ){

    //There are exactly Ni x Nj (diagonal) elements in the file and in table C
    std::string extension = filepath.substr(filepath.find(".")+1,std::string::npos);

    if( extension == "dat" ){

      std::ifstream infile(filepath);
      int i;
      while( true ){
	infile >> i >> i >> value;
	if( infile.eof() ) break;
	this->C.tri.push_back({i,i,value});
      }
      infile.close();

    } else if( extension == "fits" ){

      double* fits = (double*) malloc(this->Nm*sizeof(double));
      FitsInterface::readFits(this->Nx,this->Ny,fits,filepath);
      for(int i=0;i<this->Ny*this->Ny;i++){
	this->C.tri.push_back({i,i,fits[i]});
	//	  std::cout << 1./pow(contents[i*this->Nj+j],2) << std::endl;
      }
      free(fits);

    }

  } else if( flag == "correlated" ){

    //There are at least Ni x Nj (diagonal plus off-diagonal) elements in the file and in table C
    std::ifstream infile(filepath);
    int i,j;
    while( true ){
      infile >> i >> j >> value;
      if( infile.eof() ) break;
      this->C.tri.push_back({i,j,value});
    }
    infile.close();

  } else {

  }

}

void ImagePlane::readB(const std::string filepath,int i,int j,int ci,int cj){
  if( filepath == "0" ){

    //Create a diagonal matrix of the given dimensions.
    for(int i=0;i<this->B.Ti;i++){
      this->B.tri.push_back({i,i,1.});
    }

  } else {

    int Pi     = i;
    int Pj     = j;
    int Ncropx = ci;
    int Ncropy = cj;

    if( Pi%2 == 0 ){
      if( ci%2 != 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH EVEN" << std::endl << std::endl;
      }
    } else {
      if( ci%2 == 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH ODD" << std::endl << std::endl;
      }
    }
    if( Pj%2 == 0 ){
      if( cj%2 != 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH EVEN" << std::endl << std::endl;
      }
    } else {
      if( cj%2 == 0 ){
	std::cout << std::endl << std::endl << "PROBLEM!!! PSF AND CROPPED PSF DIMENSIONS MUST BE BOTH ODD" << std::endl << std::endl;
      }
    }

    // Read PSF from file
    double* fits = (double*) malloc(this->Nm*sizeof(double));
    FitsInterface::readFits(this->Nx,this->Ny,fits,filepath);

    // Report on the peak location of the PSF
    double vmax = fits[0];
    double imax = 0;
    for(int i=1;i<Pj*Pi;i++){
      if( fits[i] > vmax ){
	vmax = fits[i];
	imax = i;
      }
    }

    if( Pi%2 == 0 && Pj%2 == 0 ){
      int Pi2 = Pi/2;
      int Pj2 = Pj/2;
      int indices[] = {(Pi2-1)*Pj+Pj2-1,(Pi2-1)*Pj+Pj2,Pi2*Pj+Pj2-1,Pi2*Pj+Pj2};
      int loc[]     = {0,0,0,0};
      bool flag     =  true;
      for(int i=0;i<4;i++){
	if( indices[i] == imax ){
	  loc[i] = 1;
	  flag = false;
	}
      }
      if( flag ){
	//	std::cout << "Particular form of the PSF: it is even-even but the peak lies outside the four central pixels" << std::endl;
      } else {
	//	std::cout << std::endl << "PSF is even-even with the peak located at: " << std::endl;
	//	printf("%2d%2d\n%2d%2d\n",loc[0],loc[1],loc[2],loc[3]);
      }
    } else if( Pi%2 != 0 && Pj%2 != 0 ){
      int Pi2 = (Pi-1)/2;
      int Pj2 = (Pj-1)/2;
      int centre = Pi2*Pj + Pj2;
      if( centre != imax ){
	//	std::cout << "Particular form of the PSF: it is odd-odd but the peak is not the central pixel" << std::endl;	
      }
    } else {
      //      std::cout << "Particular form of the PSF: it it neither even-even nor odd-odd" << std::endl;
    }

    // Set cropped PSF and normalize
    double* blur = (double*) calloc(Ncropx*Ncropy,sizeof(double));
    int offset_y = (Pi - Ncropy)/2;
    int offset_x = (Pj - Ncropx)/2;
    int offset_tot = (offset_y)*Pj + offset_x;
    double sum = 0.0;
    for(int i=0;i<Ncropy;i++){
      for (int j=0;j<Ncropx;j++){
	blur[i*Ncropx+j] = fits[offset_tot+i*Pj+j];
	sum += blur[i*Ncropx+j];
      }
    }
    free(fits);
    double fac = 1.0/sum;
    for(int i=0;i<Ncropx*Ncropy;i++){
      blur[i] *= fac;
    }    
    
    // Set quad-kernel and odd/even functions
    int Nquadx,Nquady;
    void (ImagePlane::*xfunc)(int,int,int,int,int&,int&,int&);
    void (ImagePlane::*yfunc)(int,int,int,int,int&,int&,int&);
    if( Ncropx%2 == 0){
      Nquadx = Ncropx/2;
      xfunc = &ImagePlane::setCroppedLimitsEven;
    } else {
      Nquadx = ceil(Ncropx/2.0);
      xfunc = &ImagePlane::setCroppedLimitsOdd;
    }
    if( Ncropy%2 == 0){
      Nquady = Ncropy/2;
      yfunc = &ImagePlane::setCroppedLimitsEven;
    } else {
      Nquady = ceil(Ncropy/2.0);
      yfunc = &ImagePlane::setCroppedLimitsOdd;
    }

    // Create the blurring matrix
    int Nleft,Nright,Ntop,Nbottom,crop_offsetx,crop_offsety,crop_offset;
    int ic,jc;
    double val;
    for(int i=0;i<this->Ny;i++){
      for(int j=0;j<this->Nx;j++){

	(this->*(xfunc))(j,Ncropx,this->Nx,Nquadx,Nleft,Nright,crop_offsetx);
	(this->*(yfunc))(i,Ncropy,this->Ny,Nquady,Ntop,Nbottom,crop_offsety);
	crop_offset = crop_offsety*Ncropx + crop_offsetx;

	for(int ii=i-Ntop;ii<i+Nbottom;ii++){
	  ic = ii - i + Ntop;
	  for(int jj=j-Nleft;jj<j+Nright;jj++){
	    jc = jj - j + Nleft;
  
	    val = blur[crop_offset + ic*Ncropx + jc];

	    this->B.tri.push_back({i*this->Nx+j,     ii*this->Nx+jj,     val });
	  }
	}

      }
    }

    free(blur);
  }
}
  
void ImagePlane::setCroppedLimitsEven(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset){
  if( k < Nquad ){
    Npre   = k;
    Npost  = Nquad;
    offset = Nquad - k;
  } else if( k > (Nimg - Nquad - 1) ){
    Npre   = Nquad;
    Npost  = Nimg - k;
    offset = 0;
  } else {
    Npre   = Nquad;
    Npost  = Nquad;
    offset = 0;
  }
}

void ImagePlane::setCroppedLimitsOdd(int k,int Ncrop,int Nimg,int Nquad,int &Npre,int &Npost,int& offset){
  if( k < (Nquad - 1) ){
    Npre   = k;
    Npost  = Nquad;
    offset = Nquad - 1 - k;
  } else if( k > (Nimg - Nquad - 1) ){
    Npre   = Nquad - 1;
    Npost  = Nimg - k;
    offset = 0;
  } else {
    Npre   = Nquad-1;
    Npost  = Nquad;
    offset = 0;
  }
}




/*
void ImagePlane::maskData(std::map<int,int> lookup,ImagePlane* masked){
  typedef std::map<int,int>::iterator it;
  for(it myiterator=lookup.begin();myiterator!=lookup.end();myiterator++){
    masked->img[ myiterator->second ] = this->img[ myiterator->first ];
    masked->x[ myiterator->second ]   = this->x[ myiterator->first ];
    masked->y[ myiterator->second ]   = this->y[ myiterator->first ];
  }
}
*/
