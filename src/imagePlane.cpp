#include "imagePlane.hpp"

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iostream>

#include "tableDefinition.hpp"
#include "fitsInterface.hpp"
#include "rectGrid.hpp"

//ImagePlane class implementation
//============================================================================================
ImagePlane::ImagePlane(const std::string filepath,int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  Ni   = Ny;
  Nj   = Nx;
  grid = new RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath);

  Nm   = Ni*Nj;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  
  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      cells[i*Nj+j] = NULL;
      crosses[i*Nj+j] = NULL;
      dpsi_cells[i*Nj+j] = NULL;
    }
  }
}

ImagePlane::ImagePlane(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax){
  Ny   = Ny;
  Nx   = Nx;
  grid = new RectGrid(Nx,Ny,xmin,xmax,ymin,ymax);

  Nm   = Ni*Nj;
  S.Ti = Nm;
  S.Tj = Nm;
  C.Ti = Nm;
  C.Tj = Nm;
  B.Ti = Nm;
  B.Tj = Nm;

  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));

  for(int i=0;i<Ni;i++){
    for(int j=0;j<Nj;j++){
      cells[i*Nj+j] = NULL;
      crosses[i*Nj+j] = NULL;
      dpsi_cells[i*Nj+j] = NULL;
    }
  }
}

ImagePlane::ImagePlane(const ImagePlane& image){
  Ni = image.Ni;
  Nj = image.Nj;
  Nm = image.Nm;
  grid = image.grid;
  
  defl_x  = (double*) calloc(Nm,sizeof(double));
  defl_y  = (double*) calloc(Nm,sizeof(double));
  active  = (int*) calloc(Nm,sizeof(int));
  cells   = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  crosses = (Cross**) malloc(Nm*sizeof(Cross*));
  dpsi_cells = (InterpolationCell**) calloc(Nm,sizeof(InterpolationCell*));
  for(int i=0;i<Nm;i++){
    cells[i] = NULL;
    crosses[i] = NULL;
    dpsi_cells[i] = NULL;
  }
}

ImagePlane::~ImagePlane(){
  delete(grid);
  free(defl_x);
  free(defl_y);
  free(active);
  for(int i=0;i<this->Nm;i++){
    delete(cells[i]);
    delete(crosses[i]);
    delete(dpsi_cells[i]);
  }
  free(cells);
  free(crosses);
  free(dpsi_cells);
}

void ImagePlane::writeImage(const std::string filename){  
  std::vector<std::string> key{"WIDTH","HEIGHT"};
  std::vector<std::string> val{std::to_string(this->grid->width),std::to_string(this->grid->height)};
  std::vector<std::string> txt{"width of the image in arcsec","height of the image in arcsec"};
  FitsInterface::writeFits(this->Ni,this->Nj,this->grid->z,key,val,txt,filename);
}

void ImagePlane::readS(const std::string filepath){
  if( filepath == "0" ){
    
    for(int i=0;i<this->Nm;i++){
      this->S.tri.push_back({i,i,1.});
    }
    this->Nmask = this->Nm;

  } else {

    double* fits = (double*) malloc(this->Nm*sizeof(double));
    FitsInterface::readFits(this->Ni,this->Nj,fits,filepath);
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
      FitsInterface::readFits(this->Ni,this->Nj,fits,filepath);
      for(int i=0;i<this->Ni*this->Nj;i++){
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
    FitsInterface::readFits(this->Ni,this->Nj,fits,filepath);

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
    for(int i=0;i<this->Ni;i++){
      for(int j=0;j<this->Nj;j++){

	(this->*(xfunc))(j,Ncropx,this->Nj,Nquadx,Nleft,Nright,crop_offsetx);
	(this->*(yfunc))(i,Ncropy,this->Ni,Nquady,Ntop,Nbottom,crop_offsety);
	crop_offset = crop_offsety*Ncropx + crop_offsetx;

	for(int ii=i-Ntop;ii<i+Nbottom;ii++){
	  ic = ii - i + Ntop;
	  for(int jj=j-Nleft;jj<j+Nright;jj++){
	    jc = jj - j + Nleft;
  
	    val = blur[crop_offset + ic*Ncropx + jc];

	    this->B.tri.push_back({i*this->Nj+j,     ii*this->Nj+jj,     val });
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

//void ImagePlane::printCross(int k,mytable Ds){
void ImagePlane::printCross(int k){
  double coeffs[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double dsx = 1.0;//Ds.tri[2*k].v;
  double dsy = 1.0;//Ds.tri[2*k+1].v;

  coeffs[0]  = this->crosses[k]->coeff_y[0]*dsy;
  coeffs[1]  = this->crosses[k]->coeff_y[4]*dsy;
  coeffs[2]  = this->crosses[k]->coeff_x[0]*dsx;
  coeffs[3]  = this->crosses[k]->coeff_y[1]*dsy + this->crosses[k]->coeff_x[1]*dsx;
  coeffs[4]  = this->crosses[k]->coeff_y[5]*dsy + this->crosses[k]->coeff_x[2]*dsx;
  coeffs[5]  = this->crosses[k]->coeff_x[3]*dsx;
  coeffs[6]  = this->crosses[k]->coeff_x[4]*dsx;
  coeffs[7]  = this->crosses[k]->coeff_y[2]*dsy + this->crosses[k]->coeff_x[5]*dsx;
  coeffs[8]  = this->crosses[k]->coeff_y[6]*dsy + this->crosses[k]->coeff_x[6]*dsx;
  coeffs[9]  = this->crosses[k]->coeff_x[7]*dsx;
  coeffs[10] = this->crosses[k]->coeff_y[3]*dsy;
  coeffs[11] = this->crosses[k]->coeff_y[7]*dsy;

  printf("\n");
  printf("%10s %10.3f %10.3f %10s\n"," ",coeffs[0],coeffs[1]," ");
  printf("%10.3f %10.3f %10.3f %10.3f\n",coeffs[2],coeffs[3],coeffs[4],coeffs[5]);
  printf("%10.3f %10.3f %10.3f %10.3f\n",coeffs[6],coeffs[7],coeffs[8],coeffs[9]);
  printf("%10s %10.3f %10.3f %10s\n"," ",coeffs[10],coeffs[11]," ");
}

void ImagePlane::lowerResRebinAdditive(ImagePlane* newImage){
  double inf_dx = this->grid->width/this->Nj;
  double inf_dy = this->grid->height/this->Ni;
  double new_dx = newImage->grid->width/newImage->Nj;
  double new_dy = newImage->grid->height/newImage->Ni;
  for(int i=0;i<this->Ni;i++){
    int ii = (int) floor(i*inf_dy/new_dy);
    for(int j=0;j<this->Nj;j++){
      int jj = (int) floor(j*inf_dx/new_dx);
      newImage->grid->z[ii*newImage->Nj + jj] += this->grid->z[i*this->Nj + j];
    }
  }
}

void ImagePlane::lowerResRebinIntegrate(ImagePlane* newImage){
  // Integrating over equally sized elements is equivalent to getting the average
  int* counts = (int*) calloc(newImage->Nm,sizeof(int));
  double inf_dx = this->grid->width/this->Nj;
  double inf_dy = this->grid->height/this->Ni;
  double new_dx = newImage->grid->width/newImage->Nj;
  double new_dy = newImage->grid->height/newImage->Ni;
  for(int i=0;i<this->Ni;i++){
    int ii = (int) floor(i*inf_dy/new_dy);
    for(int j=0;j<this->Nj;j++){
      int jj = (int) floor(j*inf_dx/new_dx);
      newImage->grid->z[ii*newImage->Nj + jj] += this->grid->z[i*this->Nj + j];
      counts[ii*newImage->Nj + jj] += 1;
    }
  }
  for(int i=0;i<newImage->Nm;i++){
    newImage->grid->z[i] = newImage->grid->z[i]/counts[i];
  }
  free(counts);
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
