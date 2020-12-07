#include "sourcePlane.hpp"

#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "fitsInterface.hpp"
#include "covKernels.hpp"

//Derived class from BaseSourcePlane: FixedSource
//===============================================================================================================
FixedSource::FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,BaseCovKernel* kernel): BaseSourcePlane(kernel),RectGrid(Nx,Ny,xmin,xmax,ymin,ymax) {
  source_type = "fixed";
  Sm   = this->Nz;
}

FixedSource::FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string filepath,BaseCovKernel* kernel): BaseSourcePlane(kernel),RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath){
  source_type = "fixed";
  Sm   = this->Nz;
}

FixedSource::FixedSource(const FixedSource& other): BaseSourcePlane(other), RectGrid(other){
  source_type = "fixed";
};


//virtual
FixedSource* FixedSource::clone(){
  return new FixedSource(*this);
};

//virtual
void FixedSource::constructH(std::string reg_scheme){
  mytable newH;
  newH.Ti = this->Sm;
  newH.Tj = this->Sm;

  if( reg_scheme == "identity" ){

    this->eigenSparseMemoryAllocForH = 1;
    for(int i=0;i<newH.Ti;i++){
      newH.tri.push_back({i,i,1});
    }
    
  } else if ( reg_scheme == "gradient" ){

    this->eigenSparseMemoryAllocForH = 3;
    int i0,j0;
    double dx = this->step_x;
    double dy = this->step_y;
    std::vector<mytriplet> c;

    for(i0=0;i0<this->Ny-1;i0++){
      for(j0=0;j0<this->Nx-1;j0++){
	c = this->getChunk(i0,j0,dx,dy,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff);
	appendToH(newH,c);
      }
      // Last element in x:
      j0 = this->Nx-1;
      c = this->getChunk(i0,j0,dx,dy,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff);
      appendToH(newH,c);
    }

    // Last row in y
    i0 = this->Ny-1;
    for(j0=0;j0<this->Nx-1;j0++){
      c = this->getChunk(i0,j0,dx,dy,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff);
      appendToH(newH,c);
    }
    // Last element in x:
    j0 = this->Nx-1;
    c = this->getChunk(i0,j0,dx,dy,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff);
    appendToH(newH,c);
    
  } else if( reg_scheme == "curvature_full" ){

    this->eigenSparseMemoryAllocForH = 5;
    int i0,j0;
    double dx2 = pow(this->step_x,2);
    double dy2 = pow(this->step_y,2);
    std::vector<mytriplet> c;
    
    // First Row
    i0 = 0;
    j0 = 0;
    c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff);
    appendToH(newH,c);
    for(j0=1;j0<this->Nx-1;j0++){
      c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff);
      appendToH(newH,c);
    }
    j0 = this->Nx-1;
    c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff);
    appendToH(newH,c);
    
    // Middle rows
    for(i0=1;i0<this->Ny-1;i0++){
      j0 = 0;
      c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff);
      appendToH(newH,c);
      for(int j0=1;j0<this->Nx-1;j0++){
	c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff);
	appendToH(newH,c);
      }
      j0 = this->Nx-1;
      c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff);
      appendToH(newH,c);
    }

    // Last Row
    i0 = this->Ny-1;
    j0 = 0;
    c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_forward_1_index,Constants::derivative_2_forward_1_coeff,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff);
    appendToH(newH,c);
    for(j0=1;j0<this->Nx-1;j0++){
      c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff);
      appendToH(newH,c);
    }
    j0 = this->Nx-1;
    for(int k=0;k<Constants::derivative_2_backward_1_index.size();k++){
      c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff,Constants::derivative_2_backward_1_index,Constants::derivative_2_backward_1_coeff);
      appendToH(newH,c);
    }

  } else if( reg_scheme == "curvature" ){

    this->eigenSparseMemoryAllocForH = 5;
    int i0,j0;
    double dx2 = pow(this->step_x,2);
    double dy2 = pow(this->step_y,2);
    std::vector<mytriplet> c;
    
    // First Row
    i0 = 0;
    for(j0=0;j0<this->Nx;j0++){
      newH.tri.push_back({i0*this->Nx+j0,i0*this->Nx+j0,1.0});
    }    
    // Middle rows
    for(i0=1;i0<this->Ny-1;i0++){
      newH.tri.push_back({i0*this->Nx,i0*this->Nx,1.0});      
      for(int j0=1;j0<this->Nx-1;j0++){
	c = this->getChunk(i0,j0,dx2,dy2,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff,Constants::derivative_2_central_2_index,Constants::derivative_2_central_2_coeff);
	appendToH(newH,c);
      }
      newH.tri.push_back({i0*this->Nx+Nx-1,i0*this->Nx+Nx-1,1.0});
    }
    // Last Row
    i0 = this->Ny-1;
    for(j0=0;j0<this->Nx;j0++){
      newH.tri.push_back({i0*this->Nx+j0,i0*this->Nx+j0,1.0});
    }
    
  } else if( reg_scheme == "covariance_kernel" ){

    if( this->kernel == NULL ){
      // throw an exception
    }
    int* nonZeroRow = (int*) calloc(this->Sm,sizeof(int));
    int index1,index2;
    double x1,y1,x2,y2;
    double cov,r;
    // Match each grid point to every other point through this double loop
    for(int i1=0;i1<this->Ny-1;i1++){
      for(int j1=0;j1<this->Nx-1;j1++){
	x1 = this->center_x[j1];
	y1 = this->center_y[i1];
	index1 = i1*this->Nx+j1;

	cov = this->kernel->getCovarianceSelf();
	if( this->kernel->larger_than_cmax(cov) ){
	  newH.tri.push_back({index1,index1,cov});
	  nonZeroRow[index1]++;
	}	

	for(int i2=i1+1;i2<this->Ny;i2++){
	  for(int j2=j1+1;j2<this->Nx;j2++){
	    x2 = this->center_x[j2];
	    y2 = this->center_y[i2];
	    index2 = i2*this->Nx+j2;
	    
	    r = hypot(x1-x2,y1-y2);
	    cov = this->kernel->getCovariance(r);
	    if( this->kernel->larger_than_cmax(cov) ){
	      newH.tri.push_back({index1,index2,cov});
	      nonZeroRow[index1]++;
	    }
	  }
	}
      }
    }

    // Find maximum number of non-zero elements per row
    int maxNonZero = nonZeroRow[0];
    for(int i=1;i<this->Sm;i++){
      if( nonZeroRow[i] > maxNonZero ){
	maxNonZero = nonZeroRow[i];
      }
    }
    free(nonZeroRow);
    this->eigenSparseMemoryAllocForH = maxNonZero;

  }

  this->H.insert( std::pair<std::string,mytable> (reg_scheme,newH) );
}

//virtual
void FixedSource::outputSource(const std::string fname){
  std::vector<std::string> key{"WIDTH","HEIGHT"};
  std::vector<std::string> val{std::to_string(this->width),std::to_string(this->height)};
  std::vector<std::string> txt{"width of the source in arcsec","height of the source in arcsec"};
  FitsInterface::writeFits(this->Nx,this->Ny,this->z,key,val,txt,fname);
}


//non-virtual private
std::vector<mytriplet> FixedSource::getChunk(int i0,int j0,double hx,double hy,std::vector<int> ind_x,std::vector<double> coeff_x,std::vector<int> ind_y,std::vector<double> coeff_y){
  std::map<int,mytriplet> triplet_map;
  int ind;
  int ind0 = i0*this->Nx + j0;
  mytriplet tri;
  
  // first add the x coefficients
  for(int k=0;k<ind_x.size();k++){
    ind = i0*this->Nx+(j0+ind_x[k]);
    tri = {ind0,ind,coeff_x[k]/hx};
    triplet_map.insert(std::pair<int,mytriplet> (ind,tri));
  }
  // now add the y coefficients and check if the same i,j exists (the central pixel)
  for(int k=0;k<ind_y.size();k++){
    ind = (i0+ind_y[k])*this->Nx+j0;
    tri = {ind0,ind,coeff_y[k]/hy};
    if( triplet_map.find(ind) == triplet_map.end() ){
      triplet_map.insert(std::pair<int,mytriplet> (ind,tri));
    } else {
      triplet_map[ind].v += coeff_y[k];
    }
  }
  // check for zero entries
  std::vector<mytriplet> out;
  for(std::map<int,mytriplet>::iterator it=triplet_map.begin();it!=triplet_map.end();it++){
    if( it->second.v != 0 ){
      out.push_back(it->second);
    }
  } 
  
  return out;
}

//non-virtual private
void FixedSource::appendToH(mytable& H,std::vector<mytriplet>& chunk){
  for(int k=0;k<chunk.size();k++){
    H.tri.push_back(chunk[k]);
  }
}
