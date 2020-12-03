#include "sourcePlane.hpp"
#include "tableDefinition.hpp"

#include <vector>
#include <string>

//Derived class from BaseSourcePlane: FixedSource
//===============================================================================================================
FixedSource::FixedSource(int Nx,int Ny,double xmin,double xmax,double ymin,double ymax,std::string reg_scheme){
  type = "fixed";
  reg  = reg_scheme;
  grid = new RectGrid(Nx,Ny,xmin,xmax,ymin,ymax,filepath);
  Sm   = grid->Nz;
  H.Ti = Sm;
  H.Tj = Sm;
}

FixedSource::FixedSource(const FixedSource& source){
  type = "fixed";
  reg  = source.reg;
  grid = source.grid;
  Sm   = source.Sm;
}

//virtual
void FixedSource::constructH(){
  std::vector<mytriplet> tmp;//need to make sure that the H triplet vector is a new one

  if( this->reg == "identity" ){//---------------------------> zero order
    this->eigenSparseMemoryAllocForH = 1;
    for(int i=0;i<this->H.Ti;i++){
      tmp.push_back({i,i,1});
    }
  } else if ( this->reg == "gradient" ){//-------------------> first order
    
    // This is a combination of Tikhonov and gradient regularization.
    // Apart from the edge pixels that have a unit weight, there is also the central pixel of each row that has a unit weight.
    // For gradient regularization only (of any accuracy) the weight of the central pixel is zero.

    this->eigenSparseMemoryAllocForH = 8;
    int Si = this->Si;
    int Sj = this->Sj;
    double dx = this->x[1] - this->x[0];
    double dy = this->y[this->Sj] - this->y[0];
    double ddx =  1.0/dx;
    double ddy =  1.0/dy;

    /*
    //First pixel of first image row:  forward 1st derivative in X, forward 1st derivative in Y
    tmp.push_back({  0,    0,      -ddx - ddy   });
    tmp.push_back({  0,    1,             ddx   });
    tmp.push_back({  0,    1*Sj,          ddy   });
    //First row of image pixels: central 1st derivative in X, forward 1st derivative in Y
    for(int j=1;j<Sj-1;j++){
	tmp.push_back({  j,    j-1,           -0.5*ddx   });
	tmp.push_back({  j,    j,                 -ddy   });
	tmp.push_back({  j,    j+1,            0.5*ddx   });
	tmp.push_back({  j,    1*Sj+j,             ddy   });
    }
    //Last pixel of first image row:  backward 1st derivative in X, forward 1st derivative in Y
    tmp.push_back({  Sj-1,    Sj-2,             -ddx   });
    tmp.push_back({  Sj-1,    Sj-1,        ddx - ddy   });
    tmp.push_back({  Sj-1,    1*Sj+Sj-1,         ddy   });
    */    
    //Set all first-row pixels to unity
    for(int j=0;j<Sj;j++){
      tmp.push_back({         j,       j,               1.0   });
    }


    for(int i=1;i<Si-1;i++){
      /*
      //First pixel of each image row: forward 1st derivative in X-direction, central 1st derivative in Y
      tmp.push_back({  i*Sj,    (i-1)*Sj,          -0.5*ddy   });
      tmp.push_back({  i*Sj,        i*Sj,              -ddx   });
      tmp.push_back({  i*Sj,      i*Sj+1,               ddx   });
      tmp.push_back({  i*Sj,    (i+1)*Sj,           0.5*ddy   });
      */
      //First pixel of each image row: unity
      tmp.push_back({  i*Sj,        i*Sj,               1.0   });

      //central 2nd derivative in both X and Y directions
      for(int j=1;j<Sj-1;j++){
	tmp.push_back({  i*Sj+j,    (i-1)*Sj+j,    -0.5*ddy   });
	tmp.push_back({  i*Sj+j,      i*Sj+j-1,    -0.5*ddx   });
	tmp.push_back({  i*Sj+j,        i*Sj+j,     1.0       });
	tmp.push_back({  i*Sj+j,      i*Sj+j+1,     0.5*ddx   });
	tmp.push_back({  i*Sj+j,    (i+1)*Sj+j,     0.5*ddy   });
      }
      /*
      //Last pixel of each image row: backward 1st derivative in X-direction, central 1st derivative in Y
      tmp.push_back({  i*Sj+Sj-1,    (i-1)*Sj+Sj-1,          -0.5*ddy   });
      tmp.push_back({  i*Sj+Sj-1,        i*Sj+Sj-2,              -ddx   });
      tmp.push_back({  i*Sj+Sj-1,        i*Sj+Sj-1,               ddx   });
      tmp.push_back({  i*Sj+Sj-1,    (i+1)*Sj+Sj-1,           0.5*ddy   });
      */
      //Last pixel of each image row: unity
      tmp.push_back({  i*Sj+Sj-1,        i*Sj+Sj-1,     1.0   });
    }

    /*
    //First pixel of last image row:  forward 1st derivative in X, backward 1st derivative in Y
    tmp.push_back({  (Si-1)*Sj,     (Si-2)*Sj,         -ddy   });
    tmp.push_back({  (Si-1)*Sj,     (Si-1)*Sj,   -ddx + ddy   });
    tmp.push_back({  (Si-1)*Sj,   (Si-1)*Sj+1,          ddx   });
    //Last row of image pixels:  central 1st derivative in X, backward 1st derivative in Y
    for(int j=1;j<Sj-1;j++){
      tmp.push_back({  (Si-1)*Sj+j,     (Si-2)*Sj+j,              -ddy   });
      tmp.push_back({  (Si-1)*Sj+j,   (Si-1)*Sj+j-1,          -0.5*ddx   });
      tmp.push_back({  (Si-1)*Sj+j,     (Si-1)*Sj+j,               ddy   });
      tmp.push_back({  (Si-1)*Sj+j,   (Si-1)*Sj+j+1,           0.5*ddx   });
    }
    //Last pixel of last image row:  backward 1st derivative in X, backward 1st derivative in Y
    tmp.push_back({  (Si-1)*Sj+Sj-1,   (Si-2)*Sj+Sj-1,        -ddy   });
    tmp.push_back({  (Si-1)*Sj+Sj-1,   (Si-1)*Sj+Sj-2,        -ddx   });
    tmp.push_back({  (Si-1)*Sj+Sj-1,   (Si-1)*Sj+Sj-1,   ddx + ddy   });
    */    
    //Set all last-row pixels to unity
    for(int j=0;j<Sj;j++){
      tmp.push_back({  (Si-1)*Sj+j,    (Si-1)*Sj+j,     1.0   });
    }

  } else if( this->reg == "curvature_full" ){//-------------------> second order

    this->eigenSparseMemoryAllocForH = 8;
    int i0,j0;
    double ddx2 = pow(this->grid->step_x,2);
    double ddy2 = pow(this->grid->step_y,2);
    
    // First Row
    i0 = 0;
    j0 = 0;
    for(int k=0;k<RectGrid::derivative_2_forward_1_index.size();k++){
      tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_forward_1_index[k]),RectGrid::derivative_2_forward_1_coeff[k]/ddx2});
      tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_forward_1_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_forward_1_coeff[k]/ddy2});
    }

    for(j0=1;j0<this->grid->Nx-1;j0++){
      for(int k=0;k<RectGrid::derivative_2_central_2_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_central_2_index[k]),RectGrid::derivative_2_central_2_coeff[k]/ddx2});
      }
      for(int k=0;k<RectGrid::derivative_2_forward_1_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_forward_1_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_forward_1_coeff[k]/ddy2});
      }
    }

    j0 = this->grid->Nx-1;
    for(int k=0;k<RectGrid::derivative_2_backward_1_index.size();k++){
      tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_backward_1_index[k]),RectGrid::derivative_2_backward_1_coeff[k]/ddx2});
    }
    for(int k=0;k<RectGrid::derivative_2_forward_1_index.size();k++){
      tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_forward_1_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_forward_1_coeff[k]/ddy2});
    }
    
    // Middle rows
    for(i0=1;i0<this->grid->Ny-1;i0++){
      j0 = 0;
      for(int k=0;k<RectGrid::derivative_2_forward_1_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_forward_1_index[k]),RectGrid::derivative_2_forward_1_coeff[k]/ddx2});
      }
      for(int k=0;k<RectGrid::derivative_2_central_2_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_central_2_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_central_2_coeff[k]/ddy2});
      }
      
      for(int j0=1;j0<this->grid->Nx-1;j0++){
	for(int k=0;k<RectGrid::derivative_2_central_2_index.size();k++){
	  tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_central_2_index[k]),RectGrid::derivative_2_central_2_coeff[k]/ddx2});
	  tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_central_2_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_central_2_coeff[k]/ddy2});
	}
      }
      
      j0 = this->grid->Nx-1;
      for(int k=0;k<RectGrid::derivative_2_backward_1_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_backward_1_index[k]),RectGrid::derivative_2_backward_1_coeff[k]/ddx2});
      }
      for(int k=0;k<RectGrid::derivative_2_central_2_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_central_2_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_central_2_coeff[k]/ddy2});
      } 
    }

    // Last Row
    i0 = this->grid->Ny-1;
    j0 = 0;
    for(int k=0;k<RectGrid::derivative_2_forward_1_index.size();k++){
      tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_forward_1_index[k]),RectGrid::derivative_2_forward_1_coeff[k]/ddx2});
    }
    for(int k=0;k<RectGrid::derivative_2_backward_1_index.size();k++){
      tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_backward_1_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_backward_1_coeff[k]/ddy2});
    }

    for(j0=1;j0<this->grid->Nx-1;j0++){
      for(int k=0;k<RectGrid::derivative_2_central_2_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_central_2_index[k]),RectGrid::derivative_2_central_2_coeff[k]/ddx2});
      }
      for(int k=0;k<RectGrid::derivative_2_backward_1_index.size();k++){
	tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_backward_1_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_backward_1_coeff[k]/ddy2});
      }
    }

    j0 = this->grid->Nx-1;
    for(int k=0;k<RectGrid::derivative_2_backward_1_index.size();k++){
      tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_backward_1_index[k]),RectGrid::derivative_2_backward_1_coeff[k]/ddx2});
      tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_backward_1_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_backward_1_coeff[k]/ddy2});
    }

  } else if( this->reg == "curvature" ){//-------------------> curvature with units at the borders matrix

    this->eigenSparseMemoryAllocForH = 8;
    int i0,j0;
    double ddx2 = pow(this->grid->step_x,2);
    double ddy2 = pow(this->grid->step_y,2);
    
    // First Row
    i0 = 0;
    for(j0=0;j0<this->grid->Nx;j0++){
      tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+j0,1.0});
    }
    
    // Middle rows
    for(i0=1;i0<this->grid->Ny-1;i0++){
      tmp.push_back({i0*this->grid->Nx,i0*this->grid->Nx,1.0});
      
      for(int j0=1;j0<this->grid->Nx-1;j0++){
	for(int k=0;k<RectGrid::derivative_2_central_2_index.size();k++){
	  tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+(j0+RectGrid::derivative_2_central_2_index[k]),RectGrid::derivative_2_central_2_coeff[k]/ddx2});
	  tmp.push_back({i0*this->grid->Nx+j0,(i0+RectGrid::derivative_2_central_2_index[k])*this->grid->Nx+j0,RectGrid::derivative_2_central_2_coeff[k]/ddy2});
	}
      }
      
      tmp.push_back({i0*this->grid->Nx+Nx-1,i0*this->grid->Nx+Nx-1,1.0});
    }

    // Last Row
    i0 = this->grid->Ny-1;
    for(j0=0;j0<this->grid->Nx;j0++){
      tmp.push_back({i0*this->grid->Nx+j0,i0*this->grid->Nx+j0,1.0});
    }
    
  } else if( this->reg == "covariance_kernel" ){//-------------------> covariance matrix
    
    int* nonZeroRow = (int*) calloc(this->Sm,sizeof(int));
    int index1,index2;
    double x1,y1,x2,y2;
    double cov,r;
    // Match each grid point to every other point through this double loop
    for(int i1=0;i1<this->grid->Ny-1;i1++){
      for(int j1=0;j1<this->grid->Nx-1;j1++){
	x1 = this->grid_center_x[j1];
	y1 = this->grid_center_y[i1];
	index1 = i1*this->grid->Nx+j1;

	cov = this->kernel->getCovarianceSelf();
	if( cov > this->kernel->cmax ){
	  tmp.push_back({index1,index1,cov});
	  nonZeroRow[index1]++;
	}	

	for(int i2=i1+1;i2<this->grid->Ny;i2++){
	  for(int j2=j1+1;j2<this->grid->Nx;j2++){
	    x2 = this->grid_center_x[j2];
	    y2 = this->grid_center_y[i2];
	    index2 = i2*this->grid->Nx+j2;
	    
	    r = hypot(x1-x2,y1-y2);
	    cov = this->kernel->getCovariance(r);
	    if( cov > this->kernel->cmax ){
	      tmp.push_back({index1,index2,cov});
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

  this->H.tri.swap(tmp);
}

//virtual
void FixedSource::outputSource(const std::string fname){
  std::vector<std::string> key{"WIDTH","HEIGHT"};
  std::vector<std::string> val{std::to_string(this->grid->width),std::to_string(this->grid->height)};
  std::vector<std::string> txt{"width of the source in arcsec","height of the source in arcsec"};
  FitsInterface::writeFits(this->Si,this->Sj,this->grid->z,key,val,txt,fname);
}
