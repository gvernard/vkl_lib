#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

using namespace vkl;

//Derived class from BaseSourcePlane: AdaptiveSource
//===============================================================================================================
AdaptiveSource::AdaptiveSource(std::string mode,int Sm,double* x,double* y,std::string reg_scheme){
  type = "adaptive";
  mode = mode;
  Sm   = Sm;

  reg = reg_scheme;
  if( reg == "identity"){
    eigenSparseMemoryAllocForH = 1;
  } else if( reg == "gradient" ){
    eigenSparseMemoryAllocForH = 8;
  } else if( reg == "curvature" ){
    eigenSparseMemoryAllocForH = 8;
  } else if( reg == "covariance_kernel" || reg == "covariance_kernel_in_identity_out" ){
    // I need to call constructH() in order to set the number of non-zero entries per sparse matrix row (this number varies for a random covariance kernel).
    // This happens at the initialization of the likelihood model (BaseLikelihoodModel->initializeAlgebra()) just before setting the algebra.
    // So, do nothing here.
  }
}

AdaptiveSource::AdaptiveSource(const AdaptiveSource& other) : BaseSourcePlane(other) {
  this->type = other.type;
  this->Sm = other.Sm;

  for(int i=0;i<this->Sm;i++){
    this->src[i]  = other.src[i];
    this->x[i]    = other.x[i];
    this->y[i]    = other.y[i];
    this->s_dx[i] = other.s_dx[i];
    this->s_dy[i] = other.s_dy[i];
    this->mask_vertices[i] = other.mask_vertices[i];
    this->lambda_out[i]    = other.lambda_out[i];
  }
  
  this->reg = other.reg;
  this->eigenSparseMemoryAllocForH = other.eigenSparseMemoryAllocForH;
  this->sample_reg = other.sample_reg;

  this->n_triangles = other.n_triangles;
  this->triangles.resize( other.triangles.size() );
  for(int i=0;i<this->triangles.size();i++){
    this->triangles[i].a = other.triangles[i].a;
    this->triangles[i].b = other.triangles[i].b;
    this->triangles[i].c = other.triangles[i].c;
  }
  
  this->opposite_edges_per_vertex.resize( other.opposite_edges_per_vertex.size() );
  for(int i=0;i<this->opposite_edges_per_vertex.size();i++){
    copy(other.opposite_edges_per_vertex[i].begin(),other.opposite_edges_per_vertex[i].end(),back_inserter(this->opposite_edges_per_vertex[i])); 
  }

  if( this->reg == "covariance_kernel"  || this->reg == "covariance_kernel_in_identity_out" ){
    this->kernel = other.kernel->clone();
  }
}

AdaptiveSource::AdaptiveSource* clone(){
  return new AdaptiveSource(*this);
};

AdaptiveSource::~AdaptiveSource(){
  std::vector<a_triangle> ().swap(triangles);
}

//non-virtual
void AdaptiveSource::createAdaGrid(ImagePlane* image,CollectionMassModels* mycollection){
  // Create the adaptive grid

  if( this->mode == "random" ){
    
    double xtmp,ytmp;
    srand48(time(NULL));
    for(int i=0;i<this->Sm;i++){
      xtmp = drand48()*(image->xmax - image->xmin) + image->xmin;
      ytmp = drand48()*(image->ymax - image->ymin) + image->ymin;
      mycollection->all_defl(xtmp,ytmp,this->x[i],this->y[i]);   
    }

    for(int i=0;i<image->Nm;i++){
      image->active[i] = -1;
    }
    this->inMask(image);

  } else if( this->mode == "image" ){

    for(int i=0;i<image->Nm;i++){
      image->active[i] = -1;
    }

    int i0    = (int) floor( (this->spacing-1)/2.0 );
    int j0    = (int) floor( (this->spacing-1)/2.0 );
    int count = 0;//must go up to Sm
    for(int i=i0;i<image->Ni;i=i+this->spacing){
      for(int j=j0;j<image->Nj;j=j+this->spacing){
	this->x[count] = image->defl_x[i*image->Nj+j];
	this->y[count] = image->defl_y[i*image->Nj+j];
	image->active[i*image->Nj+j] = count;
	count++;
      }
    }

  } else if( this->mode == "grid" ){

    int Nj = image->Nj;
    int Ni = image->Ni;

    int i0    = floor(Ni/2);
    int j0    = floor(Nj/2);
    double di = image->width/(Ni);
    double dj = image->height/(Nj);
    
    int count = 0;
    double xtmp,ytmp;
    //    FILE* fh1 = fopen("ada_grid.dat","w");
    for(int ii=0;ii<Ni;ii=ii+this->spacing){
      for(int jj=0;jj<Nj;jj=jj+this->spacing){
	xtmp   =  (jj-j0)*di;
	ytmp   = -(ii-i0)*dj;//reflect y-axis
	mycollection->all_defl(xtmp,ytmp,this->x[count],this->y[count]);
	//	fprintf(fh1,"%12.4f %12.4f\n",this->x[count],this->y[count]);
	count++;
      }
    }
    //    fclose(fh1);
    //    std::cout << "The number of pixels in the adaptive source is " << count << " and it must be " << this->Sm << std::endl;


    std::vector<xypoint> frame(8);
    double d_out = 5.0;
    double d_in  = 0.000001;

    frame[0].x = image->x[0]         - d_out;
    frame[0].y = image->y[0]         + d_out;
    frame[1].x = image->x[Nj-1]      + d_out;
    frame[1].y = image->y[0]         + d_out;
    frame[2].x = image->x[0]         - d_out;
    frame[2].y = image->y[(Ni-1)*Nj] - d_out;
    frame[3].x = image->x[Nj-1]      + d_out;
    frame[3].y = image->y[(Ni-1)*Nj] - d_out;

    frame[4].x = mycollection->models[0]->mpars["x0"] - d_in;
    frame[4].y = mycollection->models[0]->mpars["y0"] + d_in;
    frame[5].x = mycollection->models[0]->mpars["x0"] + d_in;
    frame[5].y = mycollection->models[0]->mpars["y0"] + d_in;
    frame[6].x = mycollection->models[0]->mpars["x0"] - d_in;
    frame[6].y = mycollection->models[0]->mpars["y0"] - d_in;
    frame[7].x = mycollection->models[0]->mpars["x0"] + d_in;
    frame[7].y = mycollection->models[0]->mpars["y0"] - d_in;


    //    FILE* fh = fopen("outer_grid.dat","w");
    for(int i=0;i<frame.size();i++){
      mycollection->all_defl(frame[i].x,frame[i].y,this->x[count+i],this->y[count+i]);
      //      fprintf(fh,"%12.4f %12.4f\n",this->x[count+i],this->y[count+i]);
    }
    //    fclose(fh);

    for(int i=0;i<image->Nm;i++){
      image->active[i] = -1;
    }
    this->inMask(image);

  }

}

//non-virtual
void AdaptiveSource::inMask(ImagePlane* image){
  // Find which source pixels are within the data mask

  if( this->mode == "random" ){
    
  } else if( this->mode == "image" ){

    int i0    = (int) floor( (this->spacing-1)/2.0 );
    int j0    = (int) floor( (this->spacing-1)/2.0 );
    int count = 0;//must go up to Sm
    int mask_count = 0;
    this->in_total.clear();
    this->total_in.clear();
    for(int i=i0;i<image->Ni;i=i+this->spacing){
      for(int j=j0;j<image->Nj;j=j+this->spacing){
	if( image->S.tri[i*image->Nj+j].v == 1 ){
	  this->mask_vertices[count] = 1;
	  this->lambda_out[count] = 0.0;
	  this->in_total[mask_count] = count;
	  this->total_in[count] = mask_count;
	  mask_count++;
	} else {
	  this->mask_vertices[count] = 0;
	  int sep = 1;
	  double sum = 0.0;
	  int max = 2*sep+1;
	  for(int ii=0;ii<max;ii++){
	    for(int jj=0;jj<max;jj++){
	      sum += sqrt(1.0/image->C.tri[(i-sep+ii)*image->Nj+(j-sep+jj)].v); // the covariance matrix C holds 1/noise^2
	    }
	  }
	  this->lambda_out[count] = (double) max*max/sum;
	  //	  std::cout << this->lambda_out[count] << " " << sqrt(image->C.tri[i*image->Nj+j].v) << std::endl;
	  //	  this->lambda_out[count] = 1.09;
	}
	count++;
      }
    }
    this->Smask = mask_count;

    this->lambda_out_sum = 0.0;
    for(int i=0;i<this->Sm;i++){
      if( this->lambda_out[i] != 0 ){
	this->lambda_out_sum += log(this->lambda_out[i]);
      }
    }

  } else if( this->mode == "grid" ){

  }

}

//non-virtual
void AdaptiveSource::createDelaunay(){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel            K;
  typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int,K>    Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<unsigned int,K>      Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    Tds;
  typedef CGAL::Delaunay_triangulation_2<K,Tds>                          Delaunay;
  typedef K::Point_2                                                     Point;

  std::vector< std::pair<Point,int> > points;
  for(int i=0;i<this->Sm;i++){
    points.push_back( std::make_pair(Point(this->x[i],this->y[i]),i) );
  }

  Delaunay triangulation;
  triangulation.insert(points.begin(),points.end());


  //Get each Delaunay triangle in my own struct, and give them an ID
  //[constructing this->triangles]
  Delaunay::Finite_faces_iterator fit;
  Delaunay::Face_handle face;
  a_triangle triangle;
  int face_id = 0;
  this->n_triangles = triangulation.number_of_faces();
  this->triangles.resize( triangulation.number_of_faces() );

  for(fit=triangulation.finite_faces_begin();fit!=triangulation.finite_faces_end();fit++){
    face = fit;
    face->info() = face_id;

    triangle.a = (int) face->vertex(0)->info();
    triangle.b = (int) face->vertex(1)->info();
    triangle.c = (int) face->vertex(2)->info();
    
    this->triangles[face_id] = triangle;
    face_id++;

    /*    
    double area = this->triangleArea(triangle);
    if( area < 0.01 ){
      nsmall++;
      std::cout << area << std::endl;
    } else {
      std::cout << "normal" << std::endl;
    }
    */
  }


  //Find which triangle IDs around each vertex. This is zero for the convex-hull vertices (having 1 or more infinite incident faces) 
  //[constructing this->face_ids_per_vertex]
  Delaunay::Face_circulator fc;
  Delaunay::Vertex_circulator vc;
  Delaunay::Finite_vertices_iterator fvi;
  this->opposite_edges_per_vertex.resize( this->Sm );
  std::vector<int> opposite_indices;
  for(fvi=triangulation.finite_vertices_begin();fvi!=triangulation.finite_vertices_end();fvi++){

    int inf = 0;
    fc = triangulation.incident_faces(fvi);
    do{
      //      vertex_face_ids.push_back( fc->info() );
      vc = triangulation.incident_vertices(fvi,fc);
      opposite_indices.push_back( vc->info() );
      vc++;
      opposite_indices.push_back( vc->info() );      
      if( triangulation.is_infinite(fc) ){
	inf++;
      }
    }while( ++fc != triangulation.incident_faces(fvi) );

    if( inf > 0 ){
      opposite_indices.resize(0);
    }

    this->opposite_edges_per_vertex[ fvi->info() ] = opposite_indices;
    opposite_indices.resize(0);
  }

}

//non-virtual
double AdaptiveSource::triangleArea(a_triangle triangle){
  double xa = this->x[triangle.a];
  double ya = this->y[triangle.a];
  double xb = this->x[triangle.b];
  double yb = this->y[triangle.b];
  double xc = this->x[triangle.c];
  double yc = this->y[triangle.c];

  double ax = (xb - xa);
  double ay = (yb - ya);
  double bx = (xc - xa);
  double by = (yc - ya);

  return (ax*by - ay*bx)/2.0;
}

//non-virtual
void AdaptiveSource::createInterpolationWeights(ImagePlane* image){
  double wa,wb,wc;
  double ybc,xac,xcb,yac,xxc,yyc,den;
  a_triangle triangle;
  int flag = 1;
  
  for(int i=0;i<image->Nm;i++){
    if( image->active[i] == -1 ){

      flag = 0;
      for(int j=0;j<this->n_triangles;j++){
	triangle = this->triangles[j];
	
	ybc = this->y[triangle.b] - this->y[triangle.c];//(yb-yc)
	xac = this->x[triangle.a] - this->x[triangle.c];//(xa-xc)
	xcb = this->x[triangle.c] - this->x[triangle.b];//(xc-xb)
	yac = this->y[triangle.a] - this->y[triangle.c];//(ya-yc)
	xxc = image->defl_x[i]    - this->x[triangle.c];//(x -xc)
	yyc = image->defl_y[i]    - this->y[triangle.c];//(y -yc)
	den = ybc*xac + xcb*yac;
	
	wa = ( ybc*xxc+xcb*yyc)/den;
	wb = (-yac*xxc+xac*yyc)/den;
	wc = 1.0 - wa - wb;
	
	if( 0.0 <= wa && wa <= 1.0 && 0.0 <= wb && wb <= 1.0 && 0.0 <= wc && wc <= 1.0 ){
	  flag = 1;
	  delete(image->cells[i]);
	  InterpolationCell* cell = new InterpolationCell(3);
	  cell->ind[0] = triangle.a;
	  cell->ind[1] = triangle.b;
	  cell->ind[2] = triangle.c;
	  cell->wei[0] = wa;
	  cell->wei[1] = wb;
	  cell->wei[2] = wc;
	  image->cells[i] = cell;
	  break;
	}
      }
      
      if( flag == 0 ){
	delete(image->cells[i]);
	InterpolationCell* cell = new InterpolationCell(1);
	cell->ind[0] = 0;
	cell->wei[0] = 0.0;
	image->cells[i] = cell;
      }

    } else {

      delete(image->cells[i]);
      InterpolationCell* cell = new InterpolationCell(1);
      cell->ind[0] = image->active[i];
      cell->wei[0] = 1.0;
      image->cells[i] = cell;

    }
  }

}



//virtual
void AdaptiveSource::constructH(){
  std::vector<mytriplet> tmp;//need to make sure that the L triplet vector is a new one

  if( this->reg == "identity" ){//---------------------------> zero order
    
    for(int i=0;i<this->H.Ti;i++){
      tmp.push_back({i,i,1});
    }

  } else if ( this->reg == "gradient" || this->reg == "curvature" || this->reg == "curvature_in_identity_out" ){

    for(int i=0;i<this->Sm;i++){

      xypoint p0 = {this->x[i],this->y[i]};
      std::vector<int> indices = this->opposite_edges_per_vertex[i];
      
      if( indices.size() == 0 ){

	//we are on a convex-hull point that has no smoothing
	tmp.push_back({i,i,1.0});

      } else {
	
	std::map<int,double> weights;

	for(int j=0;j<indices.size();j=j+2){
	  double ymin = 0.0;
	  double ymax = 0.0;
	  double xmin = 0.0;
	  double xmax = 0.0;
	  
	  xypoint p1 = {this->x[indices[j]],this->y[indices[j]]};
	  xypoint p2 = {this->x[indices[j+1]],this->y[indices[j+1]]};
	  
	  //Smoothing in the x-direction
	  if( p1.y > p2.y ){
	    ymax = p1.y;
	    ymin = p2.y;
	  } else {
	    ymax = p2.y;
	    ymin = p1.y;
	  }
	  
	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( ymin <= p0.y && p0.y < ymax ){
	    xypoint pint = intersection_point_x(p0,p1,p2);
	    double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	    double l1 = 1.0 - l2;
	    //	    double d12 = pow(p1.x-p2.x,2) + pow(p1.y-p2.y,2);
	    //	    double di2 = pow(pint.x-p2.x,2) + pow(pint.y-p2.y,2);
	    //	    double l1 = sqrt(di2/d12);
	    //	    double l2 = 1.0 - l1;
	    double d  = fabs(p0.x-pint.x);
	    weights[ i            ] -= 1.0/d;
	    weights[ indices[j]   ] += l1/d;
	    weights[ indices[j+1] ] += l2/d;
	  }

	  //Smoothing in the y-direction
	  if( p1.x > p2.x ){
	    xmax = p1.x;
	    xmin = p2.x;
	  } else {
	    xmax = p2.x;
	    xmin = p1.x;
	  }

	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( xmin <= p0.x && p0.x < xmax ){
	    xypoint pint = intersection_point_y(p0,p1,p2);
	    double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	    double l1 = 1.0 - l2;
	    //	    double d12 = pow(p1.x-p2.x,2) + pow(p1.y-p2.y,2);
	    //	    double di2 = pow(pint.x-p2.x,2) + pow(pint.y-p2.y,2);
	    //	    double l1 = sqrt(di2/d12);
	    //	    double l2 = 1.0 - l1;
	    double d  = fabs(p0.y-pint.y);
	    weights[ i            ] -= 1.0/d;
	    weights[ indices[j]   ] += l1/d;
	    weights[ indices[j+1] ] += l2/d;
	  }

	}
	
	//Loop through weights and add entries to H
	for(std::map<int,double>::iterator it2=weights.begin();it2!=weights.end();it2++){
	  tmp.push_back({i,it2->first,it2->second});
	  //	  std::cout << i << " " << iterator->first << " " << iterator->second << std::endl;
	}
	//	std::cout << i << " " << weights.size() << " " << indices.size()/2.0;
	//	for(it_int_double iterator=weights.begin();iterator!=weights.end();iterator++){
	//	  printf("%20.5e",iterator->second);
	//	}
	//for(int k=0;k<indices.size();k++){
	//  std::cout << indices[k] << " ";
	//}
	//	std::cout << std::endl;
	weights.clear();
      }
      //      break;
    }
    
    /*
  } else if ( this->reg == "curvature_in_identity_out" ){

    this->H.Ti = this->Smask;
    this->H.Tj = this->Smask;
  
    for(std::map<int,int>::iterator it=this->in_total.begin();it!=in_total.end();it++){
      int ind_in = it->first;
      int ind_total = it->second;

      xypoint p0 = {this->x[ind_total],this->y[ind_total]};
      std::vector<int> indices = this->opposite_edges_per_vertex[ind_total];

      // check if the neighbours of the target pixel are also in the mask
      int flag_out_of_mask = 0;
      for(int j=0;j<indices.size();j=j+2){
	if( this->mask_vertices[indices[j]] == 0 ){
	  flag_out_of_mask = 1;
	  break;
	}
      }

      if( flag_out_of_mask == 1 ){
	
	// this pixel has neighbours that are outside the mask
	tmp.push_back({ind_in,ind_in,1.0});

      } else {

	std::map<int,double> weights;
	
	for(int j=0;j<indices.size();j=j+2){
	  double ymin = 0.0;
	  double ymax = 0.0;
	  double xmin = 0.0;
	  double xmax = 0.0;
	  
	  xypoint p1 = {this->x[indices[j]],this->y[indices[j]]};
	  xypoint p2 = {this->x[indices[j+1]],this->y[indices[j+1]]};
	  
	  //Smoothing in the x-direction
	  if( p1.y > p2.y ){
	    ymax = p1.y;
	    ymin = p2.y;
	  } else {
	    ymax = p2.y;
	    ymin = p1.y;
	  }
	  
	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( ymin <= p0.y && p0.y < ymax ){
	    xypoint pint = intersection_point_x(p0,p1,p2);
	    double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	    double l1 = 1.0 - l2;
	    double d  = fabs(p0.x-pint.x);
	    weights[ ind_in                 ] -= 1.0/d;
	    weights[ total_in[indices[j]]   ] += l1/d;
	    weights[ total_in[indices[j+1]] ] += l2/d;
	  }

	  //Smoothing in the y-direction
	  if( p1.x > p2.x ){
	    xmax = p1.x;
	    xmin = p2.x;
	  } else {
	    xmax = p2.x;
	    xmin = p1.x;
	  }

	  //Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	  if( xmin <= p0.x && p0.x < xmax ){
	    xypoint pint = intersection_point_y(p0,p1,p2);
	    double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	    double l1 = 1.0 - l2;
	    double d  = fabs(p0.y-pint.y);
	    weights[ ind_in                 ] -= 1.0/d;
	    weights[ total_in[indices[j]]   ] += l1/d;
	    weights[ total_in[indices[j+1]] ] += l2/d;
	  }

	}
	
	//Loop through weights and add entries to H
	for(std::map<int,double>::iterator it2=weights.begin();it2!=weights.end();it2++){
	  tmp.push_back({ind_in,it2->first,it2->second});
	}
	weights.clear();

      }

    }
    */

  } else if ( this->reg == "covariance_kernel" || this->reg == "covariance_kernel_in_identity_out" ){//-------------------> covariance matrix

    int* nonZeroRow = (int*) calloc(this->Sm,sizeof(int));
    double cov,r;
    for(int i=0;i<this->Sm;i++){
      for(int j=0;j<this->Sm;j++){
	if( i != j ){
	  r = hypot(this->x[j]-this->x[i],this->y[j]-this->y[i]);
	  cov = this->kernel->getCovariance(r);
	} else {
	  cov = this->kernel->getCovarianceSelf();
	}
	if( cov > this->kernel->cmax ){
	  tmp.push_back({i,j,cov});
	  nonZeroRow[i]++;
	}
      }
    }

    int maxNonZero = nonZeroRow[0];
    for(int i=1;i<this->Sm;i++){
      if( nonZeroRow[i] > maxNonZero ){
	maxNonZero = nonZeroRow[i];
      }
    }
    free(nonZeroRow);
    this->eigenSparseMemoryAllocForH = maxNonZero;

  }

    /*
  } else if ( this->reg == "covariance_kernel_in_identity_out" ){

    int* nonZeroRow = (int*) calloc(this->Sm,sizeof(int));
    for(int i=0;i<this->Sm;i++){
      if( this->mask_vertices[i] == 0 ){
	
	tmp.push_back({i,i,1.0});
	
      } else {

	double cov,r;
	for(int j=0;j<this->Sm;j++){
	  if( i != j ){
	    r = hypot(this->x[j]-this->x[i],this->y[j]-this->y[i]);
	    cov = this->kernel->getCovariance(r);
	  } else {
	    cov = this->kernel->getCovarianceSelf();
	  }
	  if( cov > this->kernel->cmax ){
	    tmp.push_back({i,j,cov});
	    nonZeroRow[i]++;
	  }
	}

      }
    }

    int maxNonZero = nonZeroRow[0];
    for(int i=1;i<this->Sm;i++){
      if( nonZeroRow[i] > maxNonZero ){
	maxNonZero = nonZeroRow[i];
      }
    }
    free(nonZeroRow);
    this->eigenSparseMemoryAllocForH = maxNonZero;

  }
    */

  this->H.tri.swap(tmp);
}



//virtual
void AdaptiveSource::constructDerivatives(){
  // Calculate the derivative at every source grid point
  double* dev_x_val   = (double*) malloc(2*sizeof(double));
  double* dev_x_coord = (double*) malloc(2*sizeof(double));
  double* dev_y_val   = (double*) malloc(2*sizeof(double));
  double* dev_y_coord = (double*) malloc(2*sizeof(double));
  for(int i=0;i<this->Sm;i++){
    xypoint p0 = {this->x[i],this->y[i]};
    std::vector<int> indices = this->opposite_edges_per_vertex[i];

    if( indices.size() == 0 ){
      //we are on a convex-hull point that has zero source derivative
      this->s_dx[i] = 0.0;
      this->s_dy[i] = 0.0;

    } else {
      double ymin,ymax,xmin,xmax;
      int i_x = 0;
      int i_y = 0;

      for(int j=0;j<indices.size();j=j+2){	
	xypoint p1 = {this->x[indices[j]],this->y[indices[j]]};
	xypoint p2 = {this->x[indices[j+1]],this->y[indices[j+1]]};
	

	//Smoothing in the x-direction
	if( p1.y > p2.y ){
	  ymax = p1.y;
	  ymin = p2.y;
	} else {
	  ymax = p2.y;
	  ymin = p1.y;
	}
	//Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weights
	if( ymin <= p0.y && p0.y < ymax ){
	  xypoint pint = intersection_point_x(p0,p1,p2);
	  double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	  double l1 = 1.0 - l2;
	  dev_x_val[i_x]   = l1*this->src[indices[j]] + l2*this->src[indices[j+1]];
	  dev_x_coord[i_x] = pint.x - p0.x;
	  i_x++;
	}

	//Smoothing in the y-direction
	if( p1.x > p2.x ){
	  xmax = p1.x;
	  xmin = p2.x;
	} else {
	  xmax = p2.x;
	  xmin = p1.x;
	}	
	//Find intersection point (if it exists), interpolate to get the weights, and add them to <map>weight
	if( xmin <= p0.x && p0.x < xmax ){
	  xypoint pint = intersection_point_y(p0,p1,p2);
	  double l2 = ((p1.y-p0.y)*(pint.x-p1.x)+(p0.x-p1.x)*(pint.y-p1.y))/((p2.y-p1.y)*(p0.x-p1.x)+(p1.x-p2.x)*(p0.y-p1.y));
	  double l1 = 1.0 - l2;
	  dev_y_val[i_y]   = l1*this->src[indices[j]] + l2*this->src[indices[j+1]];
	  dev_y_coord[i_y] = pint.y - p0.y;
	  i_y++;
	}
	
      }

      if( dev_x_coord[0] > dev_x_coord[1] ){
	this->s_dx[i] = (dev_x_val[0] - dev_x_val[1])/(dev_x_coord[0] - dev_x_coord[1]);
      } else {
	this->s_dx[i] = (dev_x_val[1] - dev_x_val[0])/(dev_x_coord[1] - dev_x_coord[0]);
      }

      if( dev_y_coord[0] > dev_y_coord[1] ){
	this->s_dy[i] = (dev_y_val[0] - dev_y_val[1])/(dev_y_coord[0] - dev_y_coord[1]);
      } else {
	this->s_dy[i] = (dev_y_val[1] - dev_y_val[0])/(dev_y_coord[1] - dev_y_coord[0]);
      }

    }
    
  }

  free(dev_x_coord);
  free(dev_y_coord);
  free(dev_x_val);
  free(dev_y_val);
}


//virtual
void AdaptiveSource::constructDs(ImagePlane* image){
  // Deflect the image grid, find the triangle that each ray belongs to, and interpolate between the derivatives of the vertices
  // or get the derivative of the vertex if the ray is part of the adaptive grid.
  this->Ds.Ti = image->Nm;
  this->Ds.Tj = 2*image->Nm;
  std::vector<mytriplet> tmp;

  for(int i=0;i<image->Nm;i++){
    double dev_x = 0.0;
    double dev_y = 0.0;
    for(int j=0;j<image->cells[i]->size;j++){
      dev_x += image->cells[i]->wei[j]*this->s_dx[image->cells[i]->ind[j]];
      dev_y += image->cells[i]->wei[j]*this->s_dy[image->cells[i]->ind[j]];
    }
    tmp.push_back({i,2*i,dev_x});
    tmp.push_back({i,2*i+1,dev_y});
  }

  this->Ds.tri.swap(tmp);
}






//private
AdaptiveSource::xypoint AdaptiveSource::intersection_point_x(xypoint p0,xypoint p1,xypoint p2){
  xypoint pint;

  if( p1.x == p2.x ){
    pint.x = p1.x;
    pint.y = p0.y;
  } else {
    pint.x = p1.x + (p0.y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y);
    pint.y = p0.y;
  }

  return pint;
}

AdaptiveSource::xypoint AdaptiveSource::intersection_point_y(xypoint p0,xypoint p1,xypoint p2){
  xypoint pint;

  if( p1.y == p2.y ){
    pint.x = p0.x;
    pint.y = p1.y;
  } else {
    pint.x = p0.x;
    pint.y = (p2.y-p1.y)*(p0.x-p1.x)/(p2.x-p1.x) + p1.y;
  }

  return pint;
}


//non-virtual
void AdaptiveSource::writeTriangles(){
  a_triangle triangle;
  FILE* fh = fopen("triangles.dat","w");
  for(int i=0;i<this->n_triangles;i++){
    triangle = this->triangles[i];
    fprintf(fh,"%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",this->x[triangle.a],this->y[triangle.a],this->x[triangle.b],this->y[triangle.b],this->x[triangle.c],this->y[triangle.c]);
  }
  fclose(fh);
}

//non-virtual
void AdaptiveSource::writeVertexFaces(int index,std::string output){
  std::string fname = output + std::to_string(index) + "_vertex_faces.dat";
  FILE* fh = fopen(fname.c_str(),"w");
  int b,c;
  for(int j=0;j<this->opposite_edges_per_vertex[index].size();j=j+2){
    b = this->opposite_edges_per_vertex[index][j];
    c = this->opposite_edges_per_vertex[index][j+1];
    fprintf(fh,"%20.5f %20.5f %20.5f %20.5f %20.5f %20.5f\n",this->x[index],this->y[index],this->x[b],this->y[b],this->x[c],this->y[c]);
  }
  fclose(fh);
}

//virtual
void AdaptiveSource::outputSource(const std::string path){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
  typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
  typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
  typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    Voronoi;
  typedef Voronoi::Locate_result                                               Locate_result;
  typedef Voronoi::Face_handle                                                 Face_handle;
  typedef Voronoi::Ccb_halfedge_circulator                                     Ccb_halfedge_circulator;
  typedef K::Point_2                                                           Point;

  //Calculate the Voronoi graph
  Voronoi voronoi;
  std::vector<Point> points;
  for(int i=0;i<this->Sm;i++){
    points.push_back( Point(this->x[i],this->y[i]) );
  }
  voronoi.insert(points.begin(),points.end());

  //Locate the centres of the Voronoi cells and then get the vertices of the surrounding face
  //This code is based on an online example from the CGAL voronoi package
  std::string filename = path + "source_voronoi.dat";
  FILE* fh = fopen(filename.c_str(),"w");
  for(int i=0;i<this->Sm;i++){

    //    std::cout << "111" << std::endl;

    int mask_flag = 0;
    for(int k=0;k<this->opposite_edges_per_vertex[i].size();k++){
      int vertex_id = this->opposite_edges_per_vertex[i][k];
      if( this->mask_vertices[vertex_id] == 1 ){
	mask_flag = 1;
	break;
      }
    }

    //    std::cout << "222" << std::endl;
   
    if( mask_flag == 1 ){
      //    if( this->opposite_edges_per_vertex[i].size() != 0 ){

      fprintf(fh,"%20.5f",src[i]);

      Locate_result f   = voronoi.locate(Point(this->x[i],this->y[i]));
      Face_handle* face = boost::get<Face_handle>(&f);
    
      Ccb_halfedge_circulator ec_start = (*face)->ccb();
      Ccb_halfedge_circulator ec = ec_start;
      do {
	Point p = ec->source()->point();
	fprintf(fh,"%20.5f %20.5f",p.x(),p.y());
      } while( ++ec != ec_start );
      fprintf(fh,"\n");
    }

    //    std::cout << "333" << std::endl;
	  

  }
  fclose(fh);

  // Output the source vertices and values
  std::string filename2 = path + "source_irregular.dat";
  FILE* fh2 = fopen(filename2.c_str(),"w");
  for(int i=0;i<this->Sm;i++){
    //    if( this->mask_vertices[i] == 1 ){
      fprintf(fh2,"%20.5f %20.5f %20.5f\n",this->src[i],this->x[i],this->y[i]);
      //    }
  }  
  fclose(fh2);
}


//virtual
void AdaptiveSource::outputSourceErrors(double* errors,const std::string path){
  typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
  typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
  typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
  typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    Voronoi;
  typedef Voronoi::Locate_result                                               Locate_result;
  typedef Voronoi::Face_handle                                                 Face_handle;
  typedef Voronoi::Ccb_halfedge_circulator                                     Ccb_halfedge_circulator;
  typedef K::Point_2                                                           Point;

  //Calculate the Voronoi graph
  Voronoi voronoi;
  std::vector<Point> points;
  for(int i=0;i<this->Sm;i++){
    points.push_back( Point(this->x[i],this->y[i]) );
  }
  voronoi.insert(points.begin(),points.end());

  //Locate the centres of the Voronoi cells and then get the vertices of the surrounding face
  //This code is based on an online example from the CGAL voronoi package
  FILE* fh = fopen((path+"source_voronoi_errors.dat").c_str(),"w");
  for(int i=0;i<this->Sm;i++){

    int mask_flag = 0;
    for(int k=0;k<this->opposite_edges_per_vertex[i].size();k++){
      int vertex_id = this->opposite_edges_per_vertex[i][k];
      if( this->mask_vertices[vertex_id] == 1 ){
	mask_flag = 1;
	break;
      }
    }
    
    if( mask_flag == 1 ){
      //   if( this->opposite_edges_per_vertex[i].size() != 0 ){

      fprintf(fh,"%20.5f",errors[i]);

      Locate_result f   = voronoi.locate(Point(this->x[i],this->y[i]));
      Face_handle* face = boost::get<Face_handle>(&f);
    
      Ccb_halfedge_circulator ec_start = (*face)->ccb();
      Ccb_halfedge_circulator ec = ec_start;
      do {
	Point p = ec->source()->point();
	fprintf(fh,"%20.5f %20.5f",p.x(),p.y());
      } while( ++ec != ec_start );
      fprintf(fh,"\n");
    }

  }
  fclose(fh);
}
