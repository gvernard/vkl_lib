#include "irreGrid.hpp"

using namespace vkl;

// Class: Triangle
//===============================================================================================================
Triangle::Triangle(double* x,double* y){
  for(int i=0;i<3;i++){
    this->x[i] = x[i];
    this->y[i] = y[i];
  }
}

Triangle::Triangle(double* x,double* y,double* z){
  for(int i=0;i<3;i++){
    this->x[i] = x[i];
    this->y[i] = y[i];
    this->z[i] = z[i];
  }
}

double Triangle::triangle_area(){
  double ax = (this->x[1] - this->x[0]);
  double ay = (this->y[1] - this->y[0]);
  double bx = (this->x[2] - this->x[0]);
  double by = (this->y[2] - this->y[0]);
  return (ax*by - ay*bx)/2.0;
}

bool Triangle::interp_bary(double x,double y,double& value){
  double wa,wb,wc;
  double ybc,xac,xcb,yac,xxc,yyc,den;
  ybc = this->y[1] - this->y[2]; // (yb-yc)
  xac = this->x[0] - this->x[2]; // (xa-xc)
  xcb = this->x[2] - this->x[1]; // (xc-xb)
  yac = this->y[0] - this->y[2]; // (ya-yc)
  xxc = x - this->x[2];          // (x -xc)
  yyc = y - this->y[2];          // (y -yc)
  den = ybc*xac + xcb*yac;
  wa = ( ybc*xxc+xcb*yyc)/den;
  wb = (-yac*xxc+xac*yyc)/den;
  wc = 1.0 - wa - wb;
  if( 0.0 <= wa && wa <= 1.0 && 0.0 <= wb && wb <= 1.0 && 0.0 <= wc && wc <= 1.0 ){
    value = wa*z[0] + wb*z[1] + wc*z[2];
    return true;
  } else {
    return false;
  }
}

// Class: Cell
//===============================================================================================================
void Cell::find_bounding_box(){
  xmin = this->x[0];
  xmax = this->x[0];
  ymin = this->y[0];
  ymax = this->y[0];
  for(int i=0;i<this->Nc;i++){
    if( x[i] < xmin ){
      xmin = x[i];
    }
    if( x[i] > xmax ){
      xmax = x[i];
    }
    if( y[i] < ymin ){
      ymin = y[i];
    }
    if( y[i] > ymax ){
      ymax = y[i];
    }
  }
}

bool Cell::point_in_bounding_box(double x,double y){
  if( x < xmin || xmax < x || y < ymin || ymax < y ){
    return false;
  } else {
    return true;
  }
}

bool Cell::point_in_cell(double x,double y){
  if( this->point_in_bounding_box(x,y) ){
    // use ray casting
    int intersections = 0;
    for(int i=0;i<this->Nc-1;i++){
      if( this->sides_are_intersecting(x,y,xmax+xmax/2.0,y,x[i],y[i],x[i+1],y[i+1]) ){
	intersections++;
      }
    }
    if( this->sides_are_intersecting(x,y,xmax+xmax/2.0,y,x[Nc-1],y[Nc-1],x[0],y[0]) ){
      intersections++;
    }
    if( (intersections % 2) == 0 ){
      return false;
    } else {
      return true;
    }
  } else {
    return false;
  }
}

bool Cell::sides_are_intersecting(double v1x1, double v1y1, double v1x2, double v1y2,double v2x1, double v2y1, double v2x2, double v2y2){
  double d1, d2;
  double a1, a2, b1, b2, c1, c2;
  // Convert vector 1 to a line (line 1) of infinite length. We want the line in linear equation standard form: A*x + B*y + C = 0. See: http://en.wikipedia.org/wiki/Linear_equation
  a1 = v1y2 - v1y1;
  b1 = v1x1 - v1x2;
  c1 = (v1x2 * v1y1) - (v1x1 * v1y2);  
  // Every point (x,y), that solves the equation above, is on the line, every point that does not solve it, is not.
  // The equation will have a positive result if it is on one side of the line and a negative one if is on the other side of it. We insert (x1,y1) and (x2,y2) of vector 2 into the equation above.
  d1 = (a1 * v2x1) + (b1 * v2y1) + c1;
  d2 = (a1 * v2x2) + (b1 * v2y2) + c1;  
  // If d1 and d2 both have the same sign, they are both on the same side of our line 1 and in that case no intersection is possible.
  // Careful, 0 is a special case, that's why we don't test ">=" and "<=", but "<" and ">".
  if( d1 > 0 && d2 > 0 ){
    return false;
  }
  if( d1 < 0 && d2 < 0 ){
    return false;
  }  
  // The fact that vector 2 intersected the infinite line 1 above doesn't mean it also intersects the vector 1.
  // Vector 1 is only a subset of that infinite line 1, so it may have intersected that line before the vector started or after it ended.
  // To know for sure, we have to repeat the same test the other way round. We start by calculating the infinite line 2 in linear equation standard form.
  a2 = v2y2 - v2y1;
  b2 = v2x1 - v2x2;
  c2 = (v2x2 * v2y1) - (v2x1 * v2y2);  
  // Calculate d1 and d2 again, this time using points of vector 1.
  d1 = (a2 * v1x1) + (b2 * v1y1) + c2;
  d2 = (a2 * v1x2) + (b2 * v1y2) + c2;  
  // Again, if both have the same sign (and neither one is 0), no intersection is possible.
  if( d1 > 0 && d2 > 0 ){
    return false;
  }
  if( d1 < 0 && d2 < 0 ){
    return false;
  }  
  // If we get here, only two possibilities are left. Either the two vectors intersect in exactly one point or they are collinear, which means they intersect in any number of points from zero to infinite.
  if( (a1 * b2) - (a2 * b1) == 0.0 ){
    return false;
  }  
  // If they are not collinear, they must intersect in exactly one point.
  return true;
}

double Cell::area(){
  double sum = 0.0;
  for(int i=0;i<this->Nc-1;i++){
    sum += this->x[i]*this->y[i+1] - this->y[i]*this->x[i+1];
  }
  sum += this->x[Nc-1]*this->y[0] - this->y[Nc-1]*this->x[0];
  return sum/2.0;
}

