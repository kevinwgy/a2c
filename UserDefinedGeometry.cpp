#include "UserDefinedGeometry.h"
#include<cstdio>

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ A2C.
// If multiple UserDefinedGeometry need to be defined, they should
// have different names (e.g., MyGeometry1, MyGeometry2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedGeometry.h -c UserDefinedGeometry.cpp; g++ -shared UserDefinedGeometry.o -o UserDefinedGeometry.so; rm UserDefinedGeometry.o
//------------------------------------------------------------

class MyGeometry : public UserDefinedGeometry{

public:

  //! Returns a bounding box of this set (tighter is better) in a specifieed coord. system
  //! The basis vectors may not be orthogonal to each, and they are generally NOT UNIT VECTORS.
  void GetBoundingBox(double* O, double* dir0, double* dir1, double* dir2, //inputs
                      double* lmin, double *lmax); //outputs

  //! Checks whether a point is inside the set
  //! x: current position of the point (in x,y,z coords)
  //! x0: original position of the point (in x,y,z coords)
  //! O & dirs: lattice origin and vectors.
  bool IsPointInside(double *x, double *x0, double* O, double* dir0, double* dir1, double* dir2);

  //! Calculates signed distance to a point (inside: <0, outside: >0)
  double GetDistanceToPoint(double *x, double *x0, double* O, double* dir0, double* dir1, double* dir2);


};

//------------------------------------------------------------

void
MyGeometry::GetBoundingBox(double* O, double* dir0, double* dir1, double* dir2, //inputs
                           double* lmin, double *lmax); //outputs
{
  //TODO: The user should complete this function
  fprintf(stderr,"*** Error: MyGeometry::GetBoundingBox has not been implemented.\n");
  lmin[0] = lmin[1] = lmin[2] = -1.0;
  lmax[0] = lamx[1] = lamx[2] = 1.0;
    
}

//------------------------------------------------------------

bool
MyGeometry::IsPointInside(double *x, double *x0, double* O, double* dir0, double* dir1, double* dir2)
{
  //TODO: The user should complete this function
  fprintf(stderr,"*** Error: MyGeometry::IsPointInside has not been implemented.\n");
  return false;
}

//------------------------------------------------------------

double
MyGeometry::GetDistanceToPoint(double *x, double *x0, double* O, double* dir0, double* dir1, double* dir2)
{
  //TODO: The user should complete this function, if it is used in the solver. Otherwise, leave it as is.
  fprintf(stderr,"*** Error: MyGeometry::GetDistanceToPoint has not been implemented.\n");
  return 
}


//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "MyGeometry".)
extern "C" UserDefinedGeometry* Create() {
  return new MyGeometry; //TODO: If you've changed the name of the derived class, this needs to be updated.
}

extern "C" void Destroy(UserDefinedGeometry* udg) {
  delete udg;
}

//------------------------------------------------------------

