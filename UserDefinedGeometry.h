#ifndef _USER_DEFINED_GEOMETRY_H_
#define _USER_DEFINED_GEOMETRY_H_

/*****************************************************************
 * class UserDefinedGeometry is an interface for user to specify
 * a bounded open set in R^3. 
 ****************************************************************/
class UserDefinedGeometry {

public:

  UserDefinedGeometry() {}
  virtual ~UserDefinedGeometry() {}

  //! Returns a bounding box of this set (tighter is better) in a specifieed coord. system
  //! The basis vectors may not be orthogonal to each, and they are generally NOT UNIT VECTORS.
  virtual void GetBoundingBox(double* O, double* dir0, double* dir1, double* dir2, //inputs
                              double* lmin, double *lmax); //outputs

  //! Checks whether a point is inside the set
  //! x: current position of the point (in x,y,z coords)
  //! x0: original position of the point (in x,y,z coords)
  //! O & dirs: lattice origin and vectors.
  virtual bool IsPointInside(double *x, double *x0, double* O, double* dir0, double* dir1, double* dir2);

  //! Calculates signed distance to a point (inside: <0, outside: >0)
  virtual double GetDistanceToPoint(double *x, double *x0, double* O, double* dir0, double* dir1, double* dir2);

};

typedef UserDefinedGeometry* CreateUDG();
typedef void                 DestroyUDG(UserDefinedGeometry*);

#endif
