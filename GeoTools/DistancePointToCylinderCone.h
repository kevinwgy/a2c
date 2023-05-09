/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/
#pragma once

#include<Vec3D>

namespace GeoTools {

/************************************************************
 * Calculate distance from an arbitrary point in 3D to a
 * cylinder-cone..
 ***********************************************************/

class DistanceFromPointToCylinderCone {

  Vec3D x0, dir; //!< a point on the plane, and the (normalized) normal direction

I AM HERE!
public:

  //! Constructor
  DistanceFromPointToPlane(double *x0_, double *dir_) {
    for(int i=0; i<3; i++) {
      x0[i]  = x0_[i];
      dir[i] = dir_[i];
    }
    double norm = dir.norm();
    assert(norm!=0.0);
    dir /= norm;
  }

  ~Distance() {}

  //! Calculates the signed distance, including the projection point
  double Calculate(double *Q_, double *P = NULL) {
    double dist = ProjectPointToPlane(*(Vec3D *)Q_, x0, dir, true);
    if(P)
      *(Vec3D *)P = *(Vec3D *)Q_ - dist*dir;
    return dist;
  }

};





}; //end of namespace
