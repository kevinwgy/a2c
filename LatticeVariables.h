#ifndef _LATTICE_VARIABLES_H_
#define _LATTICE_VARIABLES_H_
#include <vector>
#include <Vector3D.h>

/**********************************************
 * class Variables stores all the unknown
 * spatial variables for one lattice. It is 
 * mainly a storage class, not responsible for 
 * computations.
 *********************************************/

class LatticeVariables {

public:

  Vec3D lmin, lmax; //!< min and max lattice coords

  int size; //!< size of all the vectors below

  std::vector<int> siteid; //!< the site group id of each site
  std::vector<int> matid; //!< the material id of each 

  std::vector<Vec3D> l; //!< lattice coordinates (non-dimensional)
  std::vector<Vec3D> q; //!< mean position (NOT displacement)
  std::vector<Vec3D> q0; //!< initial mean position
//  std::vector<Vec3D> p; //!< mean momentum
  std::vector<double> sigma; //!< atomic frequency
  std::vector<int> subsurf; //!< 0 ~ not within subsurface; 1 ~ within subsurface

  // Variables defined for each species
  std::vector<std::vector<int> > diffusive; //!< (1/0) participates in mass exchange or not
  std::vector<std::vector<int> > species_id; //!< species id
  std::vector<std::vector<double> > x; //!< molar fraction
  std::vector<std::vector<double> > gamma; //!< chemical potential

  std::vector<std::vector<double> > xmin; //!< min molar fraction (from input)


public:

  LatticeVariables();
  ~LatticeVariables();

  void Clear();

};

#endif
