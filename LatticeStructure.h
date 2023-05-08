#ifndef _LATTICE_INFO_H_
#define _LATTICE_INFO_H_
#include<Vector3D.h>
#include<vector>
#include<cassert>
struct LatticeData;

/****************************************************************
 * Class LatticeStructure stores the information (not variables) of
 * a lattice structure, and performs simple calculations (e.g., change
 * of coordinate systems). Ref: Chap 3 of Tadmor and Miller, ``Modeling 
 * Materials: Continuum, Atomistic and Multiscale Techniques''
 ***************************************************************/

class LatticeStructure
{

  int lattice_id;

  Vec3D o; //!< lattice origin
  Vec3D abc[3]; //!< lattice vectors (not just directions)
  Vec3D dirabc[3]; //!< normalized version of a, b, c;
  Vec3D abc0, alpha_beta_gamma; //!< lattice constants
  double vol; //!< volume of unit cell (|a.(bxc)|)

  bool a_periodic, b_periodic, c_periodic; //!< whether the periodic b.c. holds in a, b, c dir. 

  double dmin; //!< min distance between sites on this lattice and those on other lattices

  std::vector<Vec3D> site_coords; //!< lattice coords [0, 1) of each site
  std::vector<int>   site_matid; //!< material id of each site

public:

  LatticeStructure();
  ~LatticeStructure();

  void Setup(int lattice_id_, LatticeData &iod_lattice, int nMaterials);

  inline Vec3D GetLatticeVector(int p) {assert(p>=0 && p<=2); return abc[p];}
  inline Vec3D GetLatticeDirection(int p) {assert(p>=0 && p<=2); return dirabc[p];}
  inline double GetLatticeSpacing(int p) {assert(p>=0 && p<=2); return abc0[p];}
  inline double GetLatticeAngle(int p) {assert(p>=0 && p<=2); return alpha_beta_gamma[p];}
  inline double GetUnitCellVolume() {return vol;}

  //! functions that switch from (x,y,z) and (i,j,k) coords
  inline Vec3D GetXYZ(double i, double j, double k) {return i*abc[0]+j*abc[1]+k*abc[2];}
  inline Vec3D GetXYZ(Vec3D& ijk) {return GetXYZ(ijk[0], ijk[1], ijk[2]);}
  inline Vec3D GetIJK(Vec3D& xyz) {
    Vec3D ijk; for(int p=0; p<3; p++) ijk[p] = xyz*dirabc[p]/abc0[p]; return ijk;}
  inline Vec3D GetIJK(double x, double y, double z) {Vec3D xyz(x,y,z); return GetIJK(xyz);}

  //! Get cell id & local lattice coords
  inline 

}

#endif
