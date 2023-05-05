#ifndef _LATTICE_INFO_H_
#define _LATTICE_INFO_H_
#include<Vector3D.h>
#include<vector>
#include<cassert>
struct LatticeData;

/****************************************************************
 * Class LatticeInfo stores the information (not variables) of
 * a lattice. Ref: Chap 3 of Tadmor and Miller, ``Modeling Materials:
 * Continuum, Atomistic and Multiscale Techniques''
 ***************************************************************/

class LatticeInfo
{

  int lattice_id;

  Vec3D o; //!< lattice origin
  Vec3D abc[3]; //!< lattice vectors (not just directions)
  Vec3D dirabc[3]; //!< normalized version of a, b, c;
  double a0, b0, c0, alpha, beta, gamma; //!< lattice constants
  double vol; //!< volume of unit cell (|a.(bxc)|)

  bool a_periodic, b_periodic, c_periodic; //!< whether the periodic b.c. holds in a, b, c dir. 

  double dmin; //!< min distance between sites on this lattice and those on other lattices

  std::vector<Vec3D> site_coords; //!< lattice coords [0, 1) of each site
  std::vector<int>   site_matid; //!< material id of each site

public:

  LatticeInfo();
  ~LatticeInfo();

  void Setup(int lattice_id_, LatticeData &iod_lattice, int nMaterials);

  Vec3D GetLatticeVector(int i) {assert(i>=0 && i<=2); return abc[i];}
  Vec3D GetLatticeDirection(int i) {assert(i>=0 && i<=2); return dirabc[i];}
  double GetLatticeSpacing(int i) {assert(i>=0 && i<=2); return i==0 ? a0 : (i==1 ? b0 : c0);}
  double GetLatticeAngle(int i) {assert(i>=0 && i<=2); return i==0 ? alpha : (i==1 ? beta : gamma);}
  double GetUnitCellVolume() {return vol;}

  Vec3D GetXYZ(double i, double j, double k) {return i*abc[0]+j*abc[1]+k*abc[2];}
  Vec3D GetXYZ(Vec3D& ijk) {return GetXYZ(ijk[0], ijk[1], ijk[2]);}



}

#endif
