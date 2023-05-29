#ifndef _LATTICE_INFO_H_
#define _LATTICE_INFO_H_
#include<Vector3D.h>
#include<vector>
#include<cassert>
#include<cmath>
struct LatticeData;
struct MaterialOperator;

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

  int nSpecies_max; //!< max number of species among the materials on this lattice

public:
  // these vectors are "public"
  std::vector<std::vector<int> > site_species_id; //!< the (global) ID of each species of each site
  std::vector<std::vector<int> > site_species_diff; //!< (0 or 1) whether each species is diffusive
  std::vector<std::vector<double> > site_species_x0; //!< initial molar fraction
  std::vector<std::vector<double> > site_species_xmin; //!< min molar fraction
  std::vector<std::vector<double> > site_species_gamma; //!< chemical potential (default is 0)


public:

  LatticeStructure();
  ~LatticeStructure();

  void Setup(int lattice_id_, LatticeData &iod_lattice, std::vector<MaterialOperator> &mato);

  inline int GetLatticeID() {return lattice_id;}
  inline Vec3D GetLatticeOrigin() {return o;}
  inline Vec3D GetLatticeVector(int p) {assert(p>=0 && p<=2); return abc[p];}
  inline Vec3D GetLatticeDirection(int p) {assert(p>=0 && p<=2); return dirabc[p];}
  inline double GetLatticeSpacing(int p) {assert(p>=0 && p<=2); return abc0[p];}
  inline double GetLatticeAngle(int p) {assert(p>=0 && p<=2); return alpha_beta_gamma[p];}
  inline double GetUnitCellVolume() {return vol;}

  inline double GetDmin() {return dmin;}

  inline int GetNumberOfSites() {return site_coords.size();}
  inline Vec3D GetCoordsOfSite(int p) {assert(p>=0 && p<(int)site_coords.size()); return site_coords[p];}
  inline int GetMatIDOfSite(int p) {assert(p>=0 && p<(int)site_matid.size()); return site_matid[p];}

  inline int GetMaxNumberOfSpeciesPerSite() {return nSpecies_max;}


  //! functions that switch from (x,y,z) and (la,lb,lc) coords
  inline Vec3D GetXYZ(double la, double lb, double lc) {return o+la*abc[0]+lb*abc[1]+lc*abc[2];}
  inline Vec3D GetXYZ(Vec3D& labc) {return GetXYZ(labc[0], labc[1], labc[2]);}
  inline Vec3D GetLABC(Vec3D& xyz) {
    Vec3D labc; for(int p=0; p<3; p++) labc[p] = (xyz-o)*dirabc[p]/abc0[p]; return labc;}
  inline Vec3D GetLABC(double x, double y, double z) {Vec3D xyz(x,y,z); return GetLABC(xyz);}

  //! Get cell id & local lattice coords
  Vec3D GetLocalCoordsFromXYZ(Vec3D &xyz, Int3 *ijk = NULL); 
  Vec3D GetLocalCoordsFromLABC(Vec3D &labc, Int3 *ijk = NULL);

};

#endif
