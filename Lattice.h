#ifndef _SHELLS_H_
#define _SHELLS_H_

#include<Vector3D.h>
#include<IoData.h>
#include<vector>
#include<cassert>

/****************************************************************
 * Class Lattice is responsible for generating lattice sites and 
 * interstitial sites, and storing lattice-related information. 
 * Ref. Chap 3 of Tadmor and Miller, ``Modeling Materials: Continuum,
 * Atomistic and Multiscale Techniques''
 ***************************************************************/

class Lattice 
{

  IoData &iod;

  SpaceData::LatticeType lattice_type;
  SpaceData::IntersticeType inter1_type, inter2_type;
  
  Vec3D abc[3]; //!< lattice vectors (not just directions)
  Vec3D dirabc[3]; //!< normalized version of a, b, c;
  Vec3D a0_b0_c0, alpha_beta_gamma; //!< lattice constants
  double vol; //!< volume of unit cell (a.(bxc))
  
  std::vector<std::vector<Int3> > shells; //!< lattice sites in shells
  std::vector<double> radius; //!< radius (dim: length) of each shell of lattice sites
  std::vector<std::vector<Int3> > shells_inter1; //!< interstices (group 1) in shells
  std::vector<double> radius_inter1;
  std::vector<std::vector<Int3> > shells_inter2; //!< interstices (group 1) in shells
  std::vector<double> radius_inter2;

public:

  Shells(IoData &iod_);
  ~Shells();

  int GetNumberOfLatticeShells() {return shells.size();}
  int GetNumberOfInterstitialShells1() {return shells_inter1.size();}
  int GetNumberOfInterstitialShells2() {return shells_inter2.size();}

  std::vector<std::vector<Int3> >& GetLatticeShells() {return shells;}
  
  std::vector<std::vector<Int3> >& GetInterstitialShells1() {return shells_inter1;}
  std::vector<std::vector<Int3> >& GetInterstitialShells2() {return shells_inter2;}

  Vec3D GetLatticeVectorA() {return abc[0];}
  Vec3D GetLatticeVectorB() {return abc[1];}
  Vec3D GetLatticeVectorC() {return abc[2];}
  Vec3D GetLatticeSpacing() {return a0_b0_c0;}
  Vec3D GetLatticeAngles() {return alpha_beta_gamma;}
  double GetUnitCellVolume() {return vol;}

  Vec3D GetXYZ(double i, double j, double k) {return i*abc[0]+j*abc[1]+k*abc[2];}
  Vec3D GetXYZ(Vec3D& ijk) {return GetXYZ(ijk[0], ijk[1], ijk[2]);}


private:

  //! Generate shells (called by constructor)
  void GenerateShells();

  //! HCP
  void GenerateLatticeSitesHCP(int NC, std::vector<std::vector<Int3> > *shells);

  void GenerateOctahedralIntersticesHCP(int NC, std::vector<std::vector<Int3> > *shells);

  void GenerateTetrahedralIntersticesHCP(int NC, std::vector<std::vector<Int3> > *shells);

  //! FCC
  void GenerateLatticeSitesFCC(int NC, std::vector<std::vector<Int3> > *shells);

  void GenerateOctahedralIntersticesFCC(int NC, std::vector<std::vector<Int3> > *shells);

  void GenerateTetrahedralIntersticesFCC(int NC, std::vector<std::vector<Int3> > *shells);

}

#endif
