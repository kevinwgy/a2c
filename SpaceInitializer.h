#ifndef _SPACE_INITIALIZER_H_
#define _SPACE_INITIALIZER_H_

#include<IoData.h>
#include<LatticeStructure.h>
#include<MaterialOperator.h>
#include<utility> //std::pair
#include<map>
#include<mpi.h>
class LatticeVariables;
namespace GeoTools {
  class DistanceFromPointToPlane;
  class DistanceFromPointToCylinderCone;
  class DistanceFromPointToCylinderSphere;
  class DistanceFromPointToSphere;
  class DistanceFromPointToParallelepiped;
  class DistanceFromPointToSpheroid;
};
class UserDefinedGeometry;

/****************************************************************************
 * class SpaceInitializer is responsible for setting up the computational domain
 * and initializing the lattice variables.
 ****************************************************************************/

class SpaceInitializer {

  IoData& iod;
  MPI_Comm& comm;

  std::vector<LatticeStructure>& lats;
  std::vector<MaterialOperator>& mato;

public:

  SpaceInitializer(MPI_Comm& comm_, IoData& iod_, std::vector<MaterialOperator> &mato_,
                   std::vector<LatticeStructure> &lats_);

  ~SpaceInitializer();

  void SetupLatticeVariables(std::vector<LatticeVariables> &LVS);


private:

  //! Initializes lattice variables (LV)
  void CreateOneLattice(LatticeStructure &lat, LatticeData &iod_lat, LatticeVariables &LV);

  //! Apply user-specified regional IC
  void ApplyRegionalICOnOneLattice(RegionalIcData& ic, LatticeVariables &LV);

  //! Find a bounding box aligned with lattice axes.
  void FindLatticeDomainBoundingBox(LatticeStructure &lat, LatticeData &iod_lat,
                                    UserDefinedGeometry *geom_specifier, Vec3D &lmin, Vec3D &lmax);

  //! Check if sites in a lattice site (li,lj,lk) are inside domain. If yes, add to LV.
  int CheckAndAddSitesInCell(int li, int lj, int lk, LatticeStructure &lat, LatticeData &iod_lat,
                             std::vector<std::pair<int,int> > &order,
                             std::map<int,GeoTools::DistanceFromPointToPlane*> &plane_cal,
                             std::map<int,GeoTools::DistanceFromPointToCylinderCone*> &cylindercone_cal,
                             std::map<int,GeoTools::DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                             std::map<int,GeoTools::DistanceFromPointToSphere*> &sphere_cal,
                             std::map<int,GeoTools::DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                             std::map<int,GeoTools::DistanceFromPointToSpheroid*> spheroid_cal,
                             UserDefinedGeometry *geom_specifier,
                             LatticeVariables &LV);

  //! Find the order of set operations for user-specified geometries (spheres, cylinders, etc.)
  //! Also, setup (signed) distance calculators
  //! This templated function works for IoDataGeo = LatticeData & RegionalIcData
  template<class IoDataGeo>
  void SortUserSpecifiedGeometries(IoDataGeo &iod_geo,
                                   std::vector<std::pair<int,int> > &order,
                                   std::map<int,GeoTools::DistanceFromPointToPlane*> &plane_cal,
                                   std::map<int,GeoTools::DistanceFromPointToCylinderCone*> &cylindercone_cal,
                                   std::map<int,GeoTools::DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                                   std::map<int,GeoTools::DistanceFromPointToSphere*> &sphere_cal,
                                   std::map<int,GeoTools::DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                                   std::map<int,GeoTools::DistanceFromPointToSpheroid*> spheroid_cal);
};


#endif
