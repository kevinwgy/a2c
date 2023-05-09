#ifndef _SPACE_OPERATOR_H_
#define _SPACE_OPERATOR_H_

#include<IoData.h>
#include<LatticeStructure.h>
#include<MaterialOperator.h>
#include<utility> //std::pair
#include<mpi.h>
class LatticeVariables;

/****************************************************************************
 * class SpaceOperator is responsible for calculations in the spatial domain
 ***************************************************************************/

class SpaceOperator {

  IoData& iod;
  MPI_Comm& comm;

  std::vector<MaterialOperator>& mato;
  std::vector<LatticeStructure>& lats;


public:

  SpaceOperator(MPI_Comm& comm_, IoData& iod_, std::vector<MaterialOperator> &mato_,
                std::vector<LatticeStructure> &lats_);
  ~SpaceOperator();

  void SetupLatticeVariables(std::vector<LatticeVariables> &LVS);

private:

  //! Initializes lattice variables (LV)
  void CreateOneLattice(LatticeStructure &lat, LatticeData &iod_lat, LatticeVariables &LV); 

  //! Check if sites in a lattice site (li,lj,lk) are inside domain. If yes, add to LV.
  int CheckAndAddSitesInCell(int li, int lj, int lk, LatticeStructure &lat, LatticeData &iod_lat, 
                             std::vector<std::pair<int,int> > &order, LatticeVariables &LV);

  //! Finds the order of set operations for user-specified geometries (spheres, cylinders, etc.)
  void FindUserSpecifiedGeometryOrder(int lattice_id, LatticeData &iod_lat, 
                                      std::vector<std::pair<int,int> > &order);
  
};

#endif
