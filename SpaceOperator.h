#ifndef _SPACE_OPERATOR_H_
#define _SPACE_OPERATOR_H_

#include<IoData.h>
#include<LatticeInfo.h>
#include<MaterialOperator.h>
#include<mpi.h>
class LatticeVariables;

/****************************************************************************
 * class SpaceOperator is responsible for calculations in the spatial domain
 ***************************************************************************/

class SpaceOperator {

  IoData& iod;
  MPI_Comm& comm;

  std::vector<MaterialOperator>& mato;
  std::vector<LatticeInfo>& linfo;


public:

  SpaceOperator(MPI_Comm& comm_, IoData& iod_, std::vector<MaterialOperator> &mato_,
                std::vector<LatticeInfo> &linfo_);
  ~SpaceOperator();

  void SetupLatticeVariables(std::vector<LatticeVariables> &LVS);

private:

  void CreateSpecimen(std::vector<LatticeVariables> &LVS); //!< initializes "q"


};

#endif
