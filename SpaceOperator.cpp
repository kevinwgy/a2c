#include<SpaceOperator.h>
#include<SpaceInitializer.h>
#include<LatticeVariables.h>
#include<Utils.h>
#include<cmath> //isfinite, floor
using std::vector;

//---------------------------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm& comm_, IoData& iod_, vector<MaterialOperator> &mato_,
                             vector<LatticeStructure> &lats_) 
             : iod(iod_), comm(comm_), mato(mato_), lats(lats_)
{

}

//---------------------------------------------------------------------

SpaceOperator::~SpaceOperator()
{ }

//---------------------------------------------------------------------

void 
SpaceOperator::SetupLatticeVariables(vector<LatticeVariables> &LVS)
{
  SpaceInitializer initializer(comm, iod, mato, lats); //hire a contractor to do the work :P
  initializer.SetupLatticeVariables(LVS);
} 

//---------------------------------------------------------------------


