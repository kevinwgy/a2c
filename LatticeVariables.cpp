#include<LatticeVariables.h>
#include<cmath> //isfinite
using std::isfinite;

//---------------------------------------------------------------------

LatticeVariables::LatticeVariables() : lmin(0.0), lmax(0.0), size(0)
{ }

//---------------------------------------------------------------------

LatticeVariables::~LatticeVariables()
{ }

//---------------------------------------------------------------------

void
LatticeVariables::Clear()
{
  lmin = lmax = 0.0;

  size = 0;

  siteid.resize(0);
  matid.resize(0); 

  l.resize(0);
  q.resize(0);
  q0.resize(0);
//  p.resize(0);
  sigma.resize(0);
  subsurf.resize(0);

  diffusive.resize(0);
  species_id.resize(0);
  x.resize(0);
  gamma.resize(0);

  xmin.resize(0);
}

//---------------------------------------------------------------------


//---------------------------------------------------------------------


//---------------------------------------------------------------------

