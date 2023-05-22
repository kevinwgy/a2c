#include<LatticeVariables.h>
#include<KDTree.h>
#include<cmath> //isfinite
#include<cfloat>
#include<cassert>
using std::isfinite;

//---------------------------------------------------------------------

LatticeVariables::LatticeVariables() : lmin(0.0), lmax(0.0), size(0), lattice_id(-1)
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

  lattice_id = -1;

  siteid.resize(0);
  matid.resize(0); 

  l.resize(0);
  q.resize(0);
  q0.resize(0);
//  p.resize(0);
  sigma.resize(0);
  tag.resize(0);

  diffusive.resize(0);
  species_id.resize(0);
  x.resize(0);
  gamma.resize(0);
}

//---------------------------------------------------------------------

void
LatticeVariables::FindBoundingBoxAndSize()
{
  lmin = DBL_MAX;
  lmax = -DBL_MAX;
  for(auto&& l0 : l) {
    for(int i=0; i<3; i++) {
      if(l0[i]<lmin[i])
        lmin[i] = l0[i];
      if(l0[i]>lmax[i])
        lmax[i] = l0[i];
    }
  }

  size = l.size();
  assert((int)siteid.size() == size);
  assert((int)matid.size() == size);
  assert((int)q.size() == size);
  assert((int)q0.size() == size);
  assert((int)sigma.size() == size);
  assert((int)tag.size() == size);
  assert((int)diffusive.size() == size);
  assert((int)species_id.size() == size);
  assert((int)x.size() == size);
  assert((int)gamma.size() == size);
}

//---------------------------------------------------------------------

int
LatticeVariables::RemoveConflictingSites(std::vector<Vec3D> &q2, double dmin)
{

  int nErased(0);

  //------------------------------------
  // Step 1: Storing q2 in a KDTree
  //------------------------------------
  std::vector<PointIn3D> q2_copy;
  q2_copy.reserve(q2.size());
  for(int i=0; i<(int)q2.size(); i++)
    q2_copy.push_back(PointIn3D(i, q2[i]));

  KDTree<PointIn3D,3> tree(q2_copy.size(), q2_copy.data());


  //------------------------------------
  // Step 2: Nearest-neighbor search
  //------------------------------------
  PointIn3D candidate;
  int found;

  auto it_siteid = siteid.begin();
  auto it_matid  = matid.begin();
  auto it_l      = l.begin();
  auto it_q      = q.begin();
  auto it_q0     = q0.begin();
  auto it_sigma  = sigma.begin();
  auto it_tag    = tag.begin();
  auto it_diff   = diffusive.begin();
  auto it_spid   = species_id.begin();
  auto it_x      = x.begin();
  auto it_gamma  = gamma.begin();

  while(it_q != q.end()) {
    found = tree.findCandidatesWithin(*it_q, &candidate, 1, dmin);
    if(found) { //erase it
      it_siteid = siteid.erase(it_siteid); //iterator moves to the next element
      it_matid  = matid.erase(it_matid);
      it_l      = l.erase(it_l);
      it_q      = q.erase(it_q);
      it_q0     = q0.erase(it_q0);
      it_sigma  = sigma.erase(it_sigma);
      it_tag    = tag.erase(it_tag);
      it_diff   = diffusive.erase(it_diff);
      it_spid   = species_id.erase(it_spid);
      it_x      = x.erase(it_x);
      it_gamma  = gamma.erase(it_gamma);

      nErased++;

    } else {
      it_siteid ++;
      it_matid  ++;
      it_l      ++;
      it_q      ++;
      it_q0     ++;
      it_sigma  ++;
      it_tag    ++;
      it_diff   ++;
      it_spid   ++;
      it_x      ++;
      it_gamma  ++; 
    }
  }

  return nErased;

}


//---------------------------------------------------------------------

