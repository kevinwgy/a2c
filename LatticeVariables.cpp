#include<LatticeVariables.h>
#include<Utils.h>
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
  std::vector<bool> stay(size, true);
  assert(size==(int)q.size());
  for(int i=0; i<size; i++) {
    if(tree.findCandidatesWithin(q[i], &candidate, 1, dmin)>0) {
      stay[i] = false;
      nErased++;
    }
  }
 
  if(nErased==0) //nothing to erase
    return nErased;

  size = size - nErased;

  // reset the vectors
  std::vector<int>                  siteid_2     = siteid;      siteid.resize(size);
  std::vector<int>                  matid_2      = matid;       matid.resize(size);
  std::vector<Vec3D>                l_2          = l;           l.resize(size);
  std::vector<Vec3D>                q_2          = q;           q.resize(size);
  std::vector<Vec3D>                q0_2         = q0;          q0.resize(size);
  //std::vector<Vec3D>                p_2          = p;           p.resize(size);
  std::vector<double>               sigma_2      = sigma;       sigma.resize(size);
  std::vector<int>                  tag_2        = tag;         tag.resize(size);
  std::vector<std::vector<int> >    diffusive_2  = diffusive;   diffusive.resize(size);
  std::vector<std::vector<int> >    species_id_2 = species_id;  species_id.resize(size);
  std::vector<std::vector<double> > x_2          = x;           x.resize(size);
  std::vector<std::vector<double> > gamma_2      = gamma;       gamma.resize(size);
  
  int counter = 0;
  for(int i=0; i<(int)stay.size(); i++) {
    if(stay[i]) {
      siteid[counter]     = siteid_2[i];
      matid[counter]      = matid_2[i];
      l[counter]          = l_2[i];
      q[counter]          = q_2[i];
      q0[counter]         = q0_2[i];
//      p[counter]          = p_2[i];
      sigma[counter]      = sigma_2[i];
      tag[counter]        = tag_2[i];
      diffusive[counter]  = diffusive_2[i];
      species_id[counter] = species_id_2[i];
      x[counter]          = x_2[i];
      gamma[counter]      = gamma_2[i];

      counter++;
    }
  }

  // sanity checks
  assert(counter==size);
  found = false;
  for(auto&& xi : x) {
    if((int)xi.size() == nSpecies_max) {
      found = true;
      break;
    } 
  }
  if(!found) {
    print_error("*** Error: After removing overlapping sites on Lattice[%d], nSpecies_max has reduced.\n",
                lattice_id);
    exit_mpi();
  } 

  FindBoundingBoxAndSize(); //update lmin, lmax

  return nErased;

}


//---------------------------------------------------------------------

