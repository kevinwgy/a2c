#include <LatticeStructure.h>
#include <MaterialOperator.h>
#include <Utils.h>
#include <iostream>
#include <cmath>
#include <set>
using std::vector;

//-------------------------------------------------------------------------

LatticeStructure::LatticeStructure() : lattice_id(-1), dmin(0.0)
{ 
  //nothing is done here. must call function "Setup" to setup.
}

//-------------------------------------------------------------------------

LatticeStructure::~LatticeStructure()
{ }

//-------------------------------------------------------------------------

void
LatticeStructure::Setup(int lattice_id_, LatticeData &iod_lattice, vector<MaterialOperator> &mato)
{
  lattice_id = lattice_id_;

  abc[0][0] = iod_lattice.ax;
  abc[0][1] = iod_lattice.ay;
  abc[0][2] = iod_lattice.az;

  abc[1][0] = iod_lattice.bx;
  abc[1][1] = iod_lattice.by;
  abc[1][2] = iod_lattice.bz;

  abc[2][0] = iod_lattice.cx;
  abc[2][1] = iod_lattice.cy;
  abc[2][2] = iod_lattice.cz;

  o[0] = iod_lattice.ox;
  o[1] = iod_lattice.oy;
  o[2] = iod_lattice.oz;

  for(int p=0; p<3; p++) {
    abc0[p] = abc[p].norm();
    assert(abc0[p]>0.0);
    dirabc[p] = abc[p]/abc0[p];
  }

  alpha_beta_gamma[0] = acos(dirabc[1]*dirabc[2]);
  alpha_beta_gamma[1] = acos(dirabc[0]*dirabc[2]);
  alpha_beta_gamma[2] = acos(dirabc[0]*dirabc[1]);

  vol = fabs(abc[0]*(abc[1]^abc[2]));
  
  dmin = iod_lattice.dmin;
  if(dmin<=0.0) //user did not specify dmin. Set it to the min lattice spacing constant 
    dmin = std::min(abc0[0], std::min(abc0[1], abc0[2]));

  a_periodic = (iod_lattice.a_periodic == LatticeData::TRUE);
  b_periodic = (iod_lattice.b_periodic == LatticeData::TRUE);
  c_periodic = (iod_lattice.c_periodic == LatticeData::TRUE);

  int nSites = iod_lattice.siteMap.dataMap.size();
  if(nSites==0) {
    print_error("*** Error: No lattice sites specified for Lattice[%d].\n", lattice_id);
    exit_mpi();
  }

  site_coords.resize(nSites);
  site_matid.resize(nSites);
  site_species_id.resize(nSites);
  site_species_diff.resize(nSites);
  site_species_x0.resize(nSites);
  site_species_xmin.resize(nSites);
  
  nSpecies_max = 0;

  std::set<int> site_tracker;
  for(auto it = iod_lattice.siteMap.dataMap.begin(); it != iod_lattice.siteMap.dataMap.end(); it++) {
    int site_id = it->first;
    if(site_id<0 || site_id>=nSites) {
      print_error("*** Error: Detected error in the specification of sites in Lattice[%d]. (id=%d)\n",
                  lattice_id, site_id);
      exit_mpi();
    }
    site_tracker.insert(site_id);

    site_coords[site_id][0] = it->second->la;
    site_coords[site_id][1] = it->second->lb;
    site_coords[site_id][2] = it->second->lc;
    for(int i=0; i<3; i++) {
      if(site_coords[site_id][i]<0.0 || site_coords[site_id][i]>=1.0) {
        print_error("*** Error: Site coords must be in [0, 1). Detected %e.\n", site_coords[site_id][i]);
        exit_mpi();
      }
    }

    int matid = it->second->materialid;
    if(matid<0 || matid>=(int)mato.size()) {
      print_error("*** Error: Detected non-existent material id (%d) in Lattice[%d]->Site[%d].\n",
                  site_matid[site_id], lattice_id, site_id);
      exit_mpi();
    }
    site_matid[site_id] = matid;

    // Initialize all the species of the material
    int nMatSpecs = mato[matid].GetNumberOfSpecies();
    if(nMatSpecs>nSpecies_max)
      nSpecies_max = nMatSpecs;

    for(int i=0; i<nMatSpecs; i++)
      site_species_id[site_id].push_back(mato[matid].LocalToGlobalSpeciesID(i));
    site_species_diff[site_id].resize(nMatSpecs, 0); //by default, not diffuisive
    site_species_x0[site_id].resize(nMatSpecs, mato[matid].GetMinMolarFraction());
    site_species_xmin[site_id].resize(nMatSpecs, mato[matid].GetMinMolarFraction());

    std::set<int> checker1, checker2; //for error detection only
    for(auto&& sp : it->second->speciesMap.dataMap) {
      int global_spid = sp.second->speciesID;
      int local_spid  = mato[matid].GlobalToLocalSpeciesID(global_spid);
      if(local_spid<0) {//couldn't find thie species
        print_error("*** Error: Lattice[%d]->Site[%d] has species %d, which is not in its material (%d).\n",
                    lattice_id, site_id, global_spid, matid);
        exit_mpi();
      }
      site_species_diff[site_id][local_spid] = sp.second->diffusive==LocalSpeciesData::TRUE ?
                                               1 : 0;
      if(sp.second->molar_fraction>=0.0) //user specified this
        site_species_x0[site_id][local_spid] = sp.second->molar_fraction;
      if(sp.second->min_molar_fraction>=0.0) //user specified this
        site_species_xmin[site_id][local_spid] = sp.second->min_molar_fraction;

      checker1.insert(global_spid);
      checker2.insert(sp.first);
    }
    if(checker1.size() != it->second->speciesMap.dataMap.size()) {
      print_error("*** Error: Lattice[%d]->Site[%d] has duplicate (global) species.\n",
                  lattice_id, site_id);
      exit_mpi();
    }
    if(checker2.size() != it->second->speciesMap.dataMap.size()) {
      print_error("*** Error: Lattice[%d]->Site[%d] has duplicate LocalSpecies indices.\n",
                  lattice_id, site_id);
      exit_mpi();
    }

  }
  if((int)site_tracker.size() != nSites) {
    print_error("*** Error: Detected duplicate Site IDs in Lattice[%d].\n", lattice_id);
    exit_mpi();
  }

}

//-------------------------------------------------------------------------

Vec3D
LatticeStructure::GetLocalCoordsFromXYZ(Vec3D &xyz, Int3 *ijk)
{
  Vec3D labc = GetLABC(xyz);
  if(ijk) {
    for(int p=0; p<3; p++) {
      ijk[p] = floor(labc[p]);
      labc[p] -= (*ijk)[p];
    }
  } else {
    for(int p=0; p<3; p++)
      labc[p] -= floor(labc[p]);
  }
  return labc;
}

//-------------------------------------------------------------------------

Vec3D
LatticeStructure::GetLocalCoordsFromLABC(Vec3D &labc, Int3 *ijk)
{
  Vec3D labc0;
  if(ijk) {
    for(int p=0; p<3; p++) {
      ijk[p] = floor(labc[p]);
      labc0[p] = labc[p] - (*ijk)[p];
    }
  }
  else { 
    for(int p=0; p<3; p++)
      labc0[p] = labc[p] - floor(labc[p]);
  }   
  return labc0;
}

//-------------------------------------------------------------------------




//-------------------------------------------------------------------------

