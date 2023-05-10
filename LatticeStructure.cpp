#include <LatticeStructure.h>
#include <IoData.h>
#include <Utils.h>
#include <iostream>
#include <cmath>
using std::vector;

//-------------------------------------------------------------------------

LatticeStructure::LatticeStructure() : lattice_id(-1)
{ 
  //nothing is done here. must call function "Setup" to setup.
}

//-------------------------------------------------------------------------

LatticeStructure::~LatticeStructure()
{ }

//-------------------------------------------------------------------------

void
LatticeStructure::Setup(int lattice_id_, LatticeData &iod_lattice, int nMaterials)
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

    site_matid[site_id] = it->second->materialid;
    if(site_matid[site_id]<0 || site_matid[site_id]>=nMaterials) {
      print_error("*** Error: Detected non-existent material id (%d) in Lattice[%d]->Site[%d].\n",
                  site_matid[site_id], lattice_id, site_id);
      exit_mpi();
    }
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
      labc[p] -= ijk[p];
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
      labc0[p] = labc[p] - ijk[p];
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

