#include <LatticeInfo.h>
#include <IoData.h>
#include <Utils.h>
#include <iostream>
#include <cmath>
using std::vector;

//-------------------------------------------------------------------------

LatticeInfo::LatticeInfo() : lattice_id(-1)
{ 
  //nothing is done here. must call function "Setup" to setup.
}

//-------------------------------------------------------------------------

LatticeInfo::~LatticeInfo()
{ }

//-------------------------------------------------------------------------

void
LatticeInfo::Setup(int lattice_id_, LatticeData &iod_lattice, int nMaterials)
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

  a0 = abc[0].norm();
  assert(a0>0.0);
  b0 = abc[1].norm();
  assert(b0>0.0);
  c0 = abc[2].norm();
  assert(c0>0.0);

  dirabc[0] = abc[0]/a0;
  dirabc[1] = abc[1]/b0;
  dirabc[2] = abc[2]/c0;

  alpha = acos(dirabc[1]*dirabc[2]);
  beta  = acos(dirabc[0]*dirabc[2]);
  gamma = acos(dirabc[0]*dirabc[1]);

  vol = fabs(abc[0]*(abc[1]^abc[2]));
  
  dmin = iod_lattice.dmin;

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

    site_coords[site_id][0] = it->second.la;
    site_coords[site_id][1] = it->second.lb;
    site_coords[site_id][2] = it->second.lc;
    for(int i=0; i<3; i++) {
      if(site_coords[site_id][i]<0.0 || site_coords[site_id][i]>=1.0) {
        print_error("*** Error: Site coords must be in [0, 1). Detected %e.\n", site_coords[site_id][i]);
        exit_mpi();
      }
    }

    site_matid[site_id] = it->second.materialid;
    if(site_matid[site_id]<0 || site_matid[site_id]>=nMaterials) {
      print_error("*** Error: Detected non-existent material id (%d) in Lattice[%d]->Site[%d].\n",
                  site_matid[site_id], lattice_id, site_id);
      exit_mpi();
    }
  }

}

//-------------------------------------------------------------------------




//-------------------------------------------------------------------------

