#include<SpaceOperator.h>
#include<LatticeVariables.h>
#include<Utils.h>
#include<GeoTools/DistancePointToParallelepiped.h>
#include<GeoTools/DistancePointToSpheroid.h>
#include<GeoTools/BoundingBoxes.h>
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

  int nLattices = iod.latticeMap.dataMap.size();
  if(LVS.size() != nLattices)
    LVS.resize(nLattices);

  for(auto&& l : iod.latticeMap.dataMap) {
    int lid = l.first; //lattice id
    assert(lid>=0 && lid<nLattices);
    CreateOneLattice(LVS[lid], lats[lid], l.second); 
  }

}

//---------------------------------------------------------------------

void
SpaceOperator::CreateOneLattice(vector<LatticeVariables> &LV, LatticeStructure &lat, LatticeData &iod_lat)
{

  //----------------------------------------------------------------
  // Step 0: Get lattice structure info
  //----------------------------------------------------------------
  int lattice_id = lat.GetLatticeID();
  Vec3D O = lat.GetLatticeOrigin();
  Vec3D dirabc[3];
  for(int i=0; i<3; i++)
    dirabc[i] = lat.GetLatticeDirection(i);


  //----------------------------------------------------------------
  // Step 1: Find a bounding box in terms of lattice coordinates
  //----------------------------------------------------------------
  Vec3D lmin(DBL_MAX), lmax(-DBL_MAX), lmin_local, lmax_local;

  // Go over cylinder-cones
  for(auto it = iod_lat.cylinderconeMap.dataMap.begin(); it != iod_lat.cylinderconeMap.dataMap.end(); it++) {
    if(it->second.inclusion == CylinderConeData::EXTERIOR)
      continue;
    
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();
    
    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = std::min(it->second->cone_height, Hmax); //cone's height
    
    GeoTools::GetBoundingBoxOfCylinderCone(x0, dir, R, L, H, O, dirabc[0], dirabc[1], dirabc[2],
                                           lmin_local, lmax_local, true); //dir have been normalized.

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i]; 
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over cylinder-spheres (i.e. possibly with end cap(s))
  for(auto it = iod_lat.cylindersphereMap.dataMap.begin(); it != iod_lat.cylindersphereMap.dataMap.end(); it++) {
    if(it->second.inclusion == CylinderSphereData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder & sphere's radius
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);

    GeoTools::GetBoundingBoxOfCylinderSphere(x0, dir, R, L, front_cap, back_cap, O, dirabc[0], dirabc[1],
                                             dirabc[2], lmin_local, lmax_local, true); //dir normalized

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i]; 
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over spheres
  for(auto it = iod_lat.sphereMap.dataMap.begin(); it != iod_lat.sphereMap.dataMap.end(); it++) {
    if(it->second->inclusion == SphereData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    double R = it->second->radius;

    GeoTools::GetBoundingBoxOfSphere(x0, R, O, dirabc[0], dirabc[1], dirabc[2],
                                     lmin_local, lmax_local, true); //dir normalized

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i]; 
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over parallelepipeds
  for(auto it = iod_lat.parallelepipedMap.dataMap.begin(); 
           it != iod_lat.parallelepipedMap.dataMap.end(); it++) {
    if(it->second->inclusion == ParallelepipedData::EXTERIOR)
      continue;

    Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
    Vec3D oa(it->second->ax, it->second->ay, it->second->az);  oa -= x0;
    Vec3D ob(it->second->bx, it->second->by, it->second->bz);  ob -= x0;
    Vec3D oc(it->second->cx, it->second->cy, it->second->cz);  oc -= x0;

    if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
      print_error("*** Error: Detected error in a user-specified parallelepiped. "
                  "Overlapping vertices or violation of right-hand rule.\n");
      exit_mpi();
    }

    GeoTools::GetBoundingBoxOfParallelepiped(x0, oa, ob, oc, O, dirabc[0], dirabc[1], dirabc[2],
                                             lmin_local, lmax_local, true); //dir normalized

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i];
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over spheroids
  for(auto it = iod_lat.spheroidMap.dataMap.begin(); it != iod_lat.spheroidMap.dataMap.end(); it++) {
    if(it->second->inclusion == SpheroidData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
    double semi_length = it->second->semi_length;
    double r = it->second->radius;

    GeoTools::GetBoundingBoxOfSpheroid(x0, axis, it->second->semi_length, it->second->radius, O,
                                       dirabc[0], dirabc[1], dirabc[2], lmin_local, lmax_local, true);

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i];
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  if(lmin[0]>=lmax[0] || lmin[1]>=lmax[1] || lmin[2]>=lmax[2]) {
    print_error("*** Error: Unable to find a bounding box for lattice domain [%d]. Domain is unbounded.\n",
                lattice_id);
    exit_mpi();
  }


  //----------------------------------------------------------------
  // Step 2: Create lattice sites (l, q), siteid, and matid
  //----------------------------------------------------------------
  I AM HERE 




}

//---------------------------------------------------------------------


//---------------------------------------------------------------------


//---------------------------------------------------------------------

initializeStateVariables(Input &input, vector<vector<Int3> > &SS1, vector<vector<Int3> > &SS2,
                              vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, vector<double> &sigma_Pd0,
                              vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, 
                              vector<int> &PdSubsurf, vector<int> &HSubsurf, 
                              double &gammabd, vector<int> &full_H) 
{
    // initialize boundary conditions
    double &T = input.file.T0;
    const double kB = 8.6173324e-5; // (eV/K)
    //const double cutoff = 1.0e-20; // The minima of hydrogen atomic fraction
    gammabd = input.file.mubd/(kB*T);
    double thickness = input.file.t_subsurf; // (angstrom) thickness of subsurface layer

  if(!(q_Pd0.empty() && q_H0.empty() && sigma_Pd0.empty() && sigma_H0.empty() && x0.empty() && gamma.empty() && PdSubsurf.empty() && HSubsurf.empty())) {
    cerr << "WARNING: some vectors are not empty before initialization!" << endl;
    q_Pd0.clear();
    q_H0.clear();
    sigma_Pd0.clear();
    sigma_H0.clear();
    x0.clear();
    gamma.clear();
    PdSubsurf.clear();
    HSubsurf.clear();
  }

  // INITIALIZE q_Pd0
  q_Pd0.push_back(Vec3D(0.0,0.0,0.0)); //first atomic site

  int p = 1; // store index 
  // loop through shells to add sites
  for(int i=0; i<(int)SS1.size(); i++)
    for(int j=0; j<(int)SS1[i].size(); j++) {
      Int3 &site = SS1[i][j];
      int N = input.file.N;
      if(input.file.sample_shape == 1) { //cube
        if(site[0]>-N && site[0]<=N && site[1]>-N && site[1]<=N && site[2]>-N && site[2]<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_Pd0.push_back(q_site);
          if(site[0]<-(N-1-2.0*thickness/input.file.a_Pd) || site[0]>(N-2.0*thickness/input.file.a_Pd) || site[1]<-(N-1-2.0*thickness/input.file.a_Pd) || site[1]>(N-2.0*thickness/input.file.a_Pd) || site[2]<-(N-1-2.0*thickness/input.file.a_Pd) || site[2]>(N-2.0*thickness/input.file.a_Pd))
            PdSubsurf.push_back(p);
          p++;
        }
      }
      else if(input.file.sample_shape == 2) { //octahedron
        if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_Pd0.push_back(q_site);
          if(!(fabs(site[0])+fabs(site[1])+fabs(site[2])<=(N-2.0*thickness/input.file.a_Pd*1.7321)))
            PdSubsurf.push_back(p);
          p++;
        }   
      }
      else if(input.file.sample_shape == 3) { // sphere
        if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_Pd0.push_back(q_site);
          if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] >= (N-2.0*thickness/input.file.a_Pd)*(N-2.0*thickness/input.file.a_Pd))
            PdSubsurf.push_back(p);
          p++;
        }
      }
    }
  
  // INITIALIZE q_H0
  p = 0; // store index 
  for(int i=0; i<(int)SS2.size(); i++)
    for(int j=0; j<(int)SS2[i].size(); j++) {
      Int3 &site = SS2[i][j];
      int N = input.file.N;
      if(input.file.sample_shape == 1) { //cube
        if(site[0]>-N && site[0]<=N && site[1]>-N && site[1]<=N && site[2]>-N && site[2]<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_H0.push_back(q_site);
          if(site[0]<-(N-1-2.0*thickness/input.file.a_Pd) || site[0]>(N-2.0*thickness/input.file.a_Pd) || site[1]<-(N-1-2.0*thickness/input.file.a_Pd) || site[1]>(N-2.0*thickness/input.file.a_Pd) || site[2]<-(N-1-2.0*thickness/input.file.a_Pd) || site[2]>(N-2.0*thickness/input.file.a_Pd))
            HSubsurf.push_back(p);
          p++;
        }
      }
      else if(input.file.sample_shape == 2) { //octahedron
        if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_H0.push_back(q_site);
          if(!(fabs(site[0])+fabs(site[1])+fabs(site[2])<=(N-2.0*thickness/input.file.a_Pd*1.7321)))
            HSubsurf.push_back(p);
          p++;
        }
      }
      else if(input.file.sample_shape == 3) { //sphere
        if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_H0.push_back(q_site);
          if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] >= (N-2.0*thickness/input.file.a_Pd)*(N-2.0*thickness/input.file.a_Pd))
            HSubsurf.push_back(p);
          p++;
        }
      }
    }

  int nPd = q_Pd0.size();
  int nH  = q_H0.size();

  // Initialize sigma, x, gamma
  sigma_Pd0.assign(nPd, input.file.sigma_Pd);
  sigma_H0.assign(nH, input.file.sigma_H);
  x0.assign(nH, input.file.xe);
  gamma.assign(nH, 0.0);
  full_H.assign(nH, 0);

/*
  // debug
  for(int i=0; i<nPd; i++) 
    sigma_Pd0[i] = (double)i;
  for(int i=0; i<nH; i++) 
    sigma_H0[i] = (double)(i+nPd);
*/

  return;
}
