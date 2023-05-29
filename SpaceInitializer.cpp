#include<SpaceInitializer.h>
#include<LatticeVariables.h>
#include<UserDefinedGeometry.h>
#include<Utils.h>
#include<DistancePointToParallelepiped.h>
#include<DistancePointToSpheroid.h>
#include<DistancePointToSphere.h>
#include<DistancePointToCylinderCone.h>
#include<DistancePointToCylinderSphere.h>
#include<DistancePointToPlane.h>
#include<BoundingBoxes.h>
#include<CommunicationTools.h>
#include<cstring>
#include<set>
#include<cmath> //isfinite, floor
#include<climits> //LLONG_MAX
#include<dlfcn.h> //dlopen, dlclose (dynamic linking)
#include<tuple>
using std::isfinite;
using std::vector;
using std::map;
using namespace GeoTools;

//---------------------------------------------------------------------

SpaceInitializer::SpaceInitializer(MPI_Comm& comm_, IoData& iod_, vector<MaterialOperator> &mato_,
                                   vector<LatticeStructure> &lats_) 
                 : iod(iod_), comm(comm_), mato(mato_), lats(lats_)
{

}

//---------------------------------------------------------------------

SpaceInitializer::~SpaceInitializer()
{ }

//---------------------------------------------------------------------

void 
SpaceInitializer::SetupLatticeVariables(vector<LatticeVariables> &LVS)
{

  int nLattices = iod.latticeMap.dataMap.size();
  if((int)LVS.size() != nLattices)
    LVS.resize(nLattices);

  // Create all the lattices
  for(auto&& l : iod.latticeMap.dataMap) {
    int lid = l.first; //lattice id
    assert(lid>=0 && lid<nLattices);
    LVS[lid].lattice_id = lid;
    CreateOneLattice(lats[lid], *l.second, LVS[lid]); 
  }


  // Apply regional IC
  int nRegions = iod.icMap.dataMap.size();
  std::vector<RegionalIcData*> regions(nRegions, NULL);
  std::set<int> region_checker;
  for(auto&& region : iod.icMap.dataMap) {
    int iRegion = region.first; 
    if(iRegion<0 || iRegion>=nRegions) {
      print_error("*** Error: Detected error in RegionalInitialCondition id (%d).\n",
                  iRegion);
      exit_mpi();
    }
    regions[iRegion] = region.second;
    region_checker.insert(iRegion);
  }
  if((int)region_checker.size() != nRegions) {
    print_error("*** Error: Detected error in RegionalInitialCondition indices.\n");
    exit_mpi();
  }
  
  for(int i=0; i<(int)regions.size(); i++) {
    int id = regions[i]->latticeID;
    if(id<0) {// not specified by user --> apply to all lattices
      for(auto&& LV : LVS)
        ApplyRegionalICOnOneLattice(*regions[i], LV); 
    } 
    else if(id<nLattices)
      ApplyRegionalICOnOneLattice(*regions[i], LVS[id]); 
    else {
      print_error("*** Error: Detected unknown lattice ID (%d) in RegionalInitialCondition[%d].\n",
                  id, i);
      exit_mpi();
    }
  }


  // Remove ``conflicting sites'' on different lattices
  for(int iLat=0; iLat<(int)lats.size(); iLat++) {
    if(lats[iLat].GetDmin()<=0.0)
      continue;

    int nRemoved(0);
    for(int jLat=0; jLat<(int)lats.size(); jLat++) {
      if(jLat == iLat)
        continue;
      if(jLat < iLat && lats[jLat].GetDmin()<=lats[iLat].GetDmin()) //already resolved
        continue;
      nRemoved += LVS[iLat].RemoveConflictingSites(LVS[jLat].q, lats[iLat].GetDmin());
    }
    if(nRemoved>0)
      print("- Removed %d sites from Lattice[%d] to maintain dmin = %e from sites on other lattices.\n",
            nRemoved, iLat, lats[iLat].GetDmin());
  }


  print("\n");
  for(int iLat=0; iLat<(int)LVS.size(); iLat++) {
    print("- Created Lattice[%d]: %d sites.\n", iLat, LVS[iLat].size);
  }
}

//---------------------------------------------------------------------

void
SpaceInitializer::CreateOneLattice(LatticeStructure &lat, LatticeData &iod_lat, LatticeVariables &LV)
{

  LV.Clear();

  int lattice_id = lat.GetLatticeID();
  print("- Creating Lattice[%d].\n", lattice_id);


  //----------------------------------------------------------------
  // Step 1: Setup the operation order for user-specified geometries.
  //         Also, create signed distance calculators for later use.
  //----------------------------------------------------------------
  //-------------------------------------
  // Internal Geometry ID:
  // Plane: 0, CylinderCone: 1, CylinderSphere: 2, Sphere: 3,
  // Parallelepiped: 4, Spheroid: 5, Custom-Geometry: 6
  //-------------------------------------
  vector<std::pair<int,int> > order;  //<geom type, geom dataMap index>
  map<int, DistanceFromPointToPlane*> plane_cal;
  map<int, DistanceFromPointToCylinderCone*> cylindercone_cal;
  map<int, DistanceFromPointToCylinderSphere*> cylindersphere_cal;
  map<int, DistanceFromPointToSphere*> sphere_cal;
  map<int, DistanceFromPointToParallelepiped*> parallelepiped_cal;
  map<int, DistanceFromPointToSpheroid*> spheroid_cal;

  SortUserSpecifiedGeometries<LatticeData>(iod_lat, order, plane_cal, cylindercone_cal,
                              cylindersphere_cal, sphere_cal, parallelepiped_cal, spheroid_cal);


  //----------------------------------------------------------------
  // Step 2: Setup user-specified custom geometry specifier
  //----------------------------------------------------------------
  std::tuple<UserDefinedGeometry*,void*,DestroyUDG*> user_geom(std::make_tuple(nullptr,nullptr,nullptr));
  if(strcmp(iod_lat.custom_geometry.specifier,"")) { //user has specified a custom geometry
    std::get<1>(user_geom) = dlopen(iod_lat.custom_geometry.specifier, RTLD_NOW);
    if(std::get<1>(user_geom)) {
      print_error("*** Error: Unable to load object %s.\n", iod_lat.custom_geometry.specifier);
      exit_mpi();
    }
    dlerror();
    CreateUDG* create = (CreateUDG*) dlsym(std::get<1>(user_geom), "Create");
    const char* dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Create in %s.\n",
                  iod_lat.custom_geometry.specifier);
      exit_mpi();
    }
    std::get<2>(user_geom) = (DestroyUDG*) dlsym(std::get<1>(user_geom), "Destroy");
    dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Destroy in %s.\n",
                  iod_lat.custom_geometry.specifier); 
      exit_mpi();
    }

    //This is the actual object. 
    std::get<0>(user_geom) = create();
    print("  o Loaded user-defined geometry specifier from %s.\n", iod_lat.custom_geometry.specifier);
  }


  //----------------------------------------------------------------
  // Step 3: Find a bounding box in terms of lattice coordinates.
  //----------------------------------------------------------------
  Vec3D lmin, lmax;
  FindLatticeDomainBoundingBox(lat, iod_lat, std::get<0>(user_geom), lmin, lmax);

  Int3 lmin_int, lmax_int;
  for(int i=0; i<3; i++) {
    lmin_int[i] = floor(lmin[i]);
    lmax_int[i] = ceil(lmax[i]) + 1;
    //-------------------------------
    //NOTE: lmax_int stores MAX+1, so that in the "for"
    //      loops below, we check k<lmax_int[i] for stopping
    //-------------------------------
  }
    

  //----------------------------------------------------------------
  // Step 4: Check for input error in site group ID
  //----------------------------------------------------------------
  std::set<int> site_group_tracker; //for error detection only
  for(auto&& sg : iod_lat.siteMap.dataMap) {
    if(sg.first<0 || sg.first>=(int)iod_lat.siteMap.dataMap.size()) {
      print_error("*** Error: Lattice[%d] has incorrect site group ID (%d).\n", lattice_id, sg.first);
      exit_mpi();
    }
    site_group_tracker.insert(sg.first);
  }
  if(site_group_tracker.size() != iod_lat.siteMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has incorrect (conflicting) site group IDs.\n", lattice_id);
    exit_mpi();
  }
 

  //----------------------------------------------------------------
  // Step 5: Split the work. Each processor checks some cells.
  //----------------------------------------------------------------
  long long lk_size_ll = lmax_int[2] - lmin_int[2];
  long long lj_size_ll = lmax_int[1] - lmin_int[1];
  long long li_size_ll = lmax_int[0] - lmin_int[0];
  long long sample_size = lk_size_ll*lj_size_ll*li_size_ll;

  if(!isfinite(sample_size)) {
    print_error("*** Error: The bounding box for lattice domain [%d] is too large (>%lld unit cells).\n",
                lattice_id, LLONG_MAX);
    exit_mpi();
  }
  int mpi_rank(-1), mpi_size(0);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);
  
  long long block_size = sample_size/mpi_size; //floor((double)sample_size/(double)mpi_size);
  long long remainder = sample_size - block_size*mpi_size;
  assert(remainder>=0 && remainder<mpi_size);

  //prepare for comm.
  long long* counts = new long long[mpi_size];
  long long* displacements = new long long[mpi_size];
  for(int i=0; i<mpi_size; i++) {
    counts[i] = (i<remainder) ? block_size + 1 : block_size;
    displacements[i] = (i<remainder) ? (block_size+1)*i : block_size*i + remainder;
  }
  assert(displacements[mpi_size-1]+counts[mpi_size-1] == sample_size);

  long long my_start_id = displacements[mpi_rank];
  long long my_block_size = counts[mpi_rank];

  delete [] counts;
  delete [] displacements;


  //----------------------------------------------------------------
  // Step 6: Each processor completes its assignment
  //----------------------------------------------------------------
  int lk0 = my_start_id/(lj_size_ll*li_size_ll);
  int lj0 = (my_start_id - (long long)lk0*lj_size_ll*li_size_ll)/li_size_ll;
  int li0 = my_start_id - (long long)lk0*lj_size_ll*li_size_ll - (long long)lj0*li_size_ll; 
  lk0 += lmin_int[2];
  lj0 += lmin_int[1];
  li0 += lmin_int[0];
  assert(isfinite(lk0) && isfinite(lj0) && isfinite(li0));

  long long it_counter = 0; 
  int nSites = 0;
  for(int lk=lk0; lk<lmax_int[2]; lk++) {
    for(int lj=(lk==lk0 ? lj0 : lmin_int[1]); lj<lmax_int[1]; lj++) {
      for(int li=(lk==lk0 && lj==lj0 ? li0 : lmin_int[0]); li<lmax_int[0]; li++) {

        // add new sites to LV
        nSites += CheckAndAddSitesInCell(li, lj, lk, lat, iod_lat, order, plane_cal, cylindercone_cal,
                                         cylindersphere_cal, sphere_cal, parallelepiped_cal, spheroid_cal,
                                         std::get<0>(user_geom), LV);

        if(++it_counter==my_block_size)
          goto END_OF_ASSIGNMENT; 
      }
    }
  }

END_OF_ASSIGNMENT:


  //----------------------------------------------------------------
  // Step 7: Exchange data between processors
  //----------------------------------------------------------------
  CommunicationTools::AllGatherVector<int>(comm, LV.siteid); //include data type for safety
  CommunicationTools::AllGatherVector<int>(comm, LV.matid);
  CommunicationTools::AllGatherVector<Vec3D>(comm, LV.l);
  CommunicationTools::AllGatherVector<Vec3D>(comm, LV.q);
  CommunicationTools::AllGatherVector<Vec3D>(comm, LV.q0);
  CommunicationTools::AllGatherVector<double>(comm, LV.sigma);
  CommunicationTools::AllGatherVector<int>(comm, LV.tag);

  CommunicationTools::AllGatherVectorOfVectors<int>(comm, LV.diffusive);
  CommunicationTools::AllGatherVectorOfVectors<int>(comm, LV.species_id);
  CommunicationTools::AllGatherVectorOfVectors<double>(comm, LV.x);
  CommunicationTools::AllGatherVectorOfVectors<double>(comm, LV.gamma);
  
  LV.FindBoundingBoxAndSize();
  
  LV.nSpecies_max = lat.GetMaxNumberOfSpeciesPerSite();

  //----------------------------------------------------------------
  // Clean-up
  //----------------------------------------------------------------
  for(auto&& ob : plane_cal)          if(ob.second) delete ob.second;
  for(auto&& ob : cylindercone_cal)   if(ob.second) delete ob.second;
  for(auto&& ob : cylindersphere_cal) if(ob.second) delete ob.second;
  for(auto&& ob : sphere_cal)         if(ob.second) delete ob.second;
  for(auto&& ob : parallelepiped_cal) if(ob.second) delete ob.second;
  for(auto&& ob : spheroid_cal)       if(ob.second) delete ob.second;

  if(std::get<0>(user_geom)) {
    assert(std::get<2>(user_geom)); //this is the destruction function
    (std::get<2>(user_geom))(std::get<0>(user_geom)); //use the destruction function to destroy the calculator
    assert(std::get<1>(user_geom));
    dlclose(std::get<1>(user_geom));
  }
}

//---------------------------------------------------------------------

template<class IoDataGeo>
void
SpaceInitializer::SortUserSpecifiedGeometries(IoDataGeo &iod_geo,
                                              vector<std::pair<int,int> > &order,
                                              map<int, DistanceFromPointToPlane*> &plane_cal,
                                              map<int, DistanceFromPointToCylinderCone*> &cylindercone_cal,
                                              map<int, DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                                              map<int, DistanceFromPointToSphere*> &sphere_cal,
                                              map<int, DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                                              map<int, DistanceFromPointToSpheroid*> &spheroid_cal)
{
  //NOTE: "order" records the index of the geometry in the "dataMap", not the ID of the geometry
  //      specified by the user in the input file.

  //--------------------------------------------------------
  // Step 1: Find the number of user-specified geometries
  //--------------------------------------------------------
  int nGeom = iod_geo.planeMap.dataMap.size() //0
            + iod_geo.cylinderconeMap.dataMap.size() //1
            + iod_geo.cylindersphereMap.dataMap.size() //2
            + iod_geo.sphereMap.dataMap.size() //3
            + iod_geo.parallelepipedMap.dataMap.size() //4
            + iod_geo.spheroidMap.dataMap.size(); //5
  if(strcmp(iod_geo.custom_geometry.specifier, "")) //6 (user has specified a custom geometry)
    nGeom++;

  order.assign(nGeom, std::make_pair(-1,-1));

  //--------------------------------------------------------
  // Step 2: Loop through all of them to determine the order
  //--------------------------------------------------------
  std::set<int> geom_tracker; //for error detection only

  assert(plane_cal.empty());
  for(auto it = iod_geo.planeMap.dataMap.begin(); it != iod_geo.planeMap.dataMap.end(); it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Plane[%d] has invalid or conflicting order (%d).\n", 
                  it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != PlaneData::INTERSECTION && it->second->inclusion != PlaneData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[it->second->order].first  = 0; //planes
    order[it->second->order].second = it->first;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    plane_cal[it->first] = new DistanceFromPointToPlane(x0, dir); 

    print("  o Found Plane[%d]: (%e %e %e) | normal: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
  }
  if(geom_tracker.size() != iod_geo.planeMap.dataMap.size()) {
    print_error("*** Error: Detected invalid (possibly duplicate) planes.\n");
    exit_mpi();
  }
  geom_tracker.clear();


  assert(cylindercone_cal.empty());
  for(auto it = iod_geo.cylinderconeMap.dataMap.begin(); it != iod_geo.cylinderconeMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: CylinderAndCone[%d] has invalid or conflicting order (%d).\n", 
                  it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != CylinderConeData::INTERSECTION && 
       it->second->inclusion != CylinderConeData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[it->second->order].first  = 1; //cylindercone
    order[it->second->order].second = it->first;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = std::min(it->second->cone_height, Hmax); //cone's height
    cylindercone_cal[it->first] = new DistanceFromPointToCylinderCone(x0,dir,L,R,tan_alpha,H);

    print("  o Found Cylinder-Cone[%d]: Base center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
  }
  if(geom_tracker.size() != iod_geo.cylinderconeMap.dataMap.size()) {
    print_error("*** Error: Detected invalid (possibly duplicate) cylinder-cones.\n");
    exit_mpi();
  }
  geom_tracker.clear();


  assert(cylindersphere_cal.empty());
  for(auto it = iod_geo.cylindersphereMap.dataMap.begin(); it != iod_geo.cylindersphereMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: CylinderWidthSphericalCaps[%d] has invalid or "
                  "conflicting order (%d).\n", it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != CylinderSphereData::INTERSECTION && 
       it->second->inclusion != CylinderSphereData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[it->second->order].first  = 2; //cylindersphere
    order[it->second->order].second = it->first;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of base
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);
    cylindersphere_cal[it->first] = new DistanceFromPointToCylinderSphere(x0, dir, L, R, front_cap, back_cap);

    print("  o Found Cylinder-Sphere[%d]: Base center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
  }
  if(geom_tracker.size() != iod_geo.cylindersphereMap.dataMap.size()) {
    print_error("*** Error: Detected invalid (possibly duplicate) cylinder-spheres.\n");
    exit_mpi();
  }
  geom_tracker.clear();


  assert(sphere_cal.empty());
  for(auto it = iod_geo.sphereMap.dataMap.begin(); it != iod_geo.sphereMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Sphere[%d] has invalid or conflicting order (%d).\n",
                  it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != SphereData::INTERSECTION && 
       it->second->inclusion != SphereData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[it->second->order].first  = 3; //sphere
    order[it->second->order].second = it->first;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    double R = it->second->radius;
    sphere_cal[it->first] = new DistanceFromPointToSphere(x0, R);

    print("  o Found Sphere[%d]: Center: (%e %e %e), Radius: %e.\n", it->first,
          x0[0], x0[1], x0[2], R);
  }
  if(geom_tracker.size() != iod_geo.sphereMap.dataMap.size()) {
    print_error("*** Error: Detected invalid (possibly duplicate) spheres.\n");
    exit_mpi();
  }
  geom_tracker.clear();


  assert(parallelepiped_cal.empty());
  for(auto it = iod_geo.parallelepipedMap.dataMap.begin(); it != iod_geo.parallelepipedMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Parallelepiped[%d] has invalid or conflicting order (%d).\n",
                  it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != ParallelepipedData::INTERSECTION && 
       it->second->inclusion != ParallelepipedData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[it->second->order].first  = 4; //parallelepiped
    order[it->second->order].second = it->first;

    Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
    Vec3D oa(it->second->ax, it->second->ay, it->second->az); oa -= x0;
    Vec3D ob(it->second->bx, it->second->by, it->second->bz); ob -= x0;
    Vec3D oc(it->second->cx, it->second->cy, it->second->cz); oc -= x0;

    if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
      print_error("*** Error: Detected error in Parallelepiped[%d]. "
                  "Overlapping vertices or violation of right-hand rule.\n", it->first);
      exit_mpi();
    }

    parallelepiped_cal[it->first] = new DistanceFromPointToParallelepiped(x0, oa, ob, oc);

    print("  o Found Parallelepiped[%d]: X(%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2]);
  }
  if(geom_tracker.size() != iod_geo.parallelepipedMap.dataMap.size()) {
    print_error("*** Error: Detected invalid (possibly duplicate) parallelepipeds.\n");
    exit_mpi();
  }
  geom_tracker.clear();


  assert(spheroid_cal.empty());
  for(auto it = iod_geo.spheroidMap.dataMap.begin(); it != iod_geo.spheroidMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Spheroid[%d] has invalid or "
                  "conflicting order (%d).\n", it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != SpheroidData::INTERSECTION && 
       it->second->inclusion != SpheroidData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[it->second->order].first  = 5; //spheroid
    order[it->second->order].second = it->first;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
    spheroid_cal[it->first] = new DistanceFromPointToSpheroid(x0, axis, it->second->semi_length, it->second->radius);

    print("  o Found Spheroid[%d]: Center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], axis[0], axis[1], axis[2]);
  }
  if(geom_tracker.size() != iod_geo.spheroidMap.dataMap.size()) {
    print_error("*** Error: Detected invalid (possibly duplicate) spheroids.\n");
    exit_mpi();
  }
  geom_tracker.clear();


  if(strcmp(iod_geo.custom_geometry.specifier, "")) {
    if(iod_geo.custom_geometry.order<0 || iod_geo.custom_geometry.order>=nGeom ||
       order[iod_geo.custom_geometry.order].first!=-1) {
      print_error("*** Error: CustomGeometry has invalid or "
                  "conflicting order (%d).\n", iod_geo.custom_geometry.order);
      exit_mpi();
    }

    if(iod_geo.custom_geometry.inclusion != CustomGeometryData::INTERSECTION && 
       iod_geo.custom_geometry.inclusion != CustomGeometryData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[iod_geo.custom_geometry.order].first  = 6; //custom geometry
    order[iod_geo.custom_geometry.order].second = 0;
  }
 
  for(auto&& oo : order)
    assert(oo.first>=0 && oo.second>=0);

}

//---------------------------------------------------------------------

int
SpaceInitializer::CheckAndAddSitesInCell(int li, int lj, int lk, LatticeStructure &lat,
                                         LatticeData &iod_lat, vector<std::pair<int,int> > &order,
                                         map<int,DistanceFromPointToPlane*> &plane_cal,
                                         map<int,DistanceFromPointToCylinderCone*> &cylindercone_cal,
                                         map<int,DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                                         map<int,DistanceFromPointToSphere*> &sphere_cal,
                                         map<int,DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                                         map<int,DistanceFromPointToSpheroid*> &spheroid_cal,
                                         UserDefinedGeometry *geom_specifier,
                                         LatticeVariables &LV)
{

  // Loop through site groups
  Vec3D l, q;
  bool in, intersect;
  int counter = 0;
  int nSites = lat.GetNumberOfSites();
 
  for(int sid = 0; sid < nSites; sid++) {
    l = lat.GetCoordsOfSite(sid) + Vec3D(li, lj, lk); //lattice coords
    q = lat.GetXYZ(l); //(x,y,z) coords

    in = true; //always in R^3...

    for(auto&& oo : order) { //loop through the geometries in user-specified order

      if(oo.first == 0) {//plane
        auto it = iod_lat.planeMap.dataMap[oo.second];
        intersect = it->inclusion == PlaneData::INTERSECTION; //INTERSECTION or UNION
        if((in && intersect) || (!in && !intersect)) //otherwise, no need to check, "in" remains the same
          in = plane_cal[oo.second]->Calculate(q)>0.0;
      }
      else if(oo.first == 1) {//cylindercone
        auto it = iod_lat.cylinderconeMap.dataMap[oo.second];
        intersect = it->inclusion == CylinderConeData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = cylindercone_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==CylinderConeData::EXTERIOR) ||
               (!outside && it->side==CylinderConeData::INTERIOR);
        }
      } 
      else if(oo.first == 2) {//cylindersphere
        auto it = iod_lat.cylindersphereMap.dataMap[oo.second];
        intersect = it->inclusion == CylinderSphereData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = cylindersphere_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==CylinderSphereData::EXTERIOR) ||
               (!outside && it->side==CylinderSphereData::INTERIOR);
        }
      }
      else if(oo.first == 3) {//sphere
        auto it = iod_lat.sphereMap.dataMap[oo.second];
        intersect = it->inclusion == SphereData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = sphere_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==SphereData::EXTERIOR) ||
               (!outside && it->side==SphereData::INTERIOR);
        }
      }
      else if(oo.first == 4) {//parallelepiped
        auto it = iod_lat.parallelepipedMap.dataMap[oo.second];
        intersect = it->inclusion == ParallelepipedData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = parallelepiped_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==ParallelepipedData::EXTERIOR) ||
               (!outside && it->side==ParallelepipedData::INTERIOR);
        }
      }
      else if(oo.first == 5) {//spheroid
        auto it = iod_lat.spheroidMap.dataMap[oo.second];
        intersect = it->inclusion == SpheroidData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = spheroid_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==SpheroidData::EXTERIOR) ||
               (!outside && it->side==SpheroidData::INTERIOR);
        }
      }
      else if(oo.first == 6) {//user-specified custom geometry
        intersect = iod_lat.custom_geometry.inclusion == CustomGeometryData::INTERSECTION;
        if((in && intersect) || (!in && !intersect))
          in = geom_specifier->IsPointInside(q, q, lat.GetLatticeOrigin(), lat.GetLatticeVector(0),
                                             lat.GetLatticeVector(1), lat.GetLatticeVector(2));
      }

    }

    if(in) {//add to LV

      LV.l.push_back(l);
      LV.q0.push_back(q);
      LV.siteid.push_back(sid);

      int matid = lat.GetMatIDOfSite(sid);
      LV.matid.push_back(matid);

      LatticeSiteData &iod_site(*iod_lat.siteMap.dataMap[sid]);
      LV.sigma.push_back(iod_site.atomic_frequency);
      q[0] += iod_site.mean_displacement_x;
      q[1] += iod_site.mean_displacement_y;
      q[2] += iod_site.mean_displacement_z;
      LV.q.push_back(q);

      LV.tag.push_back(iod_site.tag); //TODO: This tag can be generalized

      LV.species_id.push_back(lat.site_species_id[sid]);
      LV.diffusive.push_back(lat.site_species_diff[sid]);
      LV.x.push_back(lat.site_species_x0[sid]);

      LV.gamma.push_back(lat.site_species_gamma[sid]);
      
      LV.size++;
      counter++;

    }
  }

  return counter;

}

//---------------------------------------------------------------------

void
SpaceInitializer::FindLatticeDomainBoundingBox(LatticeStructure &lat, LatticeData &iod_lat,
                                            UserDefinedGeometry *geom_specifier, Vec3D &lmin, Vec3D &lmax)
{

  lmin = DBL_MAX;
  lmax = -DBL_MAX;

  Vec3D O = lat.GetLatticeOrigin();
  Vec3D abc[3];
  for(int i=0; i<3; i++)
    abc[i] = lat.GetLatticeVector(i);

  Vec3D lmin_local, lmax_local;

  // Go over cylinder-cones
  for(auto it = iod_lat.cylinderconeMap.dataMap.begin(); it != iod_lat.cylinderconeMap.dataMap.end(); it++) {
    if(it->second->side == CylinderConeData::EXTERIOR)
      continue;
    
    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of the base disk
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();
    
    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = std::min(it->second->cone_height, Hmax); //cone's height
    
    GetBoundingBoxOfCylinderCone(x0, dir, R, L, tan_alpha, H, O, abc[0], abc[1], abc[2],
                                 lmin_local, lmax_local);

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i]; 
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over cylinder-spheres (i.e. possibly with end cap(s))
  for(auto it = iod_lat.cylindersphereMap.dataMap.begin(); it != iod_lat.cylindersphereMap.dataMap.end(); it++) {
    if(it->second->side == CylinderSphereData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of base disk
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    dir /= dir.norm();

    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder & sphere's radius
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);

    GetBoundingBoxOfCylinderSphere(x0, dir, R, L, front_cap, back_cap, O, abc[0], abc[1],
                                   abc[2], lmin_local, lmax_local);

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i]; 
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over spheres
  for(auto it = iod_lat.sphereMap.dataMap.begin(); it != iod_lat.sphereMap.dataMap.end(); it++) {
    if(it->second->side == SphereData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    double R = it->second->radius;

    GetBoundingBoxOfSphere(x0, R, O, abc[0], abc[1], abc[2],
                           lmin_local, lmax_local);

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
    if(it->second->side == ParallelepipedData::EXTERIOR)
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

    GetBoundingBoxOfParallelepiped(x0, oa, ob, oc, O, abc[0], abc[1], abc[2],
                                   lmin_local, lmax_local);

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i];
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  // Go over spheroids
  for(auto it = iod_lat.spheroidMap.dataMap.begin(); it != iod_lat.spheroidMap.dataMap.end(); it++) {
    if(it->second->side == SpheroidData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);

    GetBoundingBoxOfSpheroid(x0, axis, it->second->semi_length, it->second->radius, O,
                             abc[0], abc[1], abc[2], lmin_local, lmax_local);

    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i];
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }

  //Check user-specified custom geometry
  if(geom_specifier) {
    geom_specifier->GetBoundingBox(O, abc[0], abc[1], abc[2], lmin_local, lmax_local);
    for(int i=0; i<3; i++) {
      if(lmin_local[i] < lmin[i])
        lmin[i] = lmin_local[i];
      if(lmax_local[i] > lmax[i])
        lmax[i] = lmax_local[i];
    }
  }


  for(int i=0; i<3; i++) {
    if(!isfinite(lmin[i]) || !isfinite(lmax[i]) || lmin[i]>=lmax[i]) {
      print_error("*** Error: Unable to find a bounding box for lattice domain [%d]. Domain is unbounded.\n",
                  lat.GetLatticeID());
      exit_mpi();
    }
  }

}

//---------------------------------------------------------------------

void
SpaceInitializer::ApplyRegionalICOnOneLattice(RegionalIcData& ic, LatticeVariables &LV)
{
  if(ic.materialID>=(int)mato.size()) {
    print_error("*** Error: Detected unknown MaterialID (%d) in RegionalInitialCondition.\n",
                ic.materialID);
    exit_mpi();
  }

  //----------------------------------------------------------------
  // Step 1: Setup the operation order for user-specified geometries.
  //         Also, create signed distance calculators for later use.
  //----------------------------------------------------------------
  //-------------------------------------
  // Internal Geometry ID:
  // Plane: 0, CylinderCone: 1, CylinderSphere: 2, Sphere: 3,
  // Parallelepiped: 4, Spheroid: 5, Custom-Geometry: 6
  //-------------------------------------
  vector<std::pair<int,int> > order;  //<geom type, geom dataMap index>
  map<int,DistanceFromPointToPlane*> plane_cal;
  map<int,DistanceFromPointToCylinderCone*> cylindercone_cal;
  map<int,DistanceFromPointToCylinderSphere*> cylindersphere_cal;
  map<int,DistanceFromPointToSphere*> sphere_cal;
  map<int,DistanceFromPointToParallelepiped*> parallelepiped_cal;
  map<int,DistanceFromPointToSpheroid*> spheroid_cal;

  SortUserSpecifiedGeometries<RegionalIcData>(ic, order, plane_cal, cylindercone_cal,
                              cylindersphere_cal, sphere_cal, parallelepiped_cal, spheroid_cal);

  //----------------------------------------------------------------
  // Step 2: Setup user-specified custom geometry specifier
  //----------------------------------------------------------------
  std::tuple<UserDefinedGeometry*,void*,DestroyUDG*> user_geom(std::make_tuple(nullptr,nullptr,nullptr));
  if(strcmp(ic.custom_geometry.specifier,"")) { //user has specified a custom geometry
    std::get<1>(user_geom) = dlopen(ic.custom_geometry.specifier, RTLD_NOW);
    if(std::get<1>(user_geom)) {
      print_error("*** Error: Unable to load object %s.\n", ic.custom_geometry.specifier);
      exit_mpi();
    }
    dlerror();
    CreateUDG* create = (CreateUDG*) dlsym(std::get<1>(user_geom), "Create");
    const char* dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Create in %s.\n",
                  ic.custom_geometry.specifier);
      exit_mpi();
    }
    std::get<2>(user_geom) = (DestroyUDG*) dlsym(std::get<1>(user_geom), "Destroy");
    dlsym_error = dlerror();
    if(dlsym_error) {
      print_error("*** Error: Unable to find function Destroy in %s.\n",
                  ic.custom_geometry.specifier); 
      exit_mpi();
    }

    //This is the actual object. 
    std::get<0>(user_geom) = create();
    print("  o Loaded user-defined geometry specifier from %s.\n", ic.custom_geometry.specifier);
  }
  UserDefinedGeometry *geom_specifier = std::get<0>(user_geom);

  //----------------------------------------------------------------
  // Step 3: Loop through sites
  //----------------------------------------------------------------
  int mid0 = ic.materialID;
  int sid0 = ic.siteID;
  bool in, intersect;

  for(int i=0; i<LV.size; i++) {
    if(LV.matid[i] != mid0 || LV.siteid[i] != sid0)
      continue; //not you. 

    Vec3D& q(LV.q[i]);

    // This part is pretty much the same as CheckAndAddSitesInCell
    
    in = true; //always in R^3...

    for(auto&& oo : order) { //loop through the geometries in user-specified order

      if(oo.first == 0) {//plane
        auto it = ic.planeMap.dataMap[oo.second];
        intersect = it->inclusion == PlaneData::INTERSECTION; //INTERSECTION or UNION
        if((in && intersect) || (!in && !intersect)) //otherwise, no need to check, "in" remains the same
          in = plane_cal[oo.second]->Calculate(q)>0.0;
      }
      else if(oo.first == 1) {//cylindercone
        auto it = ic.cylinderconeMap.dataMap[oo.second];
        intersect = it->inclusion == CylinderConeData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = cylindercone_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==CylinderConeData::EXTERIOR) ||
               (!outside && it->side==CylinderConeData::INTERIOR);
        }
      } 
      else if(oo.first == 2) {//cylindersphere
        auto it = ic.cylindersphereMap.dataMap[oo.second];
        intersect = it->inclusion == CylinderSphereData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = cylindersphere_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==CylinderSphereData::EXTERIOR) ||
               (!outside && it->side==CylinderSphereData::INTERIOR);
        }
      }
      else if(oo.first == 3) {//sphere
        auto it = ic.sphereMap.dataMap[oo.second];
        intersect = it->inclusion == SphereData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = sphere_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==SphereData::EXTERIOR) ||
               (!outside && it->side==SphereData::INTERIOR);
        }
      }
      else if(oo.first == 4) {//parallelepiped
        auto it = ic.parallelepipedMap.dataMap[oo.second];
        intersect = it->inclusion == ParallelepipedData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = parallelepiped_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==ParallelepipedData::EXTERIOR) ||
               (!outside && it->side==ParallelepipedData::INTERIOR);
        }
      }
      else if(oo.first == 5) {//spheroid
        auto it = ic.spheroidMap.dataMap[oo.second];
        intersect = it->inclusion == SpheroidData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = spheroid_cal[oo.second]->Calculate(q)>=0.0;
          in = (outside  && it->side==SpheroidData::EXTERIOR) ||
               (!outside && it->side==SpheroidData::INTERIOR);
        }
      }
      else if(oo.first == 6) {//user-specified custom geometry
        intersect = ic.custom_geometry.inclusion == CustomGeometryData::INTERSECTION;
        if((in && intersect) || (!in && !intersect))
          in = geom_specifier->IsPointInside(q, q, lats[LV.lattice_id].GetLatticeOrigin(), 
                                             lats[LV.lattice_id].GetLatticeVector(0),
                                             lats[LV.lattice_id].GetLatticeVector(1),
                                             lats[LV.lattice_id].GetLatticeVector(2));
      }

    }

    if(in) {//add to LV
      q[0] += ic.mean_displacement_x;
      q[1] += ic.mean_displacement_y;
      q[2] += ic.mean_displacement_z;

      LV.sigma[i] = ic.atomic_frequency;
      LV.tag[i] = ic.tag;

      for(auto&& sp : ic.speciesMap.dataMap) {
        // ignoring sp.first (the temporary ID given by user)
        double xnew = sp.second->molar_fraction;
        if(xnew<0.0) //user did not specify it
          continue;
        int spid = sp.second->speciesID;
        for(int s=0; s<(int)LV.species_id[i].size(); s++) {
          if(LV.species_id[i][s] == spid)
            LV.x[i][s] = xnew;
        }
      }

    }
  }

  //----------------------------------------------------------------
  // Clean-up
  //----------------------------------------------------------------
  for(auto&& ob : plane_cal)          if(ob.second) delete ob.second;
  for(auto&& ob : cylindercone_cal)   if(ob.second) delete ob.second;
  for(auto&& ob : cylindersphere_cal) if(ob.second) delete ob.second;
  for(auto&& ob : sphere_cal)         if(ob.second) delete ob.second;
  for(auto&& ob : parallelepiped_cal) if(ob.second) delete ob.second;
  for(auto&& ob : spheroid_cal)       if(ob.second) delete ob.second;

  if(std::get<0>(user_geom)) {
    assert(std::get<2>(user_geom)); //this is the destruction function
    (std::get<2>(user_geom))(std::get<0>(user_geom)); //use the destruction function to destroy the calculator
    assert(std::get<1>(user_geom));
    dlclose(std::get<1>(user_geom));
  }

}

//---------------------------------------------------------------------







//---------------------------------------------------------------------

