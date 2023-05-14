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
#include<cmath> //isfinite, floor
#include<climits> //LLONG_MAX
#include<dlfcn.h> //dlopen, dlclose (dynamic linking)
#include<tuple>
using std::isfinite;
using std::vector;
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
  if(LVS.size() != nLattices)
    LVS.resize(nLattices);

  for(auto&& l : iod.latticeMap.dataMap) {
    int lid = l.first; //lattice id
    assert(lid>=0 && lid<nLattices);
    CreateOneLattice(LVS[lid], lats[lid], *l.second); 
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
  vector<DistanceFromPointToPlane*> plane_cal;
  vector<DistanceFromPointToCylinderCone*> cylindercone_cal;
  vector<DistanceFromPointToCylinderSphere*> cylindersphere_cal;
  vector<DistanceFromPointToSphere*> sphere_cal;
  vector<DistanceFromPointToParallelepiped*> parallelepiped_cal;
  vector<DistanceFromPointToSpheroid*> spheroid_cal;

  SortUserSpecifiedGeometries(lattice_id, iod_lat, order, plane_cal, cylindercone_cal,
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
    if(sg.first<0 || sg.first>=iod_lat.siteMap.dataMap.size()) {
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


  //----------------------------------------------------------------
  // Step 6: Each processor completes its assignment
  //----------------------------------------------------------------
  int lk0 = my_start_id/(lj_size_ll*li_size_ll);
  int lj0 = (my_start_id - (long long)lk0*lj_size_ll*li_size_ll)/li_size_ll;
  int li0 = my_start_id - (long long)lk0*lj_size_ll*li_size_ll - (long long)ij0*li_size_ll; 
  lk0 += lmin_int[2];
  lj0 += lmin_int[1];
  li0 += lmin_int[0];
  assert(isfinite(lk0) && isfinite(lk1) && isfinite(lk2));

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
  LV.size = LV.q.size();
  

  //----------------------------------------------------------------
  // Clean-up
  //----------------------------------------------------------------
  for(auto&& ob : plane_cal)          if(ob) delete ob;
  for(auto&& ob : cylindercone_cal)   if(ob) delete ob;
  for(auto&& ob : cylindersphere_cal) if(ob) delete ob;
  for(auto&& ob : sphere_cal)         if(ob) delete ob;
  for(auto&& ob : parallelepiped_cal) if(ob) delete ob;
  for(auto&& ob : spheroid_cal)       if(ob) delete ob;

  if(std::get<0>(user_geom)) {
    assert(std::get<2>(user_geom)); //this is the destruction function
    (std::get<2>(user_geom))(std::get<0>(user_geom)); //use the destruction function to destroy the calculator
    assert(std::get<1>(user_geom));
    dlclose(std::get<1>(user_geom));
  }
}

//---------------------------------------------------------------------

void
SpaceInitializer::SortUserSpecifiedGeometries(int lattice_id, LatticeData &iod_lat,
                                           vector<std::pair<int,int> > &order,
                                           vector<DistanceFromPointToPlane*> &plane_cal,
                                           vector<DistanceFromPointToCylinderCone*> &cylindercone_cal,
                                           vector<DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                                           vector<DistanceFromPointToSphere*> &sphere_cal,
                                           vector<DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                                           vector<DistanceFromPointToSpheroid*> spheroid_cal)
{
  //NOTE: "order" records the index of the geometry in the "dataMap", not the ID of the geometry
  //      specified by the user in the input file.

  //--------------------------------------------------------
  // Step 1: Find the number of user-specified geometries
  //--------------------------------------------------------
  int nGeom = iod_lat.planeMap.dataMap.size() //0
            + iod_lat.cylinderconeMap.dataMap.size() //1
            + iod_lat.cylindersphereMap.dataMap.size() //2
            + iod_lat.sphereMap.dataMap.dataMap.size() //3
            + iod_lat.parallelepipedMap.dataMap.size() //4
            + iod_lat.spheroidMap.dataMap.size(); //5
  if(!strcmp(iod_lat.custom_geometry.specifier, "")) //6 (user has specified a custom geometry)
    nGeom++;

  order.assign(nGeom, std::make_pair(-1,-1));

  //--------------------------------------------------------
  // Step 2: Loop through all of them to determine the order
  //--------------------------------------------------------
  int gid;
  std::set<int> geom_tracker; //for error detection only

  assert(plane_cal.empty());
  plane_cal.assign(iod_lat.planeMap.dataMap.size(), NULL);
  for(auto it = iod_lat.planeMap.dataMap.begin(); it != iod_lat.planeMap.dataMap.end(); it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Lattice[%d]->Plane[%d] has invalid or conflicting order (%d).\n", 
                  lattice_id, it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != PlaneData::INTERSECTION && it->second->inclusion != PlaneData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    gid = it - iod_lat.planeMap.dataMap.begin(); 
    order[it->second->order].first  = 0; //planes
    order[it->second->order].second = gid;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    plane_cal[gid] = new DistanceFromPointToPlane(x0, dir); 

    print("  o Found Plane[%d]: (%e %e %e) | normal: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
  }
  if(geom_tracker.size() != iod_lat.planeMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has invalid (possibly duplicate) planes.\n", lattice_id);
    exit_mpi();
  }
  geom_tracker.clear();


  assert(cylindercone_cal.empty());
  cylindercone_cal.assign(iod_lat.cylinderconeMap.dataMap.size(), NULL);
  for(auto it = iod_lat.cylinderconeMap.dataMap.begin(); it != iod_lat.cylinderconeMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Lattice[%d]->CylinderAndCone[%d] has invalid or conflicting order (%d).\n", 
                  lattice_id, it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != CylinderConeData::INTERSECTION && 
       it->second->inclusion != CylinderConeData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    gid = it - iod_lat.cylinderconeMap.dataMap.begin();
    order[it->second->order].first  = 1; //cylindercone
    order[it->second->order].second = gid;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    double tan_alpha = tan(it->second->opening_angle_degrees/180.0*acos(-1.0));//opening angle
    double Hmax = R/tan_alpha;
    double H = min(it->second->cone_height, Hmax); //cone's height
    cylindercone_cal[gid] = new DistanceFromPointToCylinderCone(x0,dir,L,R,tan_alpha,H);

    print("  o Found Cylinder-Cone[%d]: Base center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
  }
  if(geom_tracker.size() != iod_lat.cylinderconeMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has invalid (possibly duplicate) cylinder-cones.\n", lattice_id);
    exit_mpi();
  }
  geom_tracker.clear();


  assert(cylindersphere_cal.empty());
  cylindersphere_cal.assign(iod_lat.cylindersphereMap.dataMap.size(), NULL);
  for(auto it = iod_lat.cylindersphereMap.dataMap.begin(); it != iod_lat.cylindersphereMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Lattice[%d]->CylinderWidthSphericalCaps[%d] has invalid or "
                  "conflicting order (%d).\n", lattice_id, it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != CylinderSphereData::INTERSECTION && 
       it->second->inclusion != CylinderSphereData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    gid = it - iod_lat.cylindersphereMap.dataMap.begin();
    order[it->second->order].first  = 2; //cylindersphere
    order[it->second->order].second = gid;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z); //center of base
    Vec3D dir(it->second->nx, it->second->ny, it->second->nz);
    double L = it->second->L; //cylinder height
    double R = it->second->r; //cylinder radius
    bool front_cap = (it->second->front_cap == CylinderSphereData::On);
    bool back_cap = (it->second->back_cap == CylinderSphereData::On);
    cylindersphere_cal[gid] = new DistanceFromPointToCylinderSphere(x0, dir, L, R, front_cap, back_cap);

    print("  o Found Cylinder-Sphere[%d]: Base center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], dir[0], dir[1], dir[2]);
  }
  if(geom_tracker.size() != iod_lat.cylindersphereMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has invalid (possibly duplicate) cylinder-spheres.\n", lattice_id);
    exit_mpi();
  }
  geom_tracker.clear();


  assert(sphere_cal.empty());
  sphere_cal.assign(iod_lat.sphereMap.dataMap.size(), NULL);
  for(auto it = iod_lat.sphereMap.dataMap.begin(); it != iod_lat.sphereMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Lattice[%d]->Sphere[%d] has invalid or "
                  "conflicting order (%d).\n", lattice_id, it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != SphereData::INTERSECTION && 
       it->second->inclusion != SphereData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    gid = it - iod_lat.sphereMap.dataMap.begin();
    order[it->second->order].first  = 3; //sphere
    order[it->second->order].second = gid;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    double R = it->second->radius;
    sphere_cal[gid] = new DistanceFromPointToSphere(x0, R);

    print("  o Found Sphere[%d]: Center: (%e %e %e), Radius: %e.\n", it->first,
          x0[0], x0[1], x0[2], R);
  }
  if(geom_tracker.size() != iod_lat.sphereMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has invalid (possibly duplicate) spheres.\n", lattice_id);
    exit_mpi();
  }
  geom_tracker.clear();


  assert(parallelepiped_cal.empty());
  parallelepiped_cal.assign(iod_lat.parallelepipedMap.dataMap.size(), NULL);
  for(auto it = iod_lat.parallelepipedMap.dataMap.begin(); it != iod_lat.parallelepipedMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Lattice[%d]->Parallelepiped[%d] has invalid or "
                  "conflicting order (%d).\n", lattice_id, it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != ParallelepipedData::INTERSECTION && 
       it->second->inclusion != ParallelepipedData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    gid = it - iod_lat.parallelepipedMap.dataMap.begin();
    order[it->second->order].first  = 4; //parallelepiped
    order[it->second->order].second = gid;

    Vec3D x0(it->second->x0, it->second->y0, it->second->z0);
    Vec3D oa(it->second->ax, it->second->ay, it->second->az); oa -= x0;
    Vec3D ob(it->second->bx, it->second->by, it->second->bz); ob -= x0;
    Vec3D oc(it->second->cx, it->second->cy, it->second->cz); oc -= x0;

    if(oa.norm()==0 || ob.norm()==0 || oc.norm()==0 || (oa^ob)*oc<=0.0) {
      print_error("*** Error: Detected error in Lattice[%d]->Parallelepiped[%d]. "
                  "Overlapping vertices or violation of right-hand rule.\n", lattice_id, it->first);
      exit_mpi();
    }

    parallelepiped_cal[gid] = new DistanceFromPointToParallelepiped(x0, oa, ob, oc);

    print("  o Found Parallelepiped[%d]: X(%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2]);
  }
  if(geom_tracker.size() != iod_lat.parallelepipedMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has invalid (possibly duplicate) parallelepipeds.\n", lattice_id);
    exit_mpi();
  }
  geom_tracker.clear();


  assert(spheroid_cal.empty());
  spheroid_cal.assign(iod_lat.spheroidMap.dataMap.size(), NULL);
  for(auto it = iod_lat.spheroidMap.dataMap.begin(); it != iod_lat.spheroidMap.dataMap.end();
      it++) {
    geom_tracker.insert(it->first);
    if(it->second->order<0 || it->second->order>=nGeom || order[it->second->order].first!=-1) {
      print_error("*** Error: Lattice[%d]->Spheroid[%d] has invalid or "
                  "conflicting order (%d).\n", lattice_id, it->first, it->second->order);
      exit_mpi();
    }

    if(it->second->inclusion != SpheroidData::INTERSECTION && 
       it->second->inclusion != SpheroidData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    gid = it - iod_lat.spheroidMap.dataMap.begin();
    order[it->second->order].first  = 5; //spheroid
    order[it->second->order].second = gid;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
    spheroid_cal[gid] = new DistanceFromPointToSpheroid(x0, axis, it->second->semi_length, it->second->radius);

    print("  o Found Spheroid[%d]: Center: (%e %e %e), Axis: (%e %e %e).\n", it->first,
          x0[0], x0[1], x0[2], axis[0], axis[1], axis[2]);
  }
  if(geom_tracker.size() != iod_lat.spheroidMap.dataMap.size()) {
    print_error("*** Error: Lattice[%d] has invalid (possibly duplicate) spheroids.\n", lattice_id);
    exit_mpi();
  }
  geom_tracker.clear();


  if(!strcmp(iod_lat.custom_geometry.specifier, "")) {
    if(iod_lat.custom_geometry.order<0 || iod_lat.custom_geometry>=nGeom ||
       order[iod_lat.custom_geometry.order].first!=-1) {
      print_error("*** Error: Lattice[%d]->CustomGeometry has invalid or "
                  "conflicting order (%d).\n", lattice_id, iod_lat.custom_geometry.order);
      exit_mpi();
    }

    if(iod_lat.custom_geometry.inclusion != CustomGeometryData::INTERSECTION && 
       iod_lat.custom_geometry.inclusion != CustomGeometryData::UNION) {
      print_error("*** Error: Detected unsupported geometry inclusion operator "
                  "(must be Intersection or Union).\n");
      exit_mpi();
    }

    order[iod_lat.custom_geometry.order].first  = 6; //custom geometry
    order[iod_lat.custom_geometry.order].second = 0;
  }
 
  for(auto&& oo : order)
    assert(oo.first>=0 && oo.second>=0);

}

//---------------------------------------------------------------------

int
SpaceInitializer::CheckAndAddSitesInCell(int li, int lj, int lk, LatticeStructure &lat,
                                      LatticeData &iod_lat, vector<std::pair<int,int> > &order,
                                      vector<DistanceFromPointToPlane*> &plane_cal,
                                      vector<DistanceFromPointToCylinderCone*> &cylindercone_cal,
                                      vector<DistanceFromPointToCylinderSphere*> &cylindersphere_cal,
                                      vector<DistanceFromPointToSphere*> &sphere_cal,
                                      vector<DistanceFromPointToParallelepiped*> &parallelepiped_cal,
                                      vector<DistanceFromPointToSpheroid*> spheroid_cal,
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
        auto it = iod_lat.planeMap.dataMap.begin() + oo.second;
        intersect = it->second->inclusion == PlaneData::INTERSECTION; //INTERSECTION or UNION
        if((in && intersect) || (!in && !intersect)) //otherwise, no need to check, "in" remains the same
          in = plane_cal[oo.second]->Calculate(q)>0.0;
      }
      else if(oo.first == 1) {//cylindercone
        auto it = iod_lat.cylinderconeMap.dataMap.begin() + oo.second;
        intersect = it->second->inclusion == CylinderConeData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = cylindercone_cal[oo.second]->Calculate(q)>0.0;
          in = (outside  && it->second->side==CylinderConeData::EXTERIOR) ||
               (!outside && it->second->side==CylinderConeData::INTERIOR)
        }
      } 
      else if(oo.first == 2) {//cylindersphere
        auto it = iod_lat.cylindersphereMap.dataMap.begin() + oo.second;
        intersect = it->second->inclusion == CylinderSphereData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = cylindersphere_cal[oo.second]->Calculate(q)>0.0;
          in = (outside  && it->second->side==CylinderSphereData::EXTERIOR) ||
               (!outside && it->second->side==CylinderSphereData::INTERIOR)
        }
      }
      else if(oo.first == 3) {//sphere
        auto it = iod_lat.sphereMap.dataMap.begin() + oo.second;
        intersect = it->second->inclusion == SphereData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = sphere_cal[oo.second]->Calculate(q)>0.0;
          in = (outside  && it->second->side==SphereData::EXTERIOR) ||
               (!outside && it->second->side==SphereData::INTERIOR)
        }
      }
      else if(oo.first == 4) {//parallelepiped
        auto it = iod_lat.parallelepipedMap.dataMap.begin() + oo.second;
        intersect = it->second->inclusion == ParallelepipedData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = parallelepiped_cal[oo.second]->Calculate(q)>0.0;
          in = (outside  && it->second->side==ParallelepipedData::EXTERIOR) ||
               (!outside && it->second->side==ParallelepipedData::INTERIOR)
        }
      }
      else if(oo.first == 5) {//spheroid
        auto it = iod_lat.spheroidMap.dataMap.begin() + oo.second;
        intersect = it->second->inclusion == SpheroidData::INTERSECTION;
        if((in && intersect) || (!in && !intersect)) {
          bool outside = spheroid_cal[oo.second]->Calculate(q)>0.0;
          in = (outside  && it->second->side==SpheroidData::EXTERIOR) ||
               (!outside && it->second->side==SpheroidData::INTERIOR)
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
      LV.matid.push_back(lat.GetMatIDOfSite(sid));

      LatticeSiteData &iod_site(*(iod_lattice.siteMap.dataMap.begin()+lat.GetDataMapIDOfSite(sid))->second);
      LV.sigma.push_back(iod_site.atomic_frequency);
      q[0] += iod_site.mean_displacement_x;
      q[1] += iod_site.mean_displacement_y;
      q[2] += iod_site.mean_displacement_z;
      LV.q.push_back(q);

      I AM HERE! 



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
    if(it->second->inclusion == CylinderConeData::EXTERIOR)
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
    if(it->second->inclusion == CylinderSphereData::EXTERIOR)
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
    if(it->second->inclusion == SphereData::EXTERIOR)
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
    if(it->second->inclusion == SpheroidData::EXTERIOR)
      continue;

    Vec3D x0(it->second->cen_x, it->second->cen_y, it->second->cen_z);
    Vec3D axis(it->second->axis_x, it->second->axis_y, it->second->axis_z);
    double semi_length = it->second->semi_length;
    double r = it->second->radius;

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
                  lattice_id);
      exit_mpi();
    }
  }

}

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
