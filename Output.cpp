#include <Output.h>
#include <Utils.h>
#include <cstring>
#include <cassert>
using std::vector;

//--------------------------------------------------------------------------------

Output::Output(MPI_Comm &comm_, IoData& iod_) 
      : comm(comm_), iod(iod_)
{

  iFrame = 0;

  last_snapshot_time = -1.0;


  // store material names
  int nMaterials = iod.materialMap.dataMap.size();
  material_name.resize(nMaterials);
  std::set<int> material_tracker; //for error detection only
  for(auto it = iod.materialMap.dataMap.begin(); it != iod.materialMap.dataMap.end(); it++) {
    int matid = it->first;
    if(matid<0 || matid>=(int)nMaterials) {
      print_error("*** Error: Detected error in the specification of material index (id = %d).\n", matid);
      exit_mpi();
    }
    material_tracker.insert(matid);

    if(strcmp(it->second.name, "")==0) 
      material_name[matid] = "Mat" + std::to_string(matid);
    else
      material_name[matid] = it->second.name;
  }
  if(material_tracker.size() != nMaterials) {
    print_error("*** Error: Detected error in the specification of material indices.\n");
    exit_mpi();
  }


  // store lattice axes and origins
  int nLattices = iod.latticeMap.dataMap.size();
  axis_a.resize(nLattices);
  axis_b.resize(nLattices);
  axis_c.resize(nLattices);
  origin.resize(nLattices);
  periodic_a.resize(nLattices);
  periodic_b.resize(nLattices);
  periodic_c.resize(nLattices);

  std::set<int> lattice_tracker; //for error detection only
  for(auto it = iod.latticeMap.dataMap.begin(); it != iod.latticeMap.dataMap.end(); it++) {
    int lattice_id = it->first;
    if(lattice_id<0 || lattice_id>=(int)nLattices) {
      print_error("*** Error: Detected error in the specification of lattice index (id = %d).\n", lattice_id);
      exit_mpi();
    }
    lattice_tracker.insert(lattice_id);

    axis_a[lattice_id][0] = it->second.ax;
    axis_a[lattice_id][1] = it->second.ay;
    axis_a[lattice_id][2] = it->second.az;

    axis_b[lattice_id][0] = it->second.bx;
    axis_b[lattice_id][1] = it->second.by;
    axis_b[lattice_id][2] = it->second.bz;

    axis_c[lattice_id][0] = it->second.cx;
    axis_c[lattice_id][1] = it->second.cy;
    axis_c[lattice_id][2] = it->second.cz;

    origin[lattice_id][0] = it->second.ox;
    origin[lattice_id][1] = it->second.oy;
    origin[lattice_id][2] = it->second.oz;

    periodic_a[lattice_id] = (it->second.periodic_a == LatticeData::TRUE) ? "T" : "F";
    periodic_b[lattice_id] = (it->second.periodic_b == LatticeData::TRUE) ? "T" : "F";
    periodic_c[lattice_id] = (it->second.periodic_c == LatticeData::TRUE) ? "T" : "F";
  }

  if(lattice_tracker.size() != nLattices) {
    print_error("*** Error: Detected error in the specification of lattice indices.\n");
    exit_mpi();
  }


  // setup pvd file if needed
  if(iod.output.format == Output::VTP || iod.output.format == Output::VTP_and_XYZ) {
    for(int lat=0; lat<(int)VS.size(); lat++) {
      char f1[256];
      sprintf(f1, "%s%s_L%d.pvd", iod.output.prefix, iod.output.solution_filename_base, lat);

      FILE* pvdfile = fopen(f1, "w");
      if(!pvdfile) {
        print_error("*** Error: Cannot open file '%s' for output.\n", f1);
        exit_mpi();
      }   

      print(pvdfile, "<?xml version=\"1.0\"?>\n");
      print(pvdfile, "<VTKFile type=\"Collection\" version=\"0.1\"\n");
      print(pvdfile, "byte_order=\"LittleEndian\">\n");
      print(pvdfile, "  <Collection>\n");

      print(pvdfile, "  </Collection>\n");
      print(pvdfile, "</VTKFile>\n");

      fclose(pvdfile);
      delete pvdfile;
    }
  }

}

//--------------------------------------------------------------------------------

Output::~Output()
{
}

//--------------------------------------------------------------------------------

void
Output::OutputSolution(double time, double dt, int time_step, LatticeVariables &V, bool force_write)
{

  if(isTimeToWrite(time, dt, time_step, iod.output.frequency_dt, iod.output.frequency,
                   last_snapshot_time, force_write)) {
    if(iod.output.format == Output::VTP || iod.output.format == Output::VTP_and_XYZ)
      OutputSolutionVTP(time, time_step, V);
    else if(iod.output.format == Output::XYZ || iod.output.format == Output::VTP_and_XYZ)
      OutputSolutionXYZ(time, time_step, V);

    iFrame++;
    last_snapshot_time = time;
  }

}

//--------------------------------------------------------------------------------

void
Output::OutputSolutionVTP(double time, [[maybe_unused]] int time_step, vector<LatticeVariables> &VS) 
{
  for(int lat=0; lat<(int)VS.size(); lat++)
    OutputSolutionOnLatticeVTP(time, lat, VS[lat]);
}
  
//--------------------------------------------------------------------------------

void
Output::OutputSolutionOnLatticeVTP(double time, int lattice_id, LatticeVariables &V) 
{

  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);
  if(MPI_rank)
    return;

  // create solution file
  char full_fname[256];
  char fname[256];
  if(iFrame<10) {
    sprintf(fname, "%s_L%d_000%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_000%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }
  else if(iFrame<100){
    sprintf(fname, "%s_L%d_00%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_00%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }
  else if(iFrame<1000){
    sprintf(fname, "%s_L%d_0%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_0%d.vtr", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }
  else {
    sprintf(fname, "%s_L%d_%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_%d.vtp", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }


  FILE *vtpfile = fopen(full_fname, "w");
  assert(vtpfile);

  // Write header (only proc #0 will write)
  print(vtpfile, "<?xml version=\"1.0\"?>\n");
  print(vtpfile, "<VTKFile type=\"PolyData\" version=\"0.1\"\n");
  print(vtpfile, "byte_order=\"LittleEndian\">\n");
  print(vtpfile, "  <PolyData>\n");
  print(vtpfile, "    <Piece NumberOfPoints=\"%ld\">\n", V.q.size());

  // Write coords ("q")
  print(vtpfile, "      <Points>\n");
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(auto&& q : V.q)
    print(vtpfile, "          %16.10e  %16.10e  %16.10e\n", q[0], q[1], q[2]);
  print(vtpfile, "        </DataArray>\n");
  print(vtpfile, "      </Points>\n");
  
  print(vtpfile, "      <PointData Scalars=\"ScalarResults\" Vectors=\"VectorResults\">\n");

  // Write site id 
  print(vtpfile, "        <DataArray type=\"Int8\" NumberOfComponents=\"1\" Name=\"SiteID\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(V.siteid.size() == V.q.size());
  for(auto&& sid : V.siteid)
    print(vtpfile, "%d ", sid);
  print(vtpfile, "        </DataArray>\n");

  // Write material id
  print(vtpfile, "        <DataArray type=\"Int8\" NumberOfComponents=\"1\" Name=\"MaterialID\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(V.matid.size() == V.q.size());
  for(auto&& id : V.matid)
    print(vtpfile, "%d ", id);
  print(vtpfile, "        </DataArray>\n");

  // Write location (``subsurface'') id
  print(vtpfile, "        <DataArray type=\"Int8\" NumberOfComponents=\"1\" Name=\"SubSurface\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(V.subsurf.size() == V.q.size());
  for(auto&& id : V.subsurf)
    print(vtpfile, "%d ", id);
  print(vtpfile, "        </DataArray>\n");

  // Write atomic frequency (sigma)
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"AtomicFrequency\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(V.sigma.size() == V.q.size());
  for(auto&& sig : V.sigma)
    print(vtpfile, "%e ", sigma);
  print(vtpfile, "        </DataArray>\n");

  
  int ns = 0; //number of species (max)
  for(auto&& xloc : V.x) {
    if(ns<(int)xloc.size())
      ns = xloc.size();
  }
  assert(ns>0);

  // Write molar fraction
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"%d\" Name=\"MolarFraction\" format=\"ascii\">\n", ns);
  print(vtpfile, "          ");
  assert(V.x.size() == V.q.size());
  for(auto&& x : V.x) {
    int i=0; 
    while(i<(int)x.size()) 
      print(vtpfile, "%e ", x[i++]);
    while(i++<ns)
      print(vtpfile, "0 ");
  } 
  print(vtpfile, "        </DataArray>\n");

  // Write chemical potential 
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"%d\" Name=\"ChemicalPotential\" format=\"ascii\">\n", ns);
  print(vtpfile, "          ");
  assert(V.gamma.size() == V.q.size());
  for(auto&& gam : V.gamma) {
    int i=0; 
    while(i<(int)gam.size()) 
      print(vtpfile, "%e ", gam[i++]);
    while(i++<ns)
      print(vtpfile, "0 ");
  } 
  print(vtpfile, "        </DataArray>\n");

  print(vtpfile, "      </PointData>\n");
  print(vtpfile, "    </Piece>\n");
  print(vtpfile, "  </PolyData>\n");
  print(vtpfile, "</VTKFile>\n");
  fclose(vtpfile);
  delete vtpfile;
  

  // Add a line to the pvd file to record the new solutio snapshot
  char f1[256];
  sprintf(f1, "%s%s_L%d.pvd", iod.output.prefix, iod.output.solution_filename_base, lattice_id);
  FILE* pvdfile  = fopen(f1,"r+");
  if(!pvdfile) {
    print_error("*** Error: Cannot open file '%s' for output.\n", f1);
    exit_mpi();
  }
  fseek(pvdfile, -27, SEEK_END); //!< overwrite the previous end of file script
  print(pvdfile, "  <DataSet timestep=\"%e\" file=\"%s\"/>\n", time, fname);
  print(pvdfile, "  </Collection>\n");
  print(pvdfile, "</VTKFile>\n");
  fclose(pvdfile);
  delete pvdfile;

  print("- Wrote solution at %e to %s.\n", time, fname);

}

//--------------------------------------------------------------------------------

void
Output::OutputSolutionXYZ(double time, [[maybe_unused]] int time_step, vector<LatticeVariables> &VS) 
{
  for(int lat=0; lat<(int)VS.size(); lat++)
    OutputSolutionOnLatticeXYZ(time, lat, VS[lat]);
}

//--------------------------------------------------------------------------------

void
Output::OutputSolutionOnLatticeXYZ(double time, int lattice_id, LatticeVariables &V) 
{
  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);
  if(MPI_rank)
    return;

  // sanity checks
  auto np = V.q.size();
  assert(V.siteid.size() == np);
  assert(V.matid.size() == np);
  assert(V.sigma.size() == np);
  assert(V.x.size() == np);
  assert(V.gamma.size() == np);
  assert(V.subsurf.size() == np);
  assert(lattice_id>=0 && lattice_id<(int)axis_a.size());


  // create solution file
  char full_fname[256];
  if(iFrame<10)
    sprintf(full_fname, "%s%s_L%d_000%d.xyz", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  else if(iFrame<100)
    sprintf(full_fname, "%s%s_L%d_00%d.xyz", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  else if(iFrame<1000)
    sprintf(full_fname, "%s%s_L%d_0%d.xyz", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  else
    sprintf(full_fname, "%s%s_L%d_%d.xyz", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);

  FILE *xyzfile = fopen(full_fname, "w");
  assert(xyzfile);


  //Write header
  print(xyzfile, "%d\n", np);
  print(xyzfile, "Lattice=\"%e %e %e %e %e %e %e %e %e\" ", 
        axis_a[lattice_id][0], axis_a[lattice_id][1], axis_a[lattice_id][2],
        axis_b[lattice_id][0], axis_b[lattice_id][1], axis_b[lattice_id][2],
        axis_c[lattice_id][0], axis_c[lattice_id][1], axis_c[lattice_id][2]);
  print(xyzfile, "Origin=\"%e %e %e\" ",
        origin[lattice_id][0], origin[lattice_id][1], origin[lattice_id][2]);
  print(xyzfile, "pbc=\"%s %s %s\" ",
        periodic_a[lattice_id].c_str(), periodic_b[lattice_id].c_str(), periodic_c[lattice_id].c_str());

  print(xyzfile, "Properties=\"species:S:1:"); //material name
  print(xyzfile, "pos:R:3:"); //position (q)
  print(xyzfile, "disp:R:3:"); //displacement (q-q0)
  print(xyzfile, "type:I:1:"); //site id
  print(xyzfile, "atomic_freq:R:1:"); //atomic frequency
  print(xyzfile, "subsurf:I:1:"); //subsurf 

  int ns = 0; //number of species (max)
  for(auto&& xloc : V.x) {
    if(ns<(int)xloc.size())
      ns = xloc.size();
  }
  assert(ns>0);

  print(xyzfile, "molar_frac:R:%d:", ns);
  print(xyzfile, "chem_potent:R:%d\" ", ns);

  print(xyzfile, "Time=%e\n", time);


  //Write results
  for(int i=0; i<np; i++) {
    print(xyzfile, "%16s  ", material_name[V.matid[i]].c_str());
    print(xyzfile, "%16.10e  %16.10e  %16.10e  ", V.q[i][0], V.q[i][1], V.q[i][2]);
    print(xyzfile, "%16.10e  %16.10e  %16.10e  ", 
          V.q[i][0]-V.q0[i][0], V.q[i][1]-V.q0[i][1], V.q[i][2]-V.q0[i][2]);
    print(xyzfile, "%4d  %16.10e  %4d  ", V.siteid[i], V.sigma[i], V.subsurf[i]);

    int s = 0; 
    while(s<(int)V.x.size()) 
      print(xyzfile, "%16.10e  ", V.x[s++]);
    while(s++<ns)
      print(xyzfile, "%16.10e  ", 0.0);

    s = 0;
    while(s<(int)V.gamma.size()) 
      print(xyzfile, "%16.10e  ", V.gamma[s++]);
    while(s++<ns)
      print(xyzfile, "%16.10e  ", 0.0);

    print(xyzfile,"\n");
  }

  fclose(xyzfile);
  delete xyzfile;
}

//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------

