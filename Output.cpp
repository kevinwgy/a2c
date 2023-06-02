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

  precision = 2;
  if(iod.output.precision == OutputData::LOW)  
    precision = 0;
  else if(iod.output.precision == OutputData::MEDIUM)
    precision = 1;


  // store material names
  int nMaterials = iod.materialMap.dataMap.size();
  material_name.resize(nMaterials);
  for(auto it = iod.materialMap.dataMap.begin(); it != iod.materialMap.dataMap.end(); it++) {
    int matid = it->first;
    if(strcmp(it->second->name, "")==0) //user did not specify a name
      material_name[matid] = "Mat" + std::to_string(matid);
    else
      material_name[matid] = it->second->name;
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

  for(auto it = iod.latticeMap.dataMap.begin(); it != iod.latticeMap.dataMap.end(); it++) {
    int lattice_id = it->first;

    axis_a[lattice_id][0] = it->second->ax;
    axis_a[lattice_id][1] = it->second->ay;
    axis_a[lattice_id][2] = it->second->az;

    axis_b[lattice_id][0] = it->second->bx;
    axis_b[lattice_id][1] = it->second->by;
    axis_b[lattice_id][2] = it->second->bz;

    axis_c[lattice_id][0] = it->second->cx;
    axis_c[lattice_id][1] = it->second->cy;
    axis_c[lattice_id][2] = it->second->cz;

    origin[lattice_id][0] = it->second->ox;
    origin[lattice_id][1] = it->second->oy;
    origin[lattice_id][2] = it->second->oz;

    periodic_a[lattice_id] = (it->second->a_periodic == LatticeData::TRUE) ? "T" : "F";
    periodic_b[lattice_id] = (it->second->b_periodic == LatticeData::TRUE) ? "T" : "F";
    periodic_c[lattice_id] = (it->second->c_periodic == LatticeData::TRUE) ? "T" : "F";
  }


  // write the header of the pvd file if needed
  if(iod.output.format == OutputData::VTP || iod.output.format == OutputData::VTP_and_XYZ) {
    for(int lat=0; lat<(int)iod.latticeMap.dataMap.size(); lat++) {
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
    }
  }

  MPI_Barrier(comm);
}

//--------------------------------------------------------------------------------

Output::~Output()
{
}

//--------------------------------------------------------------------------------

void
Output::OutputSolution(double time, double dt, int time_step, vector<LatticeVariables> &LVS,
                       bool force_write)
{
  if(isTimeToWrite(time, dt, time_step, iod.output.frequency_dt, iod.output.frequency,
                   last_snapshot_time, force_write)) {
    if(iod.output.format == OutputData::VTP || iod.output.format == OutputData::VTP_and_XYZ)
      OutputSolutionVTP(time, time_step, LVS);
    if(iod.output.format == OutputData::XYZ || iod.output.format == OutputData::VTP_and_XYZ)
      OutputSolutionXYZ(time, time_step, LVS);

    iFrame++;
    last_snapshot_time = time;
  }
}

//--------------------------------------------------------------------------------

void
Output::OutputSolutionVTP(double time, [[maybe_unused]] int time_step, vector<LatticeVariables> &LVS) 
{
  for(int lat=0; lat<(int)LVS.size(); lat++)
    OutputSolutionOnLatticeVTP(time, lat, LVS[lat]);
}
  
//--------------------------------------------------------------------------------

void
Output::OutputSolutionOnLatticeVTP(double time, int lattice_id, LatticeVariables &LV) 
{

  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);

  // create solution file
  char full_fname[256];
  char fname[256];
  if(iFrame<10) {
    sprintf(fname, "%s_L%d_000%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_000%d.vtp", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }
  else if(iFrame<100){
    sprintf(fname, "%s_L%d_00%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_00%d.vtp", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }
  else if(iFrame<1000){
    sprintf(fname, "%s_L%d_0%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_0%d.vtp", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }
  else {
    sprintf(fname, "%s_L%d_%d.vtp",
            iod.output.solution_filename_base, lattice_id, iFrame);
    sprintf(full_fname, "%s%s_L%d_%d.vtp", iod.output.prefix,
            iod.output.solution_filename_base, lattice_id, iFrame);
  }


  FILE *vtpfile = fopen(full_fname, "w");
  if(!vtpfile) {
    print_error("*** Error: Cannot open file '%s' for output.\n", full_fname);
    exit_mpi();
  }   

  // Write header (only proc #0 will write)
  print(vtpfile, "<?xml version=\"1.0\"?>\n");
  print(vtpfile, "<VTKFile type=\"PolyData\" version=\"0.1\"\n");
  print(vtpfile, "byte_order=\"LittleEndian\">\n");
  print(vtpfile, "  <PolyData>\n");
  print(vtpfile, "    <Piece NumberOfPoints=\"%ld\">\n", LV.q.size());

  // Write coords ("q")
  print(vtpfile, "      <Points>\n");
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  if(precision==2) {
    for(auto&& q : LV.q)
      print(vtpfile, "          %16.10e %16.10e %16.10e\n", q[0], q[1], q[2]);
  } else if(precision==1) {
    for(auto&& q : LV.q)
      print(vtpfile, "          %10.6e %10.6e %10.6e\n", q[0], q[1], q[2]);
  } else if(precision==0) {
    for(auto&& q : LV.q)
      print(vtpfile, "          %7.4e %7.4e %7.4e\n", q[0], q[1], q[2]);
  }

  print(vtpfile, "        </DataArray>\n");
  print(vtpfile, "      </Points>\n");
  
  print(vtpfile, "      <PointData Scalars=\"ScalarResults\" Vectors=\"VectorResults\">\n");

  // Write site id 
  print(vtpfile, "        <DataArray type=\"Int8\" NumberOfComponents=\"1\" Name=\"SiteID\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(LV.siteid.size() == LV.q.size());
  for(auto&& sid : LV.siteid)
    print(vtpfile, "%d ", sid);
  print(vtpfile, "\n");
  print(vtpfile, "        </DataArray>\n");

  // Write material id
  print(vtpfile, "        <DataArray type=\"Int8\" NumberOfComponents=\"1\" Name=\"MaterialID\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(LV.matid.size() == LV.q.size());
  for(auto&& id : LV.matid)
    print(vtpfile, "%d ", id);
  print(vtpfile, "\n");
  print(vtpfile, "        </DataArray>\n");

  // Write tag 
  print(vtpfile, "        <DataArray type=\"Int8\" NumberOfComponents=\"1\" Name=\"Tag\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(LV.tag.size() == LV.q.size());
  for(auto&& id : LV.tag)
    print(vtpfile, "%d ", id);
  print(vtpfile, "\n");
  print(vtpfile, "        </DataArray>\n");

  // Write atomic frequency (sigma)
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"1\" Name=\"AtomicFrequency\" format=\"ascii\">\n");
  print(vtpfile, "          ");
  assert(LV.sigma.size() == LV.q.size());
  if(precision==2) {
    for(auto&& sig : LV.sigma)
      print(vtpfile, "%16.10e ", sig);
  } else if(precision==1) {
    for(auto&& sig : LV.sigma)
      print(vtpfile, "%10.6e ", sig);
  } else if(precision==0) {
    for(auto&& sig : LV.sigma)
      print(vtpfile, "%7.4e ", sig);
  }
  print(vtpfile, "\n");
  print(vtpfile, "        </DataArray>\n");

  
  int ns = LV.nSpecies_max; //number of species (max)
  assert(ns>0);

  // Write molar fraction
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"%d\" "
                 "Name=\"MolarFraction\" format=\"ascii\">\n", ns);
  print(vtpfile, "          ");
  assert(LV.x.size() == LV.q.size());
  for(auto&& x : LV.x) {
    int i=0; 
    if(precision==2) {
      while(i<(int)x.size()) 
        print(vtpfile, "%16.10e ", x[i++]);
    } else if(precision==1) {
      while(i<(int)x.size()) 
        print(vtpfile, "%10.6e ", x[i++]);
    } else if(precision==0) {
      while(i<(int)x.size()) 
        print(vtpfile, "%7.4e ", x[i++]);
    }
    while(i++<ns)
      print(vtpfile, "0 ");
  } 
  print(vtpfile, "\n");
  print(vtpfile, "        </DataArray>\n");

  // Write chemical potential 
  print(vtpfile, "        <DataArray type=\"Float32\" NumberOfComponents=\"%d\" "
                 "Name=\"ChemicalPotential\" format=\"ascii\">\n", ns);
  print(vtpfile, "          ");
  assert(LV.gamma.size() == LV.q.size());
  for(auto&& gam : LV.gamma) {
    int i=0; 
    if(precision==2) {
      while(i<(int)gam.size()) 
        print(vtpfile, "%16.10e ", gam[i++]);
    } else if(precision==1) {
      while(i<(int)gam.size()) 
        print(vtpfile, "%10.6e ", gam[i++]);
    } else if(precision==0) {
      while(i<(int)gam.size()) 
        print(vtpfile, "%7.4e ", gam[i++]);
    }
    while(i++<ns)
      print(vtpfile, "0 ");
  } 
  print(vtpfile, "\n");
  print(vtpfile, "        </DataArray>\n");

  print(vtpfile, "      </PointData>\n");
  print(vtpfile, "    </Piece>\n");
  print(vtpfile, "  </PolyData>\n");
  print(vtpfile, "</VTKFile>\n");
  fclose(vtpfile);
  

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

  print("- Wrote solution on Lattice[%d] at time %e to %s.\n", lattice_id, time, full_fname);

  MPI_Barrier(comm);
}

//--------------------------------------------------------------------------------

void
Output::OutputSolutionXYZ(double time, [[maybe_unused]] int time_step, vector<LatticeVariables> &LVS) 
{
  for(int lat=0; lat<(int)LVS.size(); lat++)
    OutputSolutionOnLatticeXYZ(time, lat, LVS[lat]);
}

//--------------------------------------------------------------------------------

void
Output::OutputSolutionOnLatticeXYZ(double time, int lattice_id, LatticeVariables &LV) 
{
  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);

  // sanity checks
  auto np = LV.q.size();
  assert(LV.siteid.size() == np);
  assert(LV.matid.size() == np);
  assert(LV.sigma.size() == np);
  assert(LV.x.size() == np);
  assert(LV.gamma.size() == np);
  assert(LV.tag.size() == np);
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
  print(xyzfile, "tag:I:1:"); //tag

  int ns = LV.nSpecies_max; //number of species (max)
  assert(ns>0);

  print(xyzfile, "molar_frac:R:%d:", ns);
  print(xyzfile, "chem_potent:R:%d\" ", ns);

  print(xyzfile, "Time=%e\n", time);


  //Write results
  if(precision==2) {
    for(int i=0; i<(int)np; i++) {
      print(xyzfile, "%16s  ", material_name[LV.matid[i]].c_str());
      print(xyzfile, "%16.10e  %16.10e  %16.10e  ", LV.q[i][0], LV.q[i][1], LV.q[i][2]);
      print(xyzfile, "%16.10e  %16.10e  %16.10e  ", 
            LV.q[i][0]-LV.q0[i][0], LV.q[i][1]-LV.q0[i][1], LV.q[i][2]-LV.q0[i][2]);
      print(xyzfile, "%4d  %16.10e  %4d  ", LV.siteid[i], LV.sigma[i], LV.tag[i]);

      int s = 0; 
      while(s<(int)LV.x[i].size()) {
        print(xyzfile, "%16.10e  ", LV.x[i][s]);
        s++;
      }
      while(s<ns) {
        print(xyzfile, "%16.10e  ", 0.0);
        s++;
      }

      s = 0;
      while(s<(int)LV.gamma[i].size()) {
        print(xyzfile, "%16.10e  ", LV.gamma[i][s]);
        s++;
      }
      while(s<ns) {
        print(xyzfile, "%16.10e  ", 0.0);
        s++;
      }

      print(xyzfile,"\n");
    }
  }
  else if(precision==1) {
    for(int i=0; i<(int)np; i++) {
      print(xyzfile, "%16s ", material_name[LV.matid[i]].c_str());
      print(xyzfile, "%10.6e %10.6e %10.6e ", LV.q[i][0], LV.q[i][1], LV.q[i][2]);
      print(xyzfile, "%10.6e %10.6e %10.6e ", 
            LV.q[i][0]-LV.q0[i][0], LV.q[i][1]-LV.q0[i][1], LV.q[i][2]-LV.q0[i][2]);
      print(xyzfile, "%4d %10.6e %4d ", LV.siteid[i], LV.sigma[i], LV.tag[i]);

      int s = 0; 
      while(s<(int)LV.x[i].size()) {
        print(xyzfile, "%10.6e ", LV.x[i][s]);
        s++;
      }
      while(s<ns) {
        print(xyzfile, "%10.6e ", 0.0);
        s++;
      }

      s = 0;
      while(s<(int)LV.gamma[i].size()) {
        print(xyzfile, "%10.6e ", LV.gamma[i][s]);
        s++;
      }
      while(s<ns) {
        print(xyzfile, "%10.6e ", 0.0);
        s++;
      }

      print(xyzfile,"\n");
    }
  }
  else if(precision==0) {
    for(int i=0; i<(int)np; i++) {
      print(xyzfile, "%16s ", material_name[LV.matid[i]].c_str());
      print(xyzfile, "%7.4e %7.4e %7.4e ", LV.q[i][0], LV.q[i][1], LV.q[i][2]);
      print(xyzfile, "%7.4e %7.4e %7.4e ", 
      LV.q[i][0]-LV.q0[i][0], LV.q[i][1]-LV.q0[i][1], LV.q[i][2]-LV.q0[i][2]);
      print(xyzfile, "%4d %7.4e %4d ", LV.siteid[i], LV.sigma[i], LV.tag[i]);

      int s = 0; 
      while(s<(int)LV.x[i].size()) {
        print(xyzfile, "%7.4e ", LV.x[i][s]);
        s++;
      }
      while(s<ns) {
        print(xyzfile, "%7.4e ", 0.0);
        s++;
      }

      s = 0;
      while(s<(int)LV.gamma[i].size()) {
        print(xyzfile, "%7.4e ", LV.gamma[i][s]);
       s++;
      }
      while(s<ns) {
        print(xyzfile, "%7.4e ", 0.0);
        s++;
      }

      print(xyzfile,"\n");
    }
  }

  fclose(xyzfile);

  print("- Wrote solution on Lattice[%d] at time %e to %s.\n", lattice_id, time, full_fname);

  MPI_Barrier(comm);
}

//--------------------------------------------------------------------------------


//--------------------------------------------------------------------------------

