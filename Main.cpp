/**********************************************************************************
 * Copyright ©Xingsheng Sun, Kevin G. Wang, 2015
 **********************************************************************************/
#include <Output.h>
#include <SpaceOperator.h>
#include <cstring>
#include <Utils.h>
#include <set>
using std::vector;

int verbose;
clock_t start_time;
MPI_Comm a2c_comm;

//--------------------------------------------------------------
// Main Function
//--------------------------------------------------------------
int main(int argc, char* argv[])
{
  start_time = clock(); //for timing purpose only
  //--------------------------------------------------------------
  // Setup MPI 
  //--------------------------------------------------------------
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  a2c_comm = MPI_COMM_WORLD; //Different if A2C is coupled w/ something else
  int mpi_rank(-1);
  MPI_Comm_rank(comm, &mpi_rank);

  //--------------------------------------------------------------
  // Prints A2C header
  //--------------------------------------------------------------
  printHeader(argc, argv);

  //--------------------------------------------------------------
  // Inputs
  //--------------------------------------------------------------
  IoData iod(argc, argv);
  iod.finalize();
  verbose = iod.output.verbose;

  //------------------------------------------------------------------------
  // Get global species ID & Name
  //--------------------------------------------------------------
  int nSpecies = iod.speciesMap.dataMap.size();
  vector<string> global_species(nSpecies);
  std::set<int> species_tracker; //for error detection only
  for(auto&& sp : iod.speciesMap.dataMap) {
    int spid = sp.first;
    if(spid<0 || spid>=nSpecies) {
      print_error("*** Error: Detected error in global species index (id = %d).\n", spid);
      exit_mpi();
    }
    species_tracker.insert(spid);
    global_species[spid] = string(sp.second->name);
  }
  if((int)species_tracker.size() != nSpecies) {
    print_error("*** Error: Detected error in global species indices.\n");
    exit_mpi();
  } 


  //--------------------------------------------------------------
  // Setup materials
  //--------------------------------------------------------------
  int nMaterials = iod.materialMap.dataMap.size();
  vector<MaterialOperator> mato(nMaterials);
  std::set<int> material_tracker; //for error detection only
  for(auto&& mat : iod.materialMap.dataMap) {
    int matid = mat.first;
    if(matid<0 || matid>=nMaterials) {
      print_error("*** Error: Detected error in material index (id = %d).\n", matid);
      exit_mpi();
    }
    material_tracker.insert(matid);

    mato[matid].Setup(matid, *mat.second, global_species);
  }
  if((int)material_tracker.size() != nMaterials) {
    print_error("*** Error: Detected error in material indices.\n");
    exit_mpi();
  } 


  //--------------------------------------------------------------
  // Setup lattice constants
  //--------------------------------------------------------------
  int nLattices = iod.latticeMap.dataMap.size();
  vector<LatticeStructure> lats(nLattices);
  std::set<int> lattice_tracker; //for error detection only
  for(auto&& lat : iod.latticeMap.dataMap) {
    int lattice_id = lat.first;
    if(lattice_id<0 || lattice_id>=nLattices) {
      print_error("*** Error: Detected error in lattice index (id = %d).\n", lattice_id);
      exit_mpi();
    }
    lattice_tracker.insert(lattice_id);
    lats[lattice_id].Setup(lattice_id, *lat.second, mato);
  }
  if((int)lattice_tracker.size() != nLattices) {
    print_error("*** Error: Detected error in lattice indices.\n");
    exit_mpi();
  }


  //--------------------------------------------------------------
  // Setup outputs
  //--------------------------------------------------------------
  Output output(comm, iod);

  //--------------------------------------------------------------
  // Create the ``specimen''
  //--------------------------------------------------------------
  vector<LatticeVariables> LVS(nLattices);
  SpaceOperator spo(comm, iod, mato, lats);
  spo.SetupLatticeVariables(LVS);



  print("\n");
  print("----------------------------\n");
  print("--       Main Loop        --\n");
  print("----------------------------\n");
  double t = 0.0, dt = 0.0;
  int time_step = 0;

  output.OutputSolution(t, dt, time_step, LVS, true); //force-write

  // trial run
  for(int i=1; i<argc; i++)
    if(!strcmp(argv[i], "-t")) 
      goto THE_END;






  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m   NORMAL TERMINATION (t = %e)  \033[0m\n", t);
  print("\033[0;32m==========================================\033[0m\n");
  print("Total Computation Time: %f sec.\n", ((double)(clock()-start_time))/CLOCKS_PER_SEC);
  print("\n");




THE_END:

  int mpi_finalized(0);
  MPI_Finalized(&mpi_finalized);
  if(!mpi_finalized)
    MPI_Finalize();

  return 0;
}


// -------------------------------------------------------------------------------
