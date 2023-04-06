#include <Lattice.h>
#include <Utils.h>
#include <iostream>
#include <cmath>
using std::vector;

//-------------------------------------------------------------------------

Lattice::Lattice(IoData &iod_)
      : iod(iod_), lattice_type(iod_.space.lattice_type)
{

  double PI = 2.0*acos(0.0);

  if(iod.space.lattice_constant_1<=0.0) {
    print_error("*** Error: Need at least one positive lattice constant.\n");
    exit_mpi();
  }

  // Get Lattice constants and vectors
  switch (lattice_type) {
    case SpaceData::SC

 I AM HERE!!!

  };


  // Create Shells
  switch (iod.space.lattice_type) {
    case SpaceData::HCP :
      GenerateLatticeShellsHCP(number_of_shells, &SS1);
      if(SS2)
        GenerateOctahedralIntersticesHCP(number_of_shells, SS2);
      if(SS3)
        GenerateTetrahedralIntersticesHCP(number_of_shells, SS3);
      break;
    case SpaceData::FCC :
      GenerateLatticeShellsFCC(number_of_shells, &SS1);
      if(SS2)
        GenerateOctahedralIntersticesFCC(number_of_shells, SS2);
      if(SS3)
        GenerateTetrahedralIntersticesFCC(number_of_shells, SS3);
      break;
    default :
      print_error("*** Error: Cannot create shells for the specified lattice type.\n");      
      exit_mpi();
  }
}
}

//-------------------------------------------------------------------------

Lattice::~Lattice()
{ }

//-------------------------------------------------------------------------

void
Lattice::GenerateShells(vector<vector<Int3> > &SS1, vector<vector<Int3> >* SS2,
                       vector<vector<Int3> >* SS3)
{
  int number_of_shells = std::min(1.5*iod.space.N*iod.space.N + 4, 100);

  switch (iod.space.lattice_type) {
    case SpaceData::HCP :
      GenerateLatticeShellsHCP(number_of_shells, &SS1);
      if(SS2)
        GenerateOctahedralIntersticesHCP(number_of_shells, SS2);
      if(SS3)
        GenerateTetrahedralIntersticesHCP(number_of_shells, SS3);
      break;
    case SpaceData::FCC :
      GenerateLatticeShellsFCC(number_of_shells, &SS1);
      if(SS2)
        GenerateOctahedralIntersticesFCC(number_of_shells, SS2);
      if(SS3)
        GenerateTetrahedralIntersticesFCC(number_of_shells, SS3);
      break;
    default :
      print_error("*** Error: Cannot create shells for the specified lattice type.\n");      
      exit_mpi();
  }
}

//-------------------------------------------------------------------------

void
Lattice::GenerateLatticeShellsHCP(int NC, vector<vector<Int3> > *shells)
{

  assert(shells && shells->empty());
  assert(NC>=0);

  int p = 1; //step size

  for(int iShell=0; iShell<NC; iShell++) {
    shells->push_back(); //start a new shell
    vector<Int3>& shell(shells->back());
    while(shell.empty()) { //find the iShell
      int lmax = (int)floor(sqrt((double)p));
      for (int l1=-lmax; l1<=lmax; l1++) 
        for (int l2=-lmax; l2<=lmax; l2++) 
          for (int l3=-lmax; l3<=lmax; l3++)
            if((l1*l1+l2*l2+l3*l3 == p) && ((l1+l2+l3)%2 == 0))
              shell.push_back(Int3(l1,l2,l3));
      if(shell.empty())
        p++;
    }
    p++;
  }
  if((int)shells->size()!=NC) {
    cerr << "ERROR: size = " << shells->size() << ", NC = " << NC << endl;
    exit(-1);
  }
}

//-------------------------------------------------------------------------

void
Lattice::GenerateOctahedralIntersticesHCP(int NC, vector<vector<Int3> > *shells)
{
  int p = 1; //step size

  if(!shells->empty()) {
    cerr << "ERROR: shells not empty!" << endl;
    shells->clear();
  }
  for(int iShell=0; iShell<NC; iShell++) {
    vector<Int3> shell;
    while(shell.empty()) { //find the iShell
      int lmax = (int)floor(sqrt((double)p));
      for (int l1=-lmax; l1<=lmax; l1++)
        for (int l2=-lmax; l2<=lmax; l2++)
          for (int l3=-lmax; l3<=lmax; l3++) 
            if((l1*l1+l2*l2+l3*l3 == p) && ((l1+l2+l3)%2 != 0)) //Note: % behaves differently with 'mod' in MATLAB for negative numbers
              shell.push_back(Int3(l1,l2,l3));
      if(shell.empty())
        p++;
    }
    shells->push_back(shell);
    p++;
  }
  if((int)shells->size()!=NC) {
    cerr << "ERROR: size = " << shells->size() << ", NC = " << NC << endl;
    exit(-1);
  }
//  for(int i=0; i<shells->size(); i++) 
//    for(int j=0; j<(*shells)[i].size(); j++)
//      cout << "shell " << i << ", id " << j << ": " << (*shells)[i][j][0] << " " << (*shells)[i][j][1] << " " << (*shells)[i][j][2] << endl;
}
