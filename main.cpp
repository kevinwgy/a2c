/**********************************************************************************
 * Copyright © Kevin G. Wang, Xingsheng Sun, 2015
 * (1) Redistribution and use in source and binary forms, with or without modification,
 *     are permitted, provided that this copyright notice is retained.
 * (2) Use at your own risk.
 **********************************************************************************/
#include <stdio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "Vector3D.h"
#include "input.h"
#include "output.h"
#include "minimizer.h"
#include "hgversion.h"
using namespace std;

//--------------------------------------------------------------
void generateShells(int NC, vector<vector<Int3> > *shells);
void generateInterstitialShells(int NC, vector<vector<Int3> > *shells);
void initializeStateVariables(Input &input, vector<vector<Int3> > &SS1, vector<vector<Int3> > &SS2, 
                              vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, vector<double> &sigma_Pd0, 
                              vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, 
                              vector<int> &PdSubsurf, vector<int> &HSubsurf, 
                              double &gammabd, vector<int> &full_H); 
int findNeighbors(double a_Pd, double rc, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, 
                  vector<vector<int> > &PdPd, vector<vector<int> > &PdH, 
                  vector<vector<int> > &HPd, vector<vector<int> > &HH);
int enforceBoundaryConditions(Input &input, double dt, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0,
                              vector<double> &sigma_Pd0, vector<double> &sigma_H0, vector<double> &x0, 
                              vector<double> &gamma, vector<int> &HSubsurf, double &gammabd, vector<int> &full_H);
int integrateDiffusion(Input &input, double dt, vector<vector<int> > &HH1, vector<double> &gamma, vector<double> &x0, 
                       vector<double> &x, vector<int> &full_H);
int calculateMacroResults(vector<vector<int> > &PdPd, vector<vector<int> > &PdPdnum,
                          vector<Vec3D> &q_Pd, vector<double> &x, vector<int> &PdSubsurf, vector<int> &HSubsurf, 
                          double &xH, double &xHinterior);
void readResults(Input &input, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, vector<double> &sigma_Pd0,
                   vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, vector<int> &full_H);
void gaussianQuadrature(int dim, vector<vector<Vec3D> > &QP);
void printLogo(MPI_Comm comm);
//--------------------------------------------------------------


//--------------------------------------------------------------
// Main Function
//--------------------------------------------------------------
int main(int argc, char* argv[])
{
  clock_t start_time = clock(); //for timing purpose only
  //--------------------------------------------------------------
  // Setup MPI 
  //--------------------------------------------------------------
  MPI_Init(&argc, &argv);
  MPI_Comm comm;
  comm = MPI_COMM_WORLD;
  int MPI_rank, MPI_size;
  MPI_Comm_rank(comm, &MPI_rank);
  MPI_Comm_size(comm, &MPI_size);
  MPI_Barrier(comm);
  printLogo(comm);

  //--------------------------------------------------------------
  // Inputs
  //--------------------------------------------------------------
  Input input;
  input.readCmdLine(argc, argv);
  input.readCmdFile();

//  input.initialize(300.0/*T0*/, 1.0e-5/*eps*/, 5.35/*rc*/, 6/*N*/, 0.3/*xe*/, 4.3/*a_Pd*/,
//                   0.1/*sigma_Pd*/, 0.1/*sigma_H*/, 1.0e-6/*dt*/, 1.0e-5/*t_final*/,
//                   1/*output_frequency*/, "results"/*folder*/, "sol"/*filename_base*/);

  Output output(&comm, &input);

  //--------------------------------------------------------------
  // Create Shells
  //--------------------------------------------------------------
  int NC = 1.5*input.file.N*input.file.N + 4; //number of shells to be computed. TODO: enough?
  vector<vector<Int3> > SS1; // PdPd
  vector<vector<Int3> > SS2; // PdH
  generateShells(NC,&SS1);
  generateInterstitialShells(NC,&SS2);
 

  if(!MPI_rank) {
    if(input.file.sample_shape == 1) {
        cout << "- Shape of nanoparticle: cube." << endl; cout.flush();}
    else if(input.file.sample_shape == 2) {
        cout << "- Shape of nanoparticle: octahedron." << endl; cout.flush();}
    else if(input.file.sample_shape == 3) {
        cout << "- Shape of nanoparticle: sphere." << endl; cout.flush();}
    else {
        cerr << "Error: SampleShapeType in input.st must be 1 (cube), 2 (octahedron) or 3 (sphere)." << endl; cout.flush();
        exit(-1);}
    cout << "- Generated neighbor shells. " << endl; cout.flush();}

  //--------------------------------------------------------------
  // Initialize state: atomic positions, vibration frequencies, H molar fractions
  //--------------------------------------------------------------
  vector<Vec3D> q_Pd0, q_H0;
  vector<double> sigma_Pd0, sigma_H0, x0, gamma;
  vector<int> PdSubsurf, HSubsurf, full_H;
  double gammabd;

  initializeStateVariables(input, SS1, SS2, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, gamma, PdSubsurf, HSubsurf, gammabd, full_H); 


  if(!MPI_rank) {
    cout << "- Initialized atomic sites and state variables. " << endl;
    cout << "- Number of Pd sites: " << q_Pd0.size() << ". " << endl;
    cout << "- Number of H sites: " << q_H0.size() << ". " << endl;
    cout << "- Boundary conditions:" << endl;
    cout << "    Nondimensional chemical potential: " << gammabd << "." << endl;
    cout << "    Chemical potential: " << input.file.mubd << " eV." << endl;
    cout << "    Subsurface thickness: " << input.file.t_subsurf << " angstrom." << endl;
    cout << "    Number of Pd sites in the subsurface: " << PdSubsurf.size() << ". " << endl;
    cout << "    Number of H sites in the subsurface: " << HSubsurf.size() << ". " << endl;
    cout.flush();
  }

  //--------------------------------------------------------------
  // Initialize coefficients in correction term
  //--------------------------------------------------------------
  vector<double> coe_Pd, coe_H;
/*
  coe_Pd.push_back(0.0); //A_Pd
  coe_Pd.push_back(1.0334e-07); //B_PdPdPd
  coe_Pd.push_back(2.1336e-09); //B_PdPdH
  coe_Pd.push_back(4.5629e-07); //B_PdHH

  coe_H.push_back(0.0); //A_H
  coe_H.push_back(1.6805e-06); //B_HPdPd
  coe_H.push_back(1.8245e-08); //B_HPdH
  coe_H.push_back(9.8190e-08); //B_HHH
*/

  coe_Pd.push_back(0.04); //A_Pd
  coe_Pd.push_back(1.0334e-07); //B_PdPdPd
  coe_Pd.push_back(2.1336e-09); //B_PdPdH
  coe_Pd.push_back(4.5629e-07); //B_PdHH

  coe_H.push_back(0.014); //A_H
  coe_H.push_back(1.6805e-06); //B_HPdPd
  coe_H.push_back(1.8245e-08); //B_HPdH
  coe_H.push_back(9.8190e-08); //B_HHH

  //--------------------------------------------------------------
  // Find neighbors of each atom within cut-off distance rc
  // Note: This algorithm does NOT rely on shells or the (i,j,k) indices. Hence allows amorphous layers
  //--------------------------------------------------------------
  double ext_dist = 1.5;
  double rc_ext = input.file.rc + ext_dist;
  vector<vector<int> > PdPd_ext;
  vector<vector<int> > PdH_ext; 
  vector<vector<int> > HPd_ext; 
  vector<vector<int> > HH_ext; 
 
  int maxNeighbors_ext = findNeighbors(input.file.a_Pd, rc_ext, q_Pd0, q_H0, PdPd_ext, PdH_ext, HPd_ext, HH_ext);

  if(!MPI_rank) {
    cout << "- Done with extended neighbor search using a KDTree. " << endl;
    cout << "  Extented radius: " << rc_ext << ". Max. number of neighbors: " << maxNeighbors_ext << ". " << endl;
    cout.flush();}

  //--------------------------------------------------------------
  // Calculate Gaussian Quadrature points and weights with 6 variables
  //--------------------------------------------------------------
  vector<vector<Vec3D> > QPp; // for function with 6 variables.
  gaussianQuadrature(3*2, QPp);
  double QWp = 0.5/(3.0*2.0);

  //--------------------------------------------------------------
  // Initialize optimization solver
  //--------------------------------------------------------------
  vector<vector<int> > PdPd;//Pd neighbors of Pd;
  vector<vector<int> > PdH; //H neighbors of Pd;
  vector<vector<int> > HPd; //Pd neighbors of H;
  vector<vector<int> > HH;  // H neighbors of H;
  vector<vector<int> > PdPdnum;  // number of Pd neighbors on different shells, for local lattice constant;
  vector<vector<int> > HH1;  // first H neighbor shell of H, for diffusion;
  int err = 0; //error code

  Minimizer minimizer(&argc, &argv);

  int maxNeighbors = 0;

/*
if(!MPI_rank) {
//    for (int i=0; i<PdPd.size(); i++)
//        cout << i << " " << PdPd[i].size()+PdH[i].size() << endl;
  //  for (int i=0; i<HPd.size(); i++)
    //    cout << i << " " << HPd[i].size()+HH[i].size() << endl;
    cout.flush();}
*/

  //--------------------------------------------------------------
  // open file to write average H/Pd ratio and lattice constant
  //--------------------------------------------------------------
  ofstream result_file("Result.txt", ios::app);
  if(!result_file){
    cerr << "WARNING! Cannot open the Result file!" <<  endl;
  }

  //--------------------------------------------------------------
  // Check start or restart for initialization
  //--------------------------------------------------------------
  int iFrame = 0;
  int iTimeStep = 0;
  int startTimeStep = 0;
  double t = 0.0; //current simulation time
  double free_entropy = 0.0;
  double &dt = input.file.dt; // by XS
 
  vector<Vec3D> q_Pd(q_Pd0), q_H(q_H0); //to be used for storing new states
  vector<double> sigma_Pd(sigma_Pd0), sigma_H(sigma_H0), x(x0); //same as above
  

  double xH = 0.0;
  double xHinterior = 0.0;

  double xH0 = 10.0;

  if (input.file.restart == 0) {
    if(!MPI_rank) {
      cout << "---------------------------------------" << endl;
      cout << "- Attention: start a new calculation! -" << endl;
      cout << "---------------------------------------" << endl;
      cout.flush();
    }

    minimizer.initialize(&input, &PdPd, &PdH, &HPd, &HH, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, &QPp, QWp, &coe_Pd, &coe_H);

    err = minimizer.findLocalNeighbors(PdPd_ext, PdH_ext, HPd_ext, HH_ext, PdPd, PdH, HPd, HH, PdPdnum, HH1, maxNeighbors);
    if(!MPI_rank) {
      cout << "- Done with local neighbor search using a KDTree. Max. number of neighbors: " << maxNeighbors << ". " << endl;
      cout.flush();}
    MPI_Barrier(comm);

    // Enforce boundary conditions
    err = enforceBoundaryConditions(input, dt, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, gamma, HSubsurf, gammabd, full_H);

    if(!MPI_rank) {
      cout << "- Done with Boundary Condition." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);

    // Optimize the problem with initial condition
    if(!MPI_rank) {
      cout << "- Optimization for initialization: " << endl; cout.flush();}

    err = minimizer.minimizeFreeEntropy(0.0, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, q_Pd, q_H, sigma_Pd, sigma_H, free_entropy);
    if(!MPI_rank){
      cout << "- Done with Max-Ent." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);

    q_Pd0 = q_Pd;
    q_H0 = q_H;
    sigma_Pd0 = sigma_Pd;
    sigma_H0 = sigma_H;

    err = minimizer.calculateChemicalPotential(x0, gammabd, gamma, full_H);
    if(!MPI_rank) {
      cout << "- Done with Chemical Potential." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);

    err = calculateMacroResults(PdPd, PdPdnum, q_Pd0, x0, PdSubsurf, HSubsurf, xH, xHinterior);
    // write time history of macroscropic results to file
    if(!MPI_rank) {
      result_file << "  " << scientific << t << "  " << scientific << xH << "  " << scientific << xHinterior << " " << scientific << free_entropy << endl;
      result_file.flush();
    }
    // write solution to file
    output.output_solution(iFrame++, iTimeStep, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, gamma, full_H);
    MPI_Barrier(comm);
  }
  else if (input.file.restart == 1) {
    if(!MPI_rank) {
      cout << "-------------------------------------------" << endl;
      cout << "- Attention: restart the old calculation! -" << endl;
      cout << "-------------------------------------------" << endl;
      cout.flush();
    }
    iFrame = input.file.restart_file_num+1;
    iTimeStep = input.file.restart_file_num*input.file.output_frequency; 
    startTimeStep = iTimeStep;
    t = input.file.restart_file_num*dt*input.file.output_frequency; 
    readResults(input, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, gamma, full_H);
    MPI_Barrier(comm);
    if(!MPI_rank) {
      cout << "- Read restart input file: " << input.file.restart_file_num << "." << endl; cout.flush();}

    minimizer.initialize(&input, &PdPd, &PdH, &HPd, &HH, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, &QPp, QWp, &coe_Pd, &coe_H);
 
    err = minimizer.findLocalNeighbors(PdPd_ext, PdH_ext, HPd_ext, HH_ext, PdPd, PdH, HPd, HH, PdPdnum, HH1, maxNeighbors);
    MPI_Barrier(comm);
    if(!MPI_rank) {
        cout << "- Done with local neighbor search using a KDTree. Max. number of neighbors: " << maxNeighbors << ". " << endl;
        cout.flush();}
  }
  else {
    if(!MPI_rank) {
      cerr << "ERREOR! Wrong Restart in the input file!" <<  endl;
      exit(-1);
    }
  }
  
  x = x0;

  //--------------------------------------------------------------
  // MAIN LOOP: Integration in time
  //--------------------------------------------------------------
  if(!MPI_rank) {
    cout << endl;
    cout << "----------------------------" << endl;
    cout << "--       Main Loop        --" << endl;
    cout << "----------------------------" << endl;
    cout.flush();
  }
  
  bool lastStep = false;//false;

  while(!lastStep) {
    t += dt;
    iTimeStep++;
    if(!MPI_rank) {
      cout << endl;
      cout << "* Time-Step " << iTimeStep << ". Time: " << t << endl; cout.flush();}

    // check if this is the last time step
    //  double dt = input.file.dt;
    if(t >= input.file.t_final - 0.001*dt) {//last time step
      lastStep = true;
      //dt = input.file.t_final - t;
    }


    //--------------------------------------------------------------
    // Enforce boundary conditions
    // Update the boundary values in x
    //--------------------------------------------------------------
    err = enforceBoundaryConditions(input, dt, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, gamma, HSubsurf, gammabd, full_H);

    if(!MPI_rank) {
      cout << "# Done with Boundary Condition." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);


    //--------------------------------------------------------------
    // Integrate discrete diffusion law for one time step
    // Update x 
    //--------------------------------------------------------------
    err = integrateDiffusion(input, dt, HH1, gamma, x0, x, full_H);

    if(!MPI_rank) {
        cout << "# Done with Diffusion."  << endl;  cout.flush();}
    MPI_Barrier(comm);
    // update old states and advance in time
    x0 = x;

    //--------------------------------------------------------------
    // Minimize grand-canonical free energy (equiv. to free entropy at constant T)  
    // Update q_Pd, q_H, sigma_Pd, sigma_H
    //--------------------------------------------------------------
    free_entropy = 0.0;
    err = minimizer.minimizeFreeEntropy(dt, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, q_Pd, q_H, sigma_Pd, sigma_H, free_entropy);
                                        //x should not get changed...
    if(!MPI_rank){
        cout << "# Done with Max-Ent." << endl;
        cout.flush();
    }
    MPI_Barrier(comm);
    // update old states and advance in time
    q_Pd0     = q_Pd;
    q_H0      = q_H;
    sigma_Pd0 = sigma_Pd;
    sigma_H0  = sigma_H;

    err = minimizer.findLocalNeighbors(PdPd_ext, PdH_ext, HPd_ext, HH_ext, PdPd, PdH, HPd, HH, PdPdnum, HH1, maxNeighbors);
          
    if(!MPI_rank) {
        cout << "# Done with local neighbor search using a KDTree. Max. number of neighbors: " << maxNeighbors << ". " << endl;
        cout.flush();}
    MPI_Barrier(comm);

    //--------------------------------------------------------------
    // Calculate chemical potential
    // Update gamma
    //--------------------------------------------------------------
    err = minimizer.calculateChemicalPotential(x0, gammabd, gamma, full_H);

    if(!MPI_rank) {
      cout << "# Done with Chemical Potential." << endl;
      cout << endl;
      cout.flush();
    }
    MPI_Barrier(comm);
/*
    if(!MPI_rank) {
     // for(int i=0; i<gamma.size();i++)
     // cout << "site" << i << ": " << gamma[i] << " " << x0[i] << " " << HH[i].size() << endl;
      cout << gammabd << " " << gamma[0] << " " << x0[0] << " " << HH[0].size() << " " << HH1[0].size() <<  endl;
      cout << gammabd << " " << gamma[gamma.size()-1] << " " << x0[gamma.size()-1] << " " << HH[gamma.size()-1].size() << " " << HH1[gamma.size()-1].size() << endl;
      cout.flush();
    }
*/    
    //--------------------------------------------------------------
    // write solution to file
    // Calculate H mole ratio, lattice constant
    //--------------------------------------------------------------
    if((iTimeStep%input.file.output_frequency == 0) || lastStep) {
      // Calculate H mole ratio, lattice constant
      err = calculateMacroResults(PdPd, PdPdnum, q_Pd0, x0, PdSubsurf, HSubsurf, xH, xHinterior);
      // check reaching equilibrium or not
      if(fabs(xH0-xH)<1.0e-6) lastStep = true; 
      else xH0 = xH;  
      // write time history of macroscropic results to file
      if(!MPI_rank) {
        result_file << "  " << scientific << t << "  " << scientific << xH << "  " << scientific << xHinterior << " " << scientific << free_entropy << endl;
        result_file.flush();
      }
      // write solution to file
      output.output_solution(iFrame++, iTimeStep, q_Pd0, q_H0, sigma_Pd0, sigma_H0, x0, gamma, full_H);
    }
    MPI_Barrier(comm);
  }
  
  if(!MPI_rank) {
    cout << "----------------------------" << endl;
    cout << "--   End of Main Loop     --" << endl;
    cout << "----------------------------" << endl;
    cout.flush();
  }

  //--------------------------------------------------------------
  // close result file
  //-------------------------------------------------------------- 
  result_file.close();

  if(!MPI_rank) {
    cout << "======================================" << endl;
    cout << "   NORMAL TERMINATION (t = " << t << ") " << endl;
    cout << "======================================" << endl;
  }
  MPI_Barrier(comm);
  if(!MPI_rank) cout << "Total Computation Time: " << ((double)(clock()-start_time))/CLOCKS_PER_SEC << " sec." << endl;
//  MPI_Finalize(); //will be called in the destructor of Tao

  return 0;
}

void printLogo(MPI_Comm comm)
{
  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);
  if(!MPI_rank) {
    cout << endl;
    cout << "  ____                 _ _      _      _    ____   ____ " << endl;
    cout << " |  _ \\ __ _ _ __ __ _| | | ___| |    / \\  |___ \\ / ___|" << endl;
    cout << " | |_) / _` | '__/ _` | | |/ _ \\ |   / _ \\   __) | |    " << endl;
    cout << " |  __/ (_| | | | (_| | | |  __/ |  / ___ \\ / __/| |___ " << endl;
    cout << " |_|   \\__,_|_|  \\__,_|_|_|\\___|_| /_/   \\_\\_____|\\____|" << endl;
    cout << endl;
    cout << " Revision: " << hgRevisionNo << " | " << hgRevisionHash <<endl;
    cout << endl;
    cout.flush();
  }
}
