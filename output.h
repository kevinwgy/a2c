#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
struct Vec3D;
struct Input;

class Output {
  MPI_Comm *comm;
  char filename_base[128];
  char full_filename_base[128];
  ofstream summaryfile;

public:
  Output(MPI_Comm *comm_, Input *input, int argc, char* argv[]);
  ~Output();
  
  void output_solution(int iFrame, int iTimeStep, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H,
                     vector<double> &sigma_Pd, vector<double> &sigma_H, vector<double> &x,
                     vector<double> &gamma, vector<int> &full_H); 

  const string getCurrentDateTime();
  void printLogo(int argc, char *argv[]);
};
#endif
