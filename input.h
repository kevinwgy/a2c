#ifndef _INPUT_H_
#define _INPUT_H_
#include <stdio.h>

struct InputFileData
{
  // constants
  double T0;  // temperature (K)
  double eps; // tolerance
  double rc;  // cutoff distance (Ang)
  double dt;  // time-step (s)  TODO: if dt drops below 1e-8, we should change unit to reduce roundoff error
  double t_final;  // final time (unit: same as dt)
  double xUb; // upper bound of hydrogen atomic fraction
  double xLb; // lower bound of hydrogen atomic fraction

  // computational domain
  int sample_shape;
  int N; //N*N*N = number of unit cells;

  // initial condition
  double xe; //hydrogen molar fraction
  double a_Pd; //Pd lattice constant (Ang)
  double sigma_Pd; //Pd frequency (Ang^-2)
  double sigma_H; //H frequency (Ang^-2)

  // minimize frequency
  int minimize_frequency;
  int linear_extrap;

  // output
  int output_frequency;
  const char* foldername;
  const char* filename_base;

  // absorption
  double v; // attempt frequency (s^-1)
  double Qm; // activation energy (eV)

  // adsorption
  double vad; // attempt frequency (s^-1)
  double Qad; // activation energy (eV) 

  // boundary condition
  int subsurf_equil;
  double t_subsurf;
  double mubd;

  // find local neighbors
  int local_nei;

  // continue setup
  int restart;
  int restart_file_num;

  // number of sites in restart file
  int nPd_restart;
  int nH_restart;

  InputFileData();
  ~InputFileData();

  void setup(const char *);

  // the following function is obsolete. 
/*
  void initialize(double T0_, double eps_, double rc_, int N_, double xe_, double a_Pd_,
                  double sigma_Pd_, double sigma_H_, double dt_, double t_final_,
                  int output_frequency_, char* foldername_, char *filename_base_);
*/
};


class Input
{
  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  InputFileData  file;
  
public:

  Input() {}
  ~Input() {}

  void readCmdLine(int, char**);
  void readCmdFile();
  void setupCmdFileVariables();
  
};
#endif
