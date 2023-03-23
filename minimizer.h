#ifndef _MINIMIZER_H_
#define _MINIMIZER_H_
#include <petsctao.h>
#include <stdio.h>
#include <vector>
struct Vec3D;
struct Input;
using namespace std;

struct AppCtx
{
  Input *input;
  double dt;
  vector<vector<int> > *PdPd;
  vector<vector<int> > *PdH;
  vector<vector<int> > *HPd;
  vector<vector<int> > *HH;
  vector<double> *x0;
  int numDOFPerSite;
  int globDim;
  int locDim;
  int first_site, last_site_plus_one;
  double *alldof;
  vector<vector<Vec3D> > *QPp;
  double QWp;
  vector<double> *coe_Pd;
  vector<double> *coe_H;

    // by XS
    double *send_Pd, *recv_Pd;
    double *send_H, *recv_H;

  // the following vectors will be populated (repeatedly) when TAO calls "formFunctionAndGradient"
  vector<Vec3D> q_Pd0;
  vector<Vec3D> q_H0;
  vector<double> sigma_Pd0;
  vector<double> sigma_H0;
  vector<Vec3D> Fq_Pd0;
  vector<Vec3D> Fq_H0;
  vector<double> Fsigma_Pd0;
  vector<double> Fsigma_H0;

  AppCtx(Input *input_, vector<vector<int> > *PdPd_, vector<vector<int> > *PdH_,
         vector<vector<int> > *HPd_, vector<vector<int> > *HH_, vector<Vec3D> *q_Pd0_, 
         vector<Vec3D> *q_H0_, vector<double> *sigma_Pd0_, vector<double> *sigma_H0_, 
         vector<double> *x0_, vector<double> *coe_Pd_, vector<double> *coe_H_, 
         int numDOFPerSite_, int globDim_, int locDim_, int first_site_, 
         int last_site_plus_one_, double *alldof_, vector<vector<Vec3D> > *QPp_, double QWp_); 
  ~AppCtx();
}; 

class Minimizer
{
  MPI_Comm      comm;
  PetscMPIInt   size, rank;
  Tao           tao;
  AppCtx        *application_context; 
  Vec           yy;    /* all-inclusive solution vector */
  Vec           Ub, Lb; // variable bounds
 
  int numDOFPerSite;
  int globDim; // total # of dofs (not sites) to be solved
  int locDim; // # of dofs (not sites) to be handled by this CPU core
  int first_site, last_site_plus_one; //global index of first element, and last_plus_one
  double *alldof;
  
  void updateApplicationContext(double dt_, vector<Vec3D> *q_Pd0_, vector<Vec3D> *q_H0_, 
                                vector<double> *sigma_Pd0_, vector<double> *sigma_H0_,
                                vector<double> *x0_);
  static PetscErrorCode applicationData2Petsc(void* ptr/*AppCtx*/,
                                              vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, 
                                              vector<double> &sigma_Pd0,
                                              vector<double> &sigma_H0, Vec &yy);
  static PetscErrorCode petscData2Application(void* ptr/*AppCtx*/,
                                              Vec &yy, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H, 
                                              vector<double> &sigma_Pd, vector<double> &sigma_H);
  static PetscErrorCode formFunctionAndGradient(Tao, Vec/*input vector*/, PetscReal* /*func value*/,
                                         Vec/*gradient*/, void* /*AppCtx*/);
  static PetscErrorCode formVariableBounds(Tao, Vec Ub, Vec Lb, void* /*AppCtx*/);

public:
  Minimizer(int *argc, char ***argv);
  ~Minimizer();
  PetscErrorCode initialize(Input *input, vector<vector<int> > *PdPd, vector<vector<int> > *PdH,
                  vector<vector<int> > *HPd, vector<vector<int> > *HH, 
                  vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, vector<double> &sigma_Pd0,
                  vector<double> &sigma_H0, vector<double> &x0, vector<vector<Vec3D> > *QPp, double &QWp,
                  vector<double> *coe_Pd, vector<double> *coe_H); 
  PetscErrorCode minimizeFreeEntropy(double dt, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, 
                          vector<double> &sigma_Pd0, vector<double> &sigma_H0, vector<double> &x0, 
                          vector<Vec3D> &q_Pd, vector<Vec3D> &q_H, vector<double> &sigma_Pd, 
                          vector<double> &sigma_H, double &free_entropy);
  PetscErrorCode calculateChemicalPotential(vector<double> &x0, double gammabd, vector<double> &gamma, vector<int> &full_H);
  PetscErrorCode findLocalNeighbors(vector<vector<int> > &PdPd_ext, vector<vector<int> > &PdH_ext,
                  vector<vector<int> > &HPd_ext, vector<vector<int> > &HH_ext,
                  vector<vector<int> > &PdPd, vector<vector<int> > &PdH,
                  vector<vector<int> > &HPd, vector<vector<int> > &HH,
                  vector<vector<int> > &PdPdnum, vector<vector<int> > &HH1, int &maxNeib);
};
#endif
