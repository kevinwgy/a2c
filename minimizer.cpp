#include <iostream>
#include <math.h>
#include <algorithm>
#include <Vector3D.h>
#include <minimizer.h>
#include <input.h>
#include <potential.h>
#include <KDTree.h>
using namespace std;
static char help[] = "minimization of grand-canonical free entropy\n";

double electronDensity(char atype, int n, double rc, vector<vector<int> > &PdPd, 
                    vector<vector<int> > &PdH, vector<vector<int> > &HPd, 
                    vector<vector<int> > &HH, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H, 
                    vector<double> &sigma_Pd, vector<double> &sigma_H, vector<double> &x, 
                    vector<vector<Vec3D> > &QPp, double &QWp);

AppCtx::AppCtx(Input *input_, vector<vector<int> > *PdPd_, vector<vector<int> > *PdH_,
               vector<vector<int> > *HPd_, vector<vector<int> > *HH_, vector<Vec3D> *q_Pd0_, 
               vector<Vec3D> *q_H0_, vector<double> *sigma_Pd0_, vector<double> *sigma_H0_, 
               vector<double> *x0_, vector<double> *coe_Pd_, vector<double> *coe_H_, int numDOFPerSite_, int globDim_, int locDim_, 
               int first_site_, int last_site_plus_one_, double *alldof_, 
               vector<vector<Vec3D> > *QPp_, double QWp_) : 
  input(input_), PdPd(PdPd_), PdH(PdH_), HPd(HPd_), HH(HH_), coe_Pd(coe_Pd_), coe_H(coe_H_), x0(x0_), QPp(QPp_), send_Pd(NULL),recv_Pd(NULL), send_H(NULL), recv_H(NULL)
{
  dt = input->file.dt; 
  numDOFPerSite = numDOFPerSite_;
  globDim = globDim_;
  locDim = locDim_;
  first_site = first_site_;
  last_site_plus_one = last_site_plus_one_;
  alldof = alldof_;
  QWp = QWp_;
  
  // allocate memory
  q_Pd0      = *q_Pd0_;
  q_H0       = *q_H0_;
  sigma_Pd0  = *sigma_Pd0_;
  sigma_H0   = *sigma_H0_;
  Fq_Pd0     = *q_Pd0_;
  Fq_H0      = *q_H0_;
  Fsigma_Pd0 = *sigma_Pd0_;
  Fsigma_H0  = *sigma_H0_;
  
    // by XS
    // allocate memory for  electron density of each site.
    send_Pd = new double[q_Pd0.size()];
    recv_Pd = new double[q_Pd0.size()];
    send_H = new double[q_H0.size()];
    recv_H = new double[q_H0.size()];
}

AppCtx::~AppCtx() {
    delete[] send_Pd;
    delete[] recv_Pd;
    delete[] send_H;
    delete[] recv_H;
}

Minimizer::Minimizer(int *argc, char ***argv)
{
  [[maybe_unused]] PetscErrorCode ierr;
  PetscInitialize(argc, argv, *argc>=3 ? (*argv)[2] : (char*)0, help); 
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  comm = PETSC_COMM_WORLD; //same as MPI_COMM_WORLD
  application_context = NULL;
  alldof = NULL;
}


Minimizer::~Minimizer()
{
  if(alldof) delete[] alldof;
  [[maybe_unused]] PetscErrorCode ierr;
  ierr = TaoDestroy(&tao);
  ierr = VecDestroy(&yy);
  ierr = VecDestroy(&Ub);
  ierr = VecDestroy(&Lb);
}


PetscErrorCode Minimizer::formFunctionAndGradient([[maybe_unused]] Tao tao, Vec X, PetscReal* fun, Vec Grad, void* ptr)
{
  /* TODO: calculate f and Grad at X */
  // THIS IS THE MOST IMPORTANT FUNCTION
  // THIS FUNCTION MUST BE PARALLELIZED, AND EFFICIENT!
  
  PetscErrorCode ierr = 0;
  AppCtx *app = (AppCtx *) ptr;
  EAM eam; // by XS
  const double kB = 8.6173324e-5;
  double &T = app->input->file.T0;
  int nPd = app->q_Pd0.size();
  int nH  = app->q_H0.size();
  vector<vector<int> > &PdPd = *(app->PdPd); //create reference
  vector<vector<int> > &PdH = *(app->PdH); 
  vector<vector<int> > &HPd = *(app->HPd); 
  vector<vector<int> > &HH = *(app->HH);
  vector<double> &x0 = *(app->x0);
  vector<vector<Vec3D> > &QPp = *(app->QPp);
  double &QWp = app->QWp;
 
  double fun_local = 0.0;
  double sum_V = 0.0; // potential energy
  double sum_sigma = 0.0;  

  // populate q_Pd0, q_H0, ..., in the application context
  ierr = petscData2Application(app, X, app->q_Pd0, app->q_H0, app->sigma_Pd0, app->sigma_H0);
  CHKERRQ(ierr);

  // clear force vectors
    for(int i=0; i<nPd; i++) {
        app->Fq_Pd0[i] = 0.0;
        app->Fsigma_Pd0[i] = 0.0;
    }
    for(int i=0; i<nH; i++) {
        app->Fq_H0[i] = 0.0;
        app->Fsigma_H0[i] = 0.0;
    } 

    // initiallization for electron density calculation
    for(int i=0; i<nPd; i++)
        app->send_Pd[i] = 0.0;
    for(int i=0; i<nH; i++)
        app->send_H[i] = 0.0;

/*
fun_local and app->F should be revised according to EAM potential.
fun_local is the total free entropy and app->Fq is its gradient with respect to positions.
app->Fsigma is the gradient with respect to sigma.
*/

  // ------------------------------------------------------------
  // Objective function = EAM potential, revised by XS.
  // ------------------------------------------------------------i
    // ------------------------------------------------------------
    // Parallelization I: calculate electron density of each site.
    // ------------------------------------------------------------
    if(app->first_site<nPd) { // the first local site is a Pd
        for(int i=app->first_site; i<min(nPd, app->last_site_plus_one); i++) {
            app->send_Pd[i] = electronDensity('P', i, app->input->file.rc, PdPd, PdH, HPd, HH, app->q_Pd0, app->q_H0, app->sigma_Pd0, app->sigma_H0, x0, QPp, QWp);
        }
        for(int ii=nPd; ii<app->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nPd; // index in the H vectors
            app->send_H[i] = electronDensity('H', i, app->input->file.rc, PdPd, PdH, HPd, HH, app->q_Pd0, app->q_H0, app->sigma_Pd0, app->sigma_H0, x0, QPp, QWp);
       }
    }
    else { //all the local sites are H
        for(int ii=app->first_site; ii<app->last_site_plus_one; ii++) {
            int i = ii - nPd; // index in the H vectors
            app->send_H[i] = electronDensity('H', i, app->input->file.rc, PdPd, PdH, HPd, HH, app->q_Pd0, app->q_H0, app->sigma_Pd0, app->sigma_H0, x0, QPp, QWp);
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(app->send_Pd, app->recv_Pd, nPd, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(app->send_H, app->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
    // ------------------------------------------------------------
    // Parallelization II: calculate objective value and gradients.
    // ------------------------------------------------------------
    if(app->first_site<nPd) { // the first local site is a Pd
        for(int i=app->first_site; i<min(nPd, app->last_site_plus_one); i++) {
        // calc. force for Site i
            Vec3D &qi = app->q_Pd0[i];
            double &sigmai = app->sigma_Pd0[i];

            // Gaussian quadrature for embedding energy
            sum_V += eam.F_Pd(app->recv_Pd[i]);
 
            // Pd-Pd pairs
            for(int j=0; j<(int)PdPd[i].size(); j++) {
                int &jj = PdPd[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_Pd0[jj];
                double &sigmaj = app->sigma_Pd0[jj];
 
                for(int k=0; k<(int)QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    app->Fq_Pd0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[i])*eam.f_Pd_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_Pd0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[i])*eam.f_Pd_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    app->Fq_Pd0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[jj])*eam.f_Pd_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_Pd0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[jj])*eam.f_Pd_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    
                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*QWp*eam.phi_Pd(r.norm());
                    
                    // calculate gradient for pair energy
                    app->Fq_Pd0[i] += QWp*eam.phi_Pd_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_Pd0[i] += QWp*eam.phi_Pd_deriv(r.norm())/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }

            // Pd-H pairs
            for(int j=0; j<(int)PdH[i].size(); j++) {
                int &jj = PdH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_H0[jj];
                double &sigmaj = app->sigma_H0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    app->Fq_Pd0[i] += x0[jj]*QWp*eam.F_Pd_deriv(app->recv_Pd[i])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_Pd0[i] += x0[jj]*QWp*eam.F_Pd_deriv(app->recv_Pd[i])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    app->Fq_Pd0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj])*eam.f_Pd_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_Pd0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj])*eam.f_Pd_deriv(r.norm())/r.norm()*r*QPp[k][0];

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*x0[jj]*QWp*eam.phi_PdH(r.norm());

                    // calculate gradient for pair energy
                    app->Fq_Pd0[i] += x0[jj]*QWp*eam.phi_PdH_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_Pd0[i] += x0[jj]*QWp*eam.phi_PdH_deriv(r.norm())/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }
            app->Fq_Pd0[i] *= 1.0/T;
            app->Fsigma_Pd0[i] *= 1.0/T*sqrt(2.0);
            app->Fsigma_Pd0[i] -= 3.0*kB/sigmai;
            sum_sigma += log(sigmai);
        }
        for(int ii=nPd; ii<app->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nPd; // index in the H vectors
            Vec3D &qi = app->q_H0[i]; // by XS
            double &sigmai = app->sigma_H0[i];

            // Gaussian quadrature for embedding energy
            sum_V += x0[i]*eam.F_H(app->recv_H[i]);

            // H-Pd pairs
            for(int j=0; j<(int)HPd[i].size(); j++) {
                int &jj = HPd[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_Pd0[jj];
                double &sigmaj = app->sigma_Pd0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    app->Fq_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_Pd_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_Pd_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[jj])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[jj])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*x0[i]*QWp*eam.phi_PdH(r.norm());

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += QWp*eam.phi_PdH_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.phi_PdH_deriv(r.norm())/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }

            // H-H pairs
            for(int j=0; j<(int)HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_H0[jj];
                double &sigmaj = app->sigma_H0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    app->Fq_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*x0[i]*x0[jj]*QWp*eam.phi_H(r.norm());

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += x0[jj]*QWp*eam.phi_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += x0[jj]*QWp*eam.phi_H_deriv(r.norm())/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }
            app->Fq_H0[i] *= 1.0/T*x0[i];
            app->Fsigma_H0[i] *= 1.0/T*sqrt(2.0)*x0[i];
            app->Fsigma_H0[i] -= 3.0*kB/sigmai*x0[i];
            sum_sigma += log(sigmai)*x0[i];
        }
    } 
    else { //all the local sites are H
        for(int ii=app->first_site; ii<app->last_site_plus_one; ii++) {
            int i = ii - nPd; // index in the H vectors
            Vec3D &qi = app->q_H0[i]; // by XS
            double &sigmai = app->sigma_H0[i];

            // Gaussian quadrature for embedding energy
            sum_V += x0[i]*eam.F_H(app->recv_H[i]);

            // H-Pd pairs
            for(int j=0; j<(int)HPd[i].size(); j++) {
                int &jj = HPd[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_Pd0[jj];
                double &sigmaj = app->sigma_Pd0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    app->Fq_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_Pd_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_Pd_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[jj])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_Pd_deriv(app->recv_Pd[jj])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*x0[i]*QWp*eam.phi_PdH(r.norm());

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += QWp*eam.phi_PdH_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.phi_PdH_deriv(r.norm())/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }

            // H-H pairs
            for(int j=0; j<(int)HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_H0[jj];
                double &sigmaj = app->sigma_H0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    app->Fq_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj])*eam.f_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj])*eam.f_H_deriv(r.norm())/r.norm()*r*QPp[k][0];

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*x0[i]*x0[jj]*QWp*eam.phi_H(r.norm());

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += x0[jj]*QWp*eam.phi_H_deriv(r.norm())/r.norm()*r;
                    app->Fsigma_H0[i] += x0[jj]*QWp*eam.phi_H_deriv(r.norm())/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }
            app->Fq_H0[i] *= 1.0/T*x0[i];
            app->Fsigma_H0[i] *= 1.0/T*sqrt(2.0)*x0[i];
            app->Fsigma_H0[i] -= 3.0*kB/sigmai*x0[i];
            sum_sigma += log(sigmai)*x0[i];
        }
    }
    
    fun_local = 1.0/T*sum_V - 3.0*kB*sum_sigma;

  // ------------------------------------------------------------
  // Revision end.
  // ------------------------------------------------------------

  // copy forces to Grad
  ierr = applicationData2Petsc(app, app->Fq_Pd0, app->Fq_H0, app->Fsigma_Pd0, app->Fsigma_H0, Grad);
  CHKERRQ(ierr);

  //VecView(Grad,PETSC_VIEWER_STDOUT_WORLD);

  // sum up fun from all the processors
  ierr = (PetscErrorCode)MPI_Allreduce((void*)&fun_local, (void*)fun, 1, MPIU_REAL, MPIU_SUM, 
                                       MPI_COMM_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


void Minimizer::updateApplicationContext(double dt_, [[maybe_unused]] vector<Vec3D> *q_Pd0_, [[maybe_unused]] vector<Vec3D> *q_H0_,
                                         [[maybe_unused]] vector<double> *sigma_Pd0_, [[maybe_unused]] vector<double> *sigma_H0_, 
                                         vector<double> *x0_)
{
  application_context->dt = dt_;
  application_context->x0 = x0_;
}


PetscErrorCode Minimizer::applicationData2Petsc(void* ptr, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, 
                                      vector<double> &sigma_Pd0, vector<double> &sigma_H0, Vec &yy)
{
  PetscErrorCode ierr; 
  AppCtx *app = (AppCtx *) ptr;

  int nPd = q_Pd0.size();
  //int nH  = q_H0.size();
  int counter = 0;

  // create locData and loc2glob
  PetscReal locData[app->locDim];
  PetscInt  loc2glob[app->locDim];
  for(int i=0; i<app->locDim; i++) 
    loc2glob[i] = app->first_site*app->numDOFPerSite + i;

  // fill locData
  if(app->first_site<nPd) { // the first site is a Pd
    for(int i=app->first_site; i<min(nPd, app->last_site_plus_one); i++) {
      locData[counter++] = q_Pd0[i][0];
      locData[counter++] = q_Pd0[i][1];
      locData[counter++] = q_Pd0[i][2];
      locData[counter++] = sigma_Pd0[i];
    }
    for(int i=nPd; i<app->last_site_plus_one; i++) {
      locData[counter++] = q_H0[i-nPd][0];
      locData[counter++] = q_H0[i-nPd][1];
      locData[counter++] = q_H0[i-nPd][2];
      locData[counter++] = sigma_H0[i-nPd];
    }
  } else { // all the sites belong to H
    for(int i=app->first_site; i<app->last_site_plus_one; i++) {
      locData[counter++] = q_H0[i-nPd][0];
      locData[counter++] = q_H0[i-nPd][1];
      locData[counter++] = q_H0[i-nPd][2];
      locData[counter++] = sigma_H0[i-nPd];
    }
  }

  if(counter!=app->locDim) {
    PetscPrintf(MPI_COMM_WORLD, "Error in applicationData2Petsc!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // copy locData to Petsc vector 
  ierr = VecSetValues(yy, app->locDim, loc2glob, locData, INSERT_VALUES);
  ierr = VecAssemblyBegin(yy);
  ierr = VecAssemblyEnd(yy);

  return ierr;
}


PetscErrorCode Minimizer::petscData2Application(void* ptr,
                                                Vec &yy, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H, 
                                                vector<double> &sigma_Pd, vector<double> &sigma_H)
{
  PetscErrorCode ierr;
  AppCtx *app = (AppCtx *)ptr;

  PetscScalar *locCopy;
  ierr = VecGetArray(yy, &locCopy);

  for(int i=0; i<app->globDim; i++)
    (app->alldof)[i] = 0.0;

//  double recv_alldof[app->globDim];
    double *recv_alldof;
    recv_alldof = new double[app->globDim];

  for(int i=0; i<app->globDim; i++)
    recv_alldof[i] = 0.0;
 
  for(int i=0; i<app->locDim; i++)
    (app->alldof)[app->first_site*app->numDOFPerSite+i] = locCopy[i];

  // collect data
  MPI_Allreduce(app->alldof, recv_alldof, app->globDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  ierr = VecRestoreArray(yy, &locCopy);

  // copy recv_alldof to application data 
  int counter = 0;
  for(int i=0; i<(int)q_Pd.size(); i++) { //Pd sites
    for(int j=0; j<3; j++)
      q_Pd[i][j] = recv_alldof[counter++];
    sigma_Pd[i] = recv_alldof[counter++];
  }
  for(int i=0; i<(int)q_H.size(); i++) {
    for(int j=0; j<3; j++)
      q_H[i][j] = recv_alldof[counter++];
    sigma_H[i] = recv_alldof[counter++];
  }
  if(counter!=app->globDim) {
    PetscPrintf(MPI_COMM_WORLD, "Error in petscData2Application!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  delete[] recv_alldof;

  return ierr;
}


PetscErrorCode Minimizer::initialize(Input *input, vector<vector<int> > *PdPd, 
                                     vector<vector<int> > *PdH, vector<vector<int> > *HPd, 
                                     vector<vector<int> > *HH, vector<Vec3D> &q_Pd0, 
                                     vector<Vec3D> &q_H0, vector<double> &sigma_Pd0, 
                                     vector<double> &sigma_H0, vector<double> &x0,
                                     vector<vector<Vec3D> > *QPp, double &QWp, 
                                     vector<double> *coe_Pd, vector<double> *coe_H) 
{
  PetscErrorCode ierr;

  // split the solution vector over all the CPU cores
  numDOFPerSite = 4;

  int globSite = q_Pd0.size() + q_H0.size();
  int locSite = PETSC_DECIDE;
  PetscSplitOwnership(comm, &locSite, &globSite); //calculates locSite
  locDim = numDOFPerSite*locSite;
  globDim = PETSC_DECIDE;
  PetscSplitOwnership(comm, &locDim, &globDim); //calculates globDim
 
  // create the solution vector
  ierr = VecCreateMPI(comm, locDim, globDim, &yy); CHKERRQ(ierr);
  int first_dof, last_dof_plus_one;
  ierr = VecSet(yy, 0.0); CHKERRQ(ierr);//initialize 0.0

//  cout << globSite<< " " << locSite << " " << globDim << " " << locDim << endl;

  // create the bound vectors
  ierr = VecDuplicate(yy, &Lb); CHKERRQ(ierr);
  ierr = VecDuplicate(yy, &Ub); CHKERRQ(ierr);
  ierr = VecSet(Lb, 0.0); CHKERRQ(ierr); //initialize 0.0
  ierr = VecSet(Ub, 0.0); CHKERRQ(ierr); //initialize 0.0

  // create a "global" data array
  alldof = new double[globDim];

  // get local site range
  VecGetOwnershipRange(yy, &first_dof, &last_dof_plus_one);
  first_site = first_dof/numDOFPerSite;
  last_site_plus_one = last_dof_plus_one/numDOFPerSite;

  //TODO: the above calculation assumes there are NO ghost sites in the domain; may need to be 
  //      modified. Ghost sites should be included in q_Pd0, ...; but excluded from globDim

  // create application context
  application_context = new AppCtx(input, PdPd, PdH, HPd, HH, &q_Pd0, &q_H0, &sigma_Pd0, &sigma_H0, 
                                   &x0, coe_Pd, coe_H, numDOFPerSite, globDim, locDim, first_site, last_site_plus_one,
                                   alldof, QPp, QWp);

  // create Tao 
  ierr = TaoCreate(PETSC_COMM_WORLD,&tao);CHKERRQ(ierr);
  ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr); //TODO: let's start with quasi-Newton
  ierr = TaoSetObjectiveAndGradientRoutine(tao,Minimizer::formFunctionAndGradient,application_context);
  CHKERRQ(ierr);
  // set bounds
  ierr = TaoSetVariableBoundsRoutine(tao,Minimizer::formVariableBounds,application_context); CHKERRQ(ierr);

  ierr = TaoSetFromOptions(tao);CHKERRQ(ierr); //apply command-line options (if any)

  ierr = PetscPrintf(MPI_COMM_WORLD,"- Initialized Petsc/Tao for optimization.\n");

  return ierr;
}

PetscErrorCode Minimizer::minimizeFreeEntropy(double dt, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0,
                          vector<double> &sigma_Pd0, vector<double> &sigma_H0, vector<double> &x0,
                          vector<Vec3D> &q_Pd, vector<Vec3D> &q_H, vector<double> &sigma_Pd,
                          vector<double> &sigma_H, double &free_entropy)
{
  PetscErrorCode ierr;

  // TODO: I don't think we need to update/rebuild Tao -- but I could be wrong...

  // set up yy 
  ierr = applicationData2Petsc(application_context, q_Pd0, q_H0, sigma_Pd0, sigma_H0, yy);CHKERRQ(ierr);

//  VecView(yy,PETSC_VIEWER_STDOUT_WORLD);

  // update context
  updateApplicationContext(dt, &q_Pd0, &q_H0, &sigma_Pd0, &sigma_H0, &x0);

  // set initial guess (yy)
  ierr = TaoSetInitialVector(tao, yy);CHKERRQ(ierr);

  // SOLVE THE OPTIMIZATION PROBLEM --> yy is updated to store the solution
  ierr = TaoSolve(tao);CHKERRQ(ierr);

  /* Get termination information */
  TaoConvergedReason reason;
  PetscReal fmin = 0.0;
  
  ierr = TaoGetSolutionStatus(tao, NULL, &fmin, NULL, NULL, NULL, &reason);
  free_entropy = (double)fmin;

  if (reason <= 0) {
    ierr = PetscPrintf(MPI_COMM_WORLD,"  WARNING: TAO did not converge!\n");CHKERRQ(ierr);
  }

  // update q_Pd, q_H, sigma_Pd, sigma_H 
  ierr = petscData2Application(application_context, yy, q_Pd, q_H, sigma_Pd, sigma_H);CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode Minimizer::calculateChemicalPotential(vector<double> &x0, double gammabd, vector<double> &gamma, vector<int> &full_H)
{
    PetscErrorCode ierr = 0;    
    EAM eam;
    const double kB = 8.6173324e-5;
    double &T = application_context->input->file.T0;

    int nPd = application_context->q_Pd0.size();
    int nH  = application_context->q_H0.size();
    vector<vector<int> > &PdPd = *(application_context->PdPd); //create reference
    vector<vector<int> > &PdH = *(application_context->PdH);
    vector<vector<int> > &HPd = *(application_context->HPd);
    vector<vector<int> > &HH = *(application_context->HH);
    vector<vector<Vec3D> > &QPp = *(application_context->QPp);
    double &QWp = application_context->QWp;

    vector<double> &coe_Pd = *(application_context->coe_Pd);
    vector<double> &coe_H = *(application_context->coe_H);

    // initiallization for electron density calculation
    for(int i=0; i<nPd; i++)
        application_context->send_Pd[i] = 0.0;
    for(int i=0; i<nH; i++)
        application_context->send_H[i] = 0.0;

    // ------------------------------------------------------------
    // Parallelization I: calculate electron density of each site at each Gaussian point.
    // ------------------------------------------------------------
    if(application_context->first_site<nPd) { // the first local site is a Pd
        for(int i=application_context->first_site; i<min(nPd, application_context->last_site_plus_one); i++) {
            application_context->send_Pd[i] = electronDensity('P', i, application_context->input->file.rc, PdPd, PdH, HPd, HH, application_context->q_Pd0, application_context->q_H0, application_context->sigma_Pd0, application_context->sigma_H0, x0, QPp, QWp);
        }
        for(int ii=nPd; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nPd; // index in the H vectors
            application_context->send_H[i] = electronDensity('H', i, application_context->input->file.rc, PdPd, PdH, HPd, HH, application_context->q_Pd0, application_context->q_H0, application_context->sigma_Pd0, application_context->sigma_H0, x0, QPp, QWp);
       }
    }
    else { //all the local sites are H
        for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
            int i = ii - nPd; // index in the H vectors
            application_context->send_H[i] = electronDensity('H', i, application_context->input->file.rc, PdPd, PdH, HPd, HH, application_context->q_Pd0, application_context->q_H0, application_context->sigma_Pd0, application_context->sigma_H0, x0, QPp, QWp);
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_Pd, application_context->recv_Pd, nPd, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Initialize chemical potential.
    for(int i=0; i<nH; i++)
        application_context->send_H[i] = 0.0;   

    // ------------------------------------------------------------
    // Parallelization II: calculate chemical potential at H sites.
    // ------------------------------------------------------------ 
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    int MPI_rank, MPI_size;
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_size);
    
    int locSize; // local size for each cores except the last one
    int lastSize; // local size for the last one
    if(nH%MPI_size == 0) {
        locSize = nH/MPI_size;
        lastSize = nH/MPI_size;
    } else {
        locSize = nH/(MPI_size-1);
        lastSize = nH%(MPI_size-1);
    }    

    if((MPI_rank+1) != MPI_size) { // not the last core
        for(int ii=0; ii<locSize; ii++) {
            int i = MPI_rank*locSize + ii;
            if(full_H[i]==1) continue; // skip the site which is fully occupied by H
            Vec3D &qi = application_context->q_H0[i]; // by XS
            double &sigmai = application_context->sigma_H0[i];
            double sum_Vx = 0.0; // the derivative of potential energy with respect to H fraction on site i

            // Gaussian quadrature for embedding energy
            sum_Vx += eam.F_H(application_context->recv_H[i]);
            
            // correction term
            double tmp = coe_H[1] * HPd[i].size() * HPd[i].size();
            for(int j=0; j<(int)HH[i].size(); j++) {
                tmp += 2.0 * coe_H[2] * HPd[i].size() * x0[HH[i][j]];
                for(int l=0; l<(int)HH[i].size(); l++)
                    tmp += coe_H[3] * x0[HH[i][j]] * x0[HH[i][l]];
            }
            sum_Vx += 0.5 * T * coe_H[0] * tmp;
 
            // H-Pd pairs
            for(int j=0; j<(int)HPd[i].size(); j++) {
                int &jj = HPd[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_Pd0[jj];
                double &sigmaj = application_context->sigma_Pd0[jj];
                
                for(int k=0; k<(int)QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    sum_Vx += QWp*eam.F_Pd_deriv(application_context->recv_Pd[jj])*eam.f_H(r.norm());

                    // Gaussian quadrature for pair energy
                    sum_Vx += QWp*eam.phi_PdH(r.norm());
                }
            
                // correction term
                sum_Vx += T * coe_Pd[0] * coe_Pd[2] * PdPd[jj].size();
                for(int l=0; l<(int)PdH[jj].size(); l++)
                    sum_Vx += T * coe_Pd[0] * coe_Pd[3] * x0[PdH[jj][l]];

            }

            // H-H pairs
            for(int j=0; j<(int)HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_H0[jj];
                double &sigmaj = application_context->sigma_H0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    sum_Vx += x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj])*eam.f_H(r.norm());

                    // Gaussian quadrature for pair energy
                    sum_Vx += x0[jj]*QWp*eam.phi_H(r.norm());
                }

                // correction term
                sum_Vx += x0[jj] * T * coe_H[0] * coe_H[2] * HPd[jj].size();
                for(int l=0; l<(int)HH[jj].size(); l++)
                    sum_Vx += x0[jj] * T * coe_H[0] * coe_H[3] * x0[HH[jj][l]];
            }

            if (x0[i]>0.0 && x0[i]<1.0)
                application_context->send_H[i] = log(x0[i]/(1.0-x0[i])) + sum_Vx/(kB*T);
        }
    } else { // the last core
        for(int ii=0; ii<lastSize; ii++) {
            int i = MPI_rank*locSize + ii;
            if(full_H[i]==1) continue; // skip the site which is fully occupied by H
            Vec3D &qi = application_context->q_H0[i]; // by XS
            double &sigmai = application_context->sigma_H0[i];
            double sum_Vx = 0.0; // the derivative of potential energy with respect to H fraction on site i

            // Gaussian quadrature for embedding energy
            sum_Vx += eam.F_H(application_context->recv_H[i]);

            // correction term
            double tmp = coe_H[1] * HPd[i].size() * HPd[i].size();
            for(int j=0; j<(int)HH[i].size(); j++) {
                tmp += 2.0 * coe_H[2] * HPd[i].size() * x0[HH[i][j]];
                for(int l=0; l<(int)HH[i].size(); l++)
                    tmp += coe_H[3] * x0[HH[i][j]] * x0[HH[i][l]];
            }
            sum_Vx += 0.5 * T * coe_H[0] * tmp;

            // H-Pd pairs
            for(int j=0; j<(int)HPd[i].size(); j++) {
                int &jj = HPd[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_Pd0[jj];
                double &sigmaj = application_context->sigma_Pd0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    sum_Vx += QWp*eam.F_Pd_deriv(application_context->recv_Pd[jj])*eam.f_H(r.norm());

                    // Gaussian quadrature for pair energy
                    sum_Vx += QWp*eam.phi_PdH(r.norm());
                }

                // correction term
                sum_Vx += T * coe_Pd[0] * coe_Pd[2] * PdPd[jj].size();
                for(int l=0; l<(int)PdH[jj].size(); l++)
                    sum_Vx += T * coe_Pd[0] * coe_Pd[3] * x0[PdH[jj][l]];
            }

            // H-H pairs
            for(int j=0; j<(int)HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_H0[jj];
                double &sigmaj = application_context->sigma_H0[jj];

                for(int k=0; k<(int)QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    sum_Vx += x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj])*eam.f_H(r.norm());

                    // Gaussian quadrature for pair energy
                    sum_Vx += x0[jj]*QWp*eam.phi_H(r.norm());
                }

                // correction term
                sum_Vx += x0[jj] * T * coe_H[0] * coe_H[2] * HPd[jj].size();
                for(int l=0; l<(int)HH[jj].size(); l++)
                    sum_Vx += x0[jj] * T * coe_H[0] * coe_H[3] * x0[HH[jj][l]];
            }

            if (x0[i]>0.0 && x0[i]<1.0)
                application_context->send_H[i] = log(x0[i]/(1.0-x0[i])) + sum_Vx/(kB*T);
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    CHKERRQ(ierr);
    for (int i=0; i<nH; i++) {
        if(full_H[i]==1) continue; // skip the site which is fully occupied by H
        gamma[i] = application_context->recv_H[i] + 1.5;
        if(gamma[i]>gammabd) {
            gamma[i] = gammabd;
            if (x0[i]>0.9999) full_H[i] = 1; // find the fully occupied H boundary sites          
        }
    }
    return ierr;
}

PetscErrorCode Minimizer::formVariableBounds([[maybe_unused]] Tao tao, Vec Lb, Vec Ub, void* ptr)
{
    PetscErrorCode ierr = 0;
    AppCtx *app = (AppCtx *) ptr;

    vector<Vec3D> Ub_q_Pd(app->q_Pd0), Lb_q_Pd(app->q_Pd0), Ub_q_H(app->q_H0), Lb_q_H(app->q_H0);
    vector<double> Ub_sigma_Pd(app->sigma_Pd0), Lb_sigma_Pd(app->sigma_Pd0), Ub_sigma_H(app->sigma_H0), Lb_sigma_H(app->sigma_H0);

    for (int i=0; i<(int)Ub_q_Pd.size(); i++) {
        Lb_q_Pd[i] = Vec3D(PETSC_NINFINITY,PETSC_NINFINITY,PETSC_NINFINITY);
        Ub_q_Pd[i] = Vec3D(PETSC_INFINITY,PETSC_INFINITY,PETSC_INFINITY);
        Lb_sigma_Pd[i] = 1.0e-6;
        Ub_sigma_Pd[i] = PETSC_INFINITY;
    }

    for (int i=0; i<(int)Ub_q_H.size(); i++) {
        Lb_q_H[i] = Vec3D(PETSC_NINFINITY,PETSC_NINFINITY,PETSC_NINFINITY);
        Ub_q_H[i] = Vec3D(PETSC_INFINITY,PETSC_INFINITY,PETSC_INFINITY);
        Lb_sigma_H[i] = 1.0e-6;
        Ub_sigma_H[i] = PETSC_INFINITY;
    }

    ierr = applicationData2Petsc(app, Ub_q_Pd, Ub_q_H, Ub_sigma_Pd, Ub_sigma_H, Ub); CHKERRQ(ierr);
    ierr = applicationData2Petsc(app, Lb_q_Pd, Lb_q_H, Lb_sigma_Pd, Lb_sigma_H, Lb); CHKERRQ(ierr);
    
    return ierr;
}


PetscErrorCode Minimizer::findLocalNeighbors(vector<vector<int> > &PdPd_ext, vector<vector<int> > &PdH_ext,
                  vector<vector<int> > &HPd_ext, vector<vector<int> > &HH_ext,
                  vector<vector<int> > &PdPd, vector<vector<int> > &PdH,
                  vector<vector<int> > &HPd, vector<vector<int> > &HH,
                  vector<vector<int> > &PdPdnum, vector<vector<int> > &HH1, int &maxNeib)

{
    PetscErrorCode ierr = 0;

    int maxNeibTmp = 100;    
    int nPd = application_context->q_Pd0.size();
    int nH  = application_context->q_H0.size();

    vector<Vec3D> &q_Pd0 = application_context->q_Pd0;
    vector<Vec3D> &q_H0 = application_context->q_H0;
    double rc = application_context->input->file.rc;
    double a_Pd = application_context->input->file.a_Pd; 

    // temporary vectors for MPI
    int *send_Pd1, *recv_Pd1;
    int *send_Pd2, *recv_Pd2;
    int *send_Pd3, *recv_Pd3;
    int *send_H1, *recv_H1;
    int *send_H2, *recv_H2;
    int *send_H3, *recv_H3;
    int send_PdH = 0;

    send_Pd1 = new int[nPd*maxNeibTmp];
    recv_Pd1 = new int[nPd*maxNeibTmp];
    send_Pd2 = new int[nPd*maxNeibTmp];
    recv_Pd2 = new int[nPd*maxNeibTmp];
    send_Pd3 = new int[nPd*3];
    recv_Pd3 = new int[nPd*3];
    send_H1 = new int[nH*maxNeibTmp];
    recv_H1 = new int[nH*maxNeibTmp];
    send_H2 = new int[nH*maxNeibTmp];
    recv_H2 = new int[nH*maxNeibTmp];
    send_H3 = new int[nH*maxNeibTmp];
    recv_H3 = new int[nH*maxNeibTmp];


    // initiallization
    for(int i=0; i<nPd*maxNeibTmp; i++) {
        send_Pd1[i] = 0;
        send_Pd2[i] = 0;
    }

    for(int i=0; i<nPd*3; i++)
        send_Pd3[i] = 0;

    for(int i=0; i<nH*maxNeibTmp; i++) {
        send_H1[i] = 0;
        send_H2[i] = 0;
        send_H3[i] = 0;
    }

    // Find and store neighbors
    maxNeib = 0;
    int maxNei = 500; //max number of neighbors for each atom. 
    double qLocal[3];

    // ------------------------------------------------------------
    // Parallelization 
    // ------------------------------------------------------------
    if(application_context->first_site<nPd) { // the first local site is a Pd
        for(int i=application_context->first_site; i<min(nPd, application_context->last_site_plus_one); i++) {
            // Store local sites in a KD-Tree with K = 3.
            PointIn3D *atoms = new PointIn3D[PdPd_ext[i].size()+PdH_ext[i].size()];
            for(int j=0; j<(int)PdPd_ext[i].size(); j++)
                atoms[j] = PointIn3D(j, q_Pd0[PdPd_ext[i][j]]);
            for(int j=0; j<(int)PdH_ext[i].size(); j++)
                atoms[PdPd_ext[i].size() + j] = PointIn3D(PdPd_ext[i].size() + j, q_H0[PdH_ext[i][j]]);
            KDTree<PointIn3D> atomTree(PdPd_ext[i].size()+PdH_ext[i].size(), atoms); //Note: atoms are re-ordered
            
            int nNeib = 0;
            PointIn3D candidates[maxNei];
            for(int j=0; j<3; j++)
                qLocal[j] = q_Pd0[i][j]; // position of i
            // find neighbors of i within rc --> populates "candidates"
            int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
            // debug
            if(nFound>maxNei) {
                cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei
                    << ". Increase maxNei! " << endl;
                exit(-1);
            }

            int iPd1 = 0;
            int iPd2 = 0;
            // populates PdPd, PdH  
            for(int j=0; j<nFound; j++) {
                int nei = candidates[j].pid(); //one neighbor
                //filter the outside ones -- KDTree works with "pseudo distance"
                if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                   (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                   (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) > rc*rc)
                    continue;
                nNeib++;

                if(nei<(int)PdPd_ext[i].size()) {// this neighbor is Pd
                    send_Pd1[i*maxNeibTmp+iPd1] = PdPd_ext[i][nei] + 1;
                    iPd1++;
                    if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                       (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                       (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) < a_Pd*a_Pd*0.75)
                        send_Pd3[i*3]++;
                    else if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                       (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                       (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) > a_Pd*a_Pd*0.75 &&
                       (candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                       (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                       (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) < a_Pd*a_Pd*1.25)
                        send_Pd3[i*3+1]++;
                    else
                        send_Pd3[i*3+2]++;
                } else { // this neighbor is H
                    send_Pd2[i*maxNeibTmp+iPd2] = PdH_ext[i][nei-PdPd_ext[i].size()]+1;
                    iPd2++;
                }
            }
            if(nNeib>send_PdH) send_PdH = nNeib;
            delete[] atoms;
        }
        for(int ii=nPd; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nPd; // index in the H vectors
            // Store local sites in a KD-Tree with K = 3.
            PointIn3D *atoms = new PointIn3D[HPd_ext[i].size()+HH_ext[i].size()];
            for(int j=0; j<(int)HPd_ext[i].size(); j++)
                atoms[j] = PointIn3D(j, q_Pd0[HPd_ext[i][j]]);
            for(int j=0; j<(int)HH_ext[i].size(); j++)
                atoms[HPd_ext[i].size() + j] = PointIn3D(HPd_ext[i].size() + j, q_H0[HH_ext[i][j]]);
            KDTree<PointIn3D> atomTree(HPd_ext[i].size()+HH_ext[i].size(), atoms); //Note: atoms are re-ordered

            int nNeib = 0;
            PointIn3D candidates[maxNei];
            for(int j=0; j<3; j++)
                qLocal[j] = q_H0[i][j]; // position of i
            // find neighbors of i within rc --> populates "candidates"
            int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
            // debug
            if(nFound>maxNei) {
                cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei
                    << ". Increase maxNei! " << endl;
                exit(-1);
            }

            int iH1 = 0;
            int iH2 = 0;
            int iH3 = 0;
            // populates HPd, HH, HH1
            for(int j=0; j<nFound; j++) {
                int nei = candidates[j].pid(); //one neighbor
                //filter the outside ones -- KDTree works with "pseudo distance"
                if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                   (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                   (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) > rc*rc)
                  continue;            

                nNeib++;

                if(nei<(int)HPd_ext[i].size()) {//this neighbor is Pd
                    send_H2[i*maxNeibTmp+iH2] = HPd_ext[i][nei]+1;
                    iH2++;
                } else {//this neighbor is H
                    send_H1[i*maxNeibTmp+iH1] = HH_ext[i][nei-HPd_ext[i].size()]+1;
                    iH1++;
                    if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                       (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                       (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) < a_Pd*a_Pd*0.75){
                        send_H3[i*maxNeibTmp+iH3] = HH_ext[i][nei-HPd_ext[i].size()]+1;
                        iH3++;
                    } 
                }
            }
            if(nNeib>send_PdH) send_PdH = nNeib;
            delete[] atoms;
        }
    }
    else { //all the local sites are H
        for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
            int i = ii - nPd; // index in the H vectors
            // Store local sites in a KD-Tree with K = 3.
            PointIn3D *atoms = new PointIn3D[HPd_ext[i].size()+HH_ext[i].size()];
            for(int j=0; j<(int)HPd_ext[i].size(); j++)
                atoms[j] = PointIn3D(j, q_Pd0[HPd_ext[i][j]]);
            for(int j=0; j<(int)HH_ext[i].size(); j++)
                atoms[HPd_ext[i].size() + j] = PointIn3D(HPd_ext[i].size() + j, q_H0[HH_ext[i][j]]);
            KDTree<PointIn3D> atomTree(HPd_ext[i].size()+HH_ext[i].size(), atoms); //Note: atoms are re-ordered

            int nNeib = 0;
            PointIn3D candidates[maxNei];
            for(int j=0; j<3; j++)
                qLocal[j] = q_H0[i][j]; // position of i
            // find neighbors of i within rc --> populates "candidates"
            int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
            // debug
            if(nFound>maxNei) {
                cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei
                    << ". Increase maxNei! " << endl;
                exit(-1);
            }

            int iH1 = 0;
            int iH2 = 0;
            int iH3 = 0;
            // populates HPd, HH, HH1
            for(int j=0; j<nFound; j++) {
                int nei = candidates[j].pid(); //one neighbor
                //filter the outside ones -- KDTree works with "pseudo distance"
                if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                   (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                   (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) > rc*rc)
                  continue;

                nNeib++;

                if(nei<(int)HPd_ext[i].size()) {//this neighbor is Pd
                    send_H2[i*maxNeibTmp+iH2] = HPd_ext[i][nei]+1;
                    iH2++;
                } else {//this neighbor is H
                    send_H1[i*maxNeibTmp+iH1] = HH_ext[i][nei-HPd_ext[i].size()]+1;
                    iH1++;
                    if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
                       (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
                       (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) < a_Pd*a_Pd*0.75){
                        send_H3[i*maxNeibTmp+iH3] = HH_ext[i][nei-HPd_ext[i].size()]+1;
                        iH3++;
                    }
                }
            }
            if(nNeib>send_PdH) send_PdH = nNeib;
            delete[] atoms;
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(send_Pd1, recv_Pd1, nPd*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_Pd2, recv_Pd2, nPd*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_Pd3, recv_Pd3, nPd*3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_H1, recv_H1, nH*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_H2, recv_H2, nH*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_H3, recv_H3, nH*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    ierr = (PetscErrorCode)MPI_Allreduce(&send_PdH, &maxNeib, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    delete[] send_Pd1;
    delete[] send_Pd2;
    delete[] send_Pd3;
    delete[] send_H1;
    delete[] send_H2;
    delete[] send_H3;

    if(!(PdPd.empty() && PdH.empty() && HPd.empty() && HH.empty() && PdPdnum.empty() && HH1.empty())) {
        PdPd.clear();
        PdH.clear();
        HPd.clear();
        HH.clear();
        PdPdnum.clear();
        HH1.clear();
    }

    // Initialize neighbor lists
    vector<int> tmp; //an empty vector
    tmp.reserve(maxNeibTmp);
    PdPd.assign(nPd, tmp); //Pd neighbors of Pd
    PdH.assign(nPd, tmp); //H neighbors of Pd
    HPd.assign(nH, tmp); //Pd neighbors of H
    HH.assign(nH, tmp); //H neighbors of H
    HH1.assign(nH, tmp); //first H neighbor shell of H

    vector<int> tmp1; //a zero vector
    tmp1.assign(3, 0);
    PdPdnum.assign(nPd, tmp1); // number of Pd neighbors on different shells;


    for(int i=0; i<nPd; i++) {
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_Pd1[i*maxNeibTmp+j]==0) break;
            PdPd[i].push_back(recv_Pd1[i*maxNeibTmp+j]-1);
        }
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_Pd2[i*maxNeibTmp+j]==0) break;
            PdH[i].push_back(recv_Pd2[i*maxNeibTmp+j]-1);
        }
        for (int j=0; j<3; j++) {
            PdPdnum[i][j] = recv_Pd3[i*3+j];
        }
    }

    for(int i=0; i<nH; i++) {
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_H1[i*maxNeibTmp+j]==0) break;
            HH[i].push_back(recv_H1[i*maxNeibTmp+j]-1);
        }
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_H2[i*maxNeibTmp+j]==0) break;
            HPd[i].push_back(recv_H2[i*maxNeibTmp+j]-1);
        }
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_H3[i*maxNeibTmp+j]==0) break;
            HH1[i].push_back(recv_H3[i*maxNeibTmp+j]-1);
        }
    }

    delete[] recv_Pd1;
    delete[] recv_Pd2;
    delete[] recv_Pd3;
    delete[] recv_H1;
    delete[] recv_H2;
    delete[] recv_H3;

    return ierr;
}

