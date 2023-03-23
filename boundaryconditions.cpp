#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>
#include "Vector3D.h"
#include "input.h"
using namespace std;
/*
int enforceBoundaryConditions(Input &input, double dt, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0,
                              vector<double> &sigma_Pd0, vector<double> &sigma_H0, vector<double> &x0, 
                              vector<double> &gamma, vector<int> &HSubsurf, double &gammabd, vector<int> &full_H)
{
  // TODO: Appropriate boundary conditions need to be specified in "Input" and enforced here.

    const double xUb = 1.0-1.0e-16; // The upper bound of hydrogen atomic fraction in the sample
    const double xLb = 1.0e-3; // The lower bound of hydrogen atomic fraction in the sample

    for(int i=0; i<HSubsurf.size(); i++) {
        x0[HSubsurf[i]] = input.file.x_subsurf;
        if (input.file.fix_mu == 1) {
            full_H[HSubsurf[i]] = 1;
            gamma[HSubsurf[i]] = gammabd;
        }
    }

    return 0;
}
*/

int enforceBoundaryConditions(Input &input, [[maybe_unused]] double dt, [[maybe_unused]] vector<Vec3D> &q_Pd0,
                              [[maybe_unused]] vector<Vec3D> &q_H0,
                              [[maybe_unused]] vector<double> &sigma_Pd0, [[maybe_unused]] vector<double> &sigma_H0, vector<double> &x0, 
                              vector<double> &gamma, vector<int> &HSubsurf, double &gammabd, vector<int> &full_H)
{
  // TODO: Appropriate boundary conditions need to be specified in "Input" and enforced here.

    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    int MPI_rank, MPI_size;
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_size);

    //const double kB = 8.6173324e-5;
    double &xUb = input.file.xUb; // The upper bound of hydrogen atomic fraction in the sample
    //double &xLb = input.file.xLb; // The lower bound of hydrogen atomic fraction in the sample   
 
    for(int i=0; i<(int)HSubsurf.size(); i++) {
        gamma[HSubsurf[i]] = gammabd;
        full_H[HSubsurf[i]] = 1;
        x0[HSubsurf[i]] = xUb;
    }
 
    return 0;
}

