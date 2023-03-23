#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>
#include <Vector3D.h>
#include <input.h>
using namespace std;

int integrateDiffusion(Input &input, double dt, vector<vector<int> > &HH1,
                       vector<double> &gamma, vector<double> &x0, vector<double> &x, vector<int> &full_H)
{
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    int MPI_rank, MPI_size;
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_size);

    const double kB = 8.6173324e-5;
    double &xUb = input.file.xUb; // The upper bound of hydrogen atomic fraction in the sample
    double &xLb = input.file.xLb; // The lower bound of hydrogen atomic fraction in the sample
    
    double &T = input.file.T0; // temperature
    double &v = input.file.v; // attempt frequency
    double &Qm = input.file.Qm; // activation energy

    double prefactor = v * exp(-Qm/(kB*T));

    int nH  = x0.size();
    int locSize; // local size for each cores except the last one
    int lastSize; // local size for the last one
    
    double *send_x, *recv_x;
    send_x = new double[nH];
    recv_x = new double[nH];

    for(int i=0; i<nH; i++)
        send_x[i] = 0.0;

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
            double x_deriv = 0.0;
            double fi = gamma[i] - log(x0[i]/(1.0-x0[i])); 
            for(int j=0; j<(int)HH1[i].size(); j++) { // restrict hydrogen diffusion to the nearest first shell
                double fj = gamma[HH1[i][j]] - log(x0[HH1[i][j]]/(1.0-x0[HH1[i][j]]));
                x_deriv += x0[HH1[i][j]]*(1.0-x0[i])*exp(-(fi-fj)/2.0) - x0[i]*(1.0-x0[HH1[i][j]])*exp((fi-fj)/2.0);
            } 
            send_x[i] = x0[i] + prefactor*x_deriv*dt;
        }
    } else { // the last core
        for(int ii=0; ii<lastSize; ii++) {
            int i = MPI_rank*locSize + ii;
            if(full_H[i]==1) continue; // skip the site which is fully occupied by H
            double x_deriv = 0.0;
            double fi = gamma[i] - log(x0[i]/(1.0-x0[i]));
            for(int j=0; j<(int)HH1[i].size(); j++) { // restrict hydrogen diffusion to the nearest first shell
                double fj = gamma[HH1[i][j]] - log(x0[HH1[i][j]]/(1.0-x0[HH1[i][j]]));
                x_deriv += x0[HH1[i][j]]*(1.0-x0[i])*exp(-(fi-fj)/2.0) - x0[i]*(1.0-x0[HH1[i][j]])*exp((fi-fj)/2.0);
            }
            send_x[i] = x0[i] + prefactor*x_deriv*dt;
        }
    }

    MPI_Barrier(comm);
    MPI_Allreduce(send_x, recv_x, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
    for(int i=0; i<nH; i++) {
        if(full_H[i]==1) continue; // skip the site which is fully occupied by H
        x[i] = recv_x[i];
        if(x[i]<xLb) x[i] = xLb;
        else if(x[i]>xUb) x[i] = xUb;
    }

    delete[] send_x;
    delete[] recv_x;
    
    return 0;
}
