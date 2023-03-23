#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>
#include "Vector3D.h"
using namespace std;

int calculateMacroResults([[maybe_unused]] vector<vector<int> > &PdPd, [[maybe_unused]] vector<vector<int> > &PdPdnum,
                          [[maybe_unused]] vector<Vec3D> &q_Pd, vector<double> &x, [[maybe_unused]] vector<int> &PdSurf,
                          vector<int> &HSurf, double &xH, double &xHinterior)
{
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    int MPI_rank, MPI_size;
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_size);

    xH = 0.0;
    xHinterior = 0.0;
    //int nPd = q_Pd.size();
    int nH  = x.size();
    double *send, *recv;
    double *send3, *recv3;
    send = new double[1];
    recv = new double[1];
    send3 = new double[1];
    recv3 = new double[1];

    // calculate H mole ratio
    send[0] = 0.0;
    recv[0] = 0.0;
    send3[0] = 0.0;
    recv3[0] = 0.0;

    int locSize = 0; // local size for each cores except the last one
    int lastSize = 0; // local size for the last one

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
            send[0] += x[i];

            vector<int>::iterator it; 
            it = find(HSurf.begin(),HSurf.end(),i);
            if (it==HSurf.end()) send3[0] += x[i];
        }
    } else { // the last core
        for(int ii=0; ii<lastSize; ii++) {
            int i = MPI_rank*locSize + ii;
            send[0] += x[i];
            
            vector<int>::iterator it;                 
            it = find(HSurf.begin(),HSurf.end(),i);
            if (it==HSurf.end()) send3[0] += x[i];
        }
    }

    MPI_Barrier(comm);
    MPI_Allreduce(send, recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(send3, recv3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    xH = recv[0]/nH;
    xHinterior = recv3[0]/(nH-HSurf.size());

    delete[] send;
    delete[] recv;
    delete[] send3;
    delete[] recv3;
    
    return 0;
}
