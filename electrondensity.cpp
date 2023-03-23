/*
- calculate averaged electron density of one Pd or H atom
*/
#include <iostream>
#include <vector>
#include <math.h>
#include "Vector3D.h"
#include "potential.h"
using namespace std;

EAM eam;
double electronDensity(char atype, int n, [[maybe_unused]] double rc, vector<vector<int> > &PdPd, 
                    vector<vector<int> > &PdH, vector<vector<int> > &HPd, 
                    vector<vector<int> > &HH, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H, 
                    vector<double> &sigma_Pd, vector<double> &sigma_H, vector<double> &x, 
                    vector<vector<Vec3D> > &QPp, double &QWp)
{
    double rho = 0.0;

    if (atype=='P') {
        for(int i=0; i<(int)PdPd[n].size(); i++) {
            for(int k=0; k<(int)QPp.size(); k++) {
                Vec3D r = q_Pd[n]+sqrt(2.0)*sigma_Pd[n]*QPp[k][0]-q_Pd[PdPd[n][i]]-sqrt(2.0)*sigma_Pd[PdPd[n][i]]*QPp[k][1];
                rho += eam.f_Pd(r.norm())*QWp;    
            }
        }

        for(int i=0; i<(int)PdH[n].size(); i++) {
            for(int k=0; k<(int)QPp.size(); k++) {
                Vec3D r = q_Pd[n]+sqrt(2.0)*sigma_Pd[n]*QPp[k][0]-q_H[PdH[n][i]]-sqrt(2.0)*sigma_H[PdH[n][i]]*QPp[k][1];
                rho += x[PdH[n][i]]*eam.f_H(r.norm())*QWp;
            }
        }
    }
    else if (atype=='H') {
        for(int i=0; i<(int)HPd[n].size(); i++) {
            for(int k=0; k<(int)QPp.size(); k++) {
                Vec3D r = q_H[n]+sqrt(2.0)*sigma_H[n]*QPp[k][0]-q_Pd[HPd[n][i]]-sqrt(2.0)*sigma_Pd[HPd[n][i]]*QPp[k][1];
                rho += eam.f_Pd(r.norm())*QWp;
            }
        }

        for(int i=0; i<(int)HH[n].size(); i++) {
            for(int k=0; k<(int)QPp.size(); k++) {
                Vec3D r = q_H[n]+sqrt(2.0)*sigma_H[n]*QPp[k][0]-q_H[HH[n][i]]-sqrt(2.0)*sigma_H[HH[n][i]]*QPp[k][1];
                rho += x[HH[n][i]]*eam.f_H(r.norm())*QWp;
            }
        }
    }
    else {
        cerr << "Error: First input in electronDensity must be either 'P' or 'H'." << endl; cout.flush();
        exit(-1);
    }

    return rho;
}


