/*
- generate Gaussian points and weight values
*/
#include <iostream>
#include <vector>
#include <math.h>
#include <Vector3D.h>
using namespace std;
/*
void gaussianQuadrature(int dim, vector<vector<double> > &QP, vector<double> &QW)
{   
    if(dim%3){
        cerr << "Error: Quadrature dimension must be a multiple of 3." << endl; cout.flush();
        exit(-1);
    }

    if(!(QP.empty() && QW.empty())) {
        QP.clear();
        QW.clear();
    }
    
    for (int i=0; i<dim; i++){
        vector<double> tmp1, tmp2;
        for (int j=0; j<dim; j++){
            if (i==j){        
                tmp1.push_back(sqrt(dim/2.0));
                tmp2.push_back(-sqrt(dim/2.0));
            }
            else{
                tmp1.push_back(0);
                tmp2.push_back(0);
            }
        }
        QP.push_back(tmp1);
        QP.push_back(tmp2);
    }   

    for (int i=0; i<2*dim; i++)
        QW.push_back(0.5/dim);
}
*/

void gaussianQuadrature(int dim, vector<vector<Vec3D> > &QP)
{   
    if (dim%3) {
        cerr << "Error: Quadrature dimension must be a multiple of 3." << endl; cout.flush();
        exit(-1);
    }

    if (!(QP.empty()))
        QP.clear();
    
    for (int i=0; i<dim; i++) {
        vector<double> tmp1, tmp2;
        for (int j=0; j<dim; j++) {
            if (i==j) {        
                tmp1.push_back(sqrt(dim/2.0));
                tmp2.push_back(-sqrt(dim/2.0));
            }
            else {
                tmp1.push_back(0.0);
                tmp2.push_back(0.0);
            }
        }

        vector<Vec3D> tmp3, tmp4;
        for (int j=0; j<dim/3; j++) {
            tmp3.push_back(Vec3D(tmp1[j*3], tmp1[j*3+1], tmp1[j*3+2]));
            tmp4.push_back(Vec3D(tmp2[j*3], tmp2[j*3+1], tmp2[j*3+2]));
        }        

        QP.push_back(tmp3);
        QP.push_back(tmp4);
    }   
}


void gaussianPointsInRij(int dim, vector<Vec3D> &QPi)
{
    if (!QPi.empty()) { 
        QPi.clear();
    }

    for (int i=0; i<3; i++) {
        vector<double> tmp1, tmp2;
        for (int j=0; j<3; j++) {
            if (i==j) {
                tmp1.push_back(sqrt(dim/2.0));
                tmp2.push_back(-sqrt(dim/2.0));
            }
            else {
                tmp1.push_back(0.0);
                tmp2.push_back(0.0);
            }
        }

        QPi.push_back(Vec3D(tmp1[0], tmp1[1], tmp1[2]));
        QPi.push_back(Vec3D(tmp2[0], tmp2[1], tmp2[2]));
    }
}

