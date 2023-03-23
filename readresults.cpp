#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <Vector3D.h>
#include <input.h>
using namespace std;

// read the results for restart
void readResults(Input &input, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, vector<double> &sigma_Pd0,
                   vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, vector<int> &full_H)
{
    char full_filename_base[256];    
    sprintf(full_filename_base, "%s/%s", input.file.foldername, input.file.filename_base);

    char full_fname[512];
    
    sprintf(full_fname, "%s.%d.dat", full_filename_base, input.file.restart_file_num);
    ifstream file(full_fname, ios::in);

    if(!file){
        cerr << "ERREOR! The input solution file for restart does not exit!" <<  endl;
        exit(-1);
    }

    // read the solution file for restart line by line
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');

    double tmp1;
    for(int i=0; i<(int)q_Pd0.size(); i++)
        file >> tmp1 >> tmp1 >>  q_Pd0[i][0] >> q_Pd0[i][1] >> q_Pd0[i][2] >> sigma_Pd0[i] >> tmp1 >> tmp1 >> tmp1;
    for(int i=0; i<(int)q_H0.size(); i++)
        file >> tmp1 >> tmp1 >> q_H0[i][0] >> q_H0[i][1] >> q_H0[i][2] >> sigma_H0[i] >> x0[i] >> gamma[i] >> full_H[i];

    file.close();

    double &xUb = input.file.xUb; // The upper bound of hydrogen atomic fraction in the sample
    double &xLb = input.file.xLb; // The lower bound of hydrogen atomic fraction in the sample

    for(int i=0; i<(int)x0.size(); i++) {
        if(x0[i]<xLb) x0[i] = xLb;
        else if(x0[i]>xUb) x0[i] = xUb;
    }

    return; 
}

