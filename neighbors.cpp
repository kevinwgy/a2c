#include <vector>
#include <iostream>
#include <Vector3D.h>
#include <KDTree.h>
#include <input.h>
using namespace std;

int findNeighbors([[maybe_unused]] double a_Pd, double rc, vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0,
                  vector<vector<int> > &PdPd, vector<vector<int> > &PdH,
                  vector<vector<int> > &HPd, vector<vector<int> > &HH)
{
  if(!(PdPd.empty() && PdH.empty() && HPd.empty() && HH.empty())) {
    cout << "WARNING: Some neighbor list(s) is/are not empty. Clearing all the lists." << endl;
    PdPd.clear();
    PdH.clear();
    HPd.clear();
    HH.clear();
  }

  // Initialize neighbor lists
  vector<int> tmp; //an empty vector
  tmp.reserve(1000); 
  PdPd.assign(q_Pd0.size(), tmp); //Pd neighbors of Pd
  PdH.assign(q_Pd0.size(), tmp); //H neighbors of Pd
  HPd.assign(q_H0.size(), tmp); //Pd neighbors of H
  HH.assign(q_H0.size(), tmp); //H neighbors of H

  vector<int> tmp1; //a zero vector
  tmp1.assign(3, 0);

  // Store all the sites in a KD-Tree with K = 3.
  int Natoms = q_Pd0.size() + q_H0.size();
  PointIn3D *atoms = new PointIn3D[Natoms];
  for(int i=0; i<(int)q_Pd0.size(); i++)
    atoms[i] = PointIn3D(i, q_Pd0[i]);
  for(int i=0; i<(int)q_H0.size(); i++)
    atoms[q_Pd0.size() + i] = PointIn3D(q_Pd0.size() + i, q_H0[i]);
  KDTree<PointIn3D> atomTree(Natoms, atoms); //Note: atoms are re-ordered.

  // Find and store neighbors
  int maxNeib = 0;
  int maxNei = 1000; //max number of neighbors for each atom. 
  PointIn3D candidates[maxNei];
  double qLocal[3];
  for(int i=0; i<Natoms; i++) {
    int nNeib = 0; 
    for(int j=0; j<3; j++)
      qLocal[j] = atoms[i].val(j); // position of i (could be either Pd or H)
    // find neighbors of i within rc --> populates "candidates"
    int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
    // debug
    if(nFound>maxNei) {
      cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei 
           << ". Increase maxNei! " << endl;
      exit(-1);
    } 
/*
    for(int j=0; j<nFound; j++) {
      Vec3D thisone(qLocal[0],qLocal[1],qLocal[2]);
      Vec3D neib(candidates[j].val(0), candidates[j].val(1), candidates[j].val(2));
      Vec3D neib2;
      if (candidates[j].pid()<q_Pd0.size())
        neib2 = q_Pd0[candidates[j].pid()];
      else
        neib2 = q_H0[candidates[j].pid() - q_Pd0.size()];

      Vec3D dif = thisone - neib;
      Vec3D dif0 = neib - neib2;
      cout << "i = " << i << ", j = " << j << ", d = " << max(fabs(dif[0]), max(fabs(dif[1]), fabs(dif[2]))) << ", 0 = " << dif0.norm() << endl;
    }
*/

    // populates PdPd, PdH, HPd, HH  
    for(int j=0; j<nFound; j++) {
      int nei = candidates[j].pid(); //one neighbor
      //filter the same one and the outside ones -- KDTree works with "pseudo distance"
      if(atoms[i].pid()==nei) 
        continue;
      if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
         (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
         (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) > rc*rc)
        continue;
      nNeib++;
      if(atoms[i].pid()<(int)q_Pd0.size()) {//working on a Pd site
        if(nei<(int)q_Pd0.size()) {// this neighbor is Pd
          PdPd[atoms[i].pid()].push_back(nei);
        } else { // this neighbor is H
          PdH[atoms[i].pid()].push_back(nei-q_Pd0.size());
        }
      } else {// working on an H site
        if(nei<(int)q_Pd0.size()) {//this neighbor is Pd
          HPd[atoms[i].pid()-q_Pd0.size()].push_back(nei);
        } else {//this neighbor is H
          HH[atoms[i].pid()-q_Pd0.size()].push_back(nei-q_Pd0.size());
        }
      }
    }

    if(nNeib>maxNeib) maxNeib = nNeib;
  }
/*
  for(int i=0; i<q_Pd0.size(); i++) {
    cout << "Pd[" << i << "]:" << q_Pd0[i][0] << " " << q_Pd0[i][1] << " " << q_Pd0[i][2] << endl;
    for(int j=0; j<PdPd[i].size(); j++) {
      Vec3D delta = q_Pd0[i] - q_Pd0[PdPd[i][j]];
      double dist = delta.norm();
      cout << "  Nei " << j << ": id = " << PdPd[i][j] << ", " << q_Pd0[PdPd[i][j]][0] << " " <<  q_Pd0[PdPd[i][j]][1] << " " << q_Pd0[PdPd[i][j]][2] << ", d = " << dist << endl;
    }
  }
*/

  delete[] atoms;
  return maxNeib;
}

