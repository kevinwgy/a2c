#include<SpaceOperator.h>
#include<LatticeVariables.h>
#include<Utils.h>
#include<GeoTools/DistancePointToParallelepiped.h>
#include<GeoTools/DistancePointToSpheroid.h>
using std::vector;

//---------------------------------------------------------------------

SpaceOperator::SpaceOperator(MPI_Comm& comm_, IoData& iod_, vector<MaterialOperator> &mato_,
                             vector<LatticeInfo> &linfo_) 
             : iod(iod_), comm(comm_), mato(mato_), linfo(linfo_)
{

}

//---------------------------------------------------------------------

SpaceOperator::~SpaceOperator()
{ }

//---------------------------------------------------------------------

void 
SpaceOperator::SetupLatticeVariables(vector<LatticeVariables> &LVS)
{
  CreateSpecimen(LVS); //!< initializes LVS.


}

//---------------------------------------------------------------------

void
SpaceOperator::CreateSpecimen(vector<LatticeVariables> &LVS)
{
  



}

//---------------------------------------------------------------------


//---------------------------------------------------------------------

initializeStateVariables(Input &input, vector<vector<Int3> > &SS1, vector<vector<Int3> > &SS2,
                              vector<Vec3D> &q_Pd0, vector<Vec3D> &q_H0, vector<double> &sigma_Pd0,
                              vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, 
                              vector<int> &PdSubsurf, vector<int> &HSubsurf, 
                              double &gammabd, vector<int> &full_H) 
{
    // initialize boundary conditions
    double &T = input.file.T0;
    const double kB = 8.6173324e-5; // (eV/K)
    //const double cutoff = 1.0e-20; // The minima of hydrogen atomic fraction
    gammabd = input.file.mubd/(kB*T);
    double thickness = input.file.t_subsurf; // (angstrom) thickness of subsurface layer

  if(!(q_Pd0.empty() && q_H0.empty() && sigma_Pd0.empty() && sigma_H0.empty() && x0.empty() && gamma.empty() && PdSubsurf.empty() && HSubsurf.empty())) {
    cerr << "WARNING: some vectors are not empty before initialization!" << endl;
    q_Pd0.clear();
    q_H0.clear();
    sigma_Pd0.clear();
    sigma_H0.clear();
    x0.clear();
    gamma.clear();
    PdSubsurf.clear();
    HSubsurf.clear();
  }

  // INITIALIZE q_Pd0
  q_Pd0.push_back(Vec3D(0.0,0.0,0.0)); //first atomic site

  int p = 1; // store index 
  // loop through shells to add sites
  for(int i=0; i<(int)SS1.size(); i++)
    for(int j=0; j<(int)SS1[i].size(); j++) {
      Int3 &site = SS1[i][j];
      int N = input.file.N;
      if(input.file.sample_shape == 1) { //cube
        if(site[0]>-N && site[0]<=N && site[1]>-N && site[1]<=N && site[2]>-N && site[2]<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_Pd0.push_back(q_site);
          if(site[0]<-(N-1-2.0*thickness/input.file.a_Pd) || site[0]>(N-2.0*thickness/input.file.a_Pd) || site[1]<-(N-1-2.0*thickness/input.file.a_Pd) || site[1]>(N-2.0*thickness/input.file.a_Pd) || site[2]<-(N-1-2.0*thickness/input.file.a_Pd) || site[2]>(N-2.0*thickness/input.file.a_Pd))
            PdSubsurf.push_back(p);
          p++;
        }
      }
      else if(input.file.sample_shape == 2) { //octahedron
        if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_Pd0.push_back(q_site);
          if(!(fabs(site[0])+fabs(site[1])+fabs(site[2])<=(N-2.0*thickness/input.file.a_Pd*1.7321)))
            PdSubsurf.push_back(p);
          p++;
        }   
      }
      else if(input.file.sample_shape == 3) { // sphere
        if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_Pd0.push_back(q_site);
          if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] >= (N-2.0*thickness/input.file.a_Pd)*(N-2.0*thickness/input.file.a_Pd))
            PdSubsurf.push_back(p);
          p++;
        }
      }
    }
  
  // INITIALIZE q_H0
  p = 0; // store index 
  for(int i=0; i<(int)SS2.size(); i++)
    for(int j=0; j<(int)SS2[i].size(); j++) {
      Int3 &site = SS2[i][j];
      int N = input.file.N;
      if(input.file.sample_shape == 1) { //cube
        if(site[0]>-N && site[0]<=N && site[1]>-N && site[1]<=N && site[2]>-N && site[2]<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_H0.push_back(q_site);
          if(site[0]<-(N-1-2.0*thickness/input.file.a_Pd) || site[0]>(N-2.0*thickness/input.file.a_Pd) || site[1]<-(N-1-2.0*thickness/input.file.a_Pd) || site[1]>(N-2.0*thickness/input.file.a_Pd) || site[2]<-(N-1-2.0*thickness/input.file.a_Pd) || site[2]>(N-2.0*thickness/input.file.a_Pd))
            HSubsurf.push_back(p);
          p++;
        }
      }
      else if(input.file.sample_shape == 2) { //octahedron
        if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_H0.push_back(q_site);
          if(!(fabs(site[0])+fabs(site[1])+fabs(site[2])<=(N-2.0*thickness/input.file.a_Pd*1.7321)))
            HSubsurf.push_back(p);
          p++;
        }
      }
      else if(input.file.sample_shape == 3) { //sphere
        if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_Pd;
          q_H0.push_back(q_site);
          if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] >= (N-2.0*thickness/input.file.a_Pd)*(N-2.0*thickness/input.file.a_Pd))
            HSubsurf.push_back(p);
          p++;
        }
      }
    }

  int nPd = q_Pd0.size();
  int nH  = q_H0.size();

  // Initialize sigma, x, gamma
  sigma_Pd0.assign(nPd, input.file.sigma_Pd);
  sigma_H0.assign(nH, input.file.sigma_H);
  x0.assign(nH, input.file.xe);
  gamma.assign(nH, 0.0);
  full_H.assign(nH, 0);

/*
  // debug
  for(int i=0; i<nPd; i++) 
    sigma_Pd0[i] = (double)i;
  for(int i=0; i<nH; i++) 
    sigma_H0[i] = (double)(i+nPd);
*/

  return;
}
