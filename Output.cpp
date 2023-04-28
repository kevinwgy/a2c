#include <Output.h>
#include <Vector3D.h>
#include <Utils.h>
#include <cstring>
using std::vector;

//--------------------------------------------------------------------------------

Output::Output(MPI_Comm &comm_, IoData& iod_) : comm(comm_), iod(iod_)
{
  sprintf(filename_base, "%s%s", iod.output.prefix, iod.output.solution_filename_base);

  iFrame = 0;

  last_snapshot_time = -1.0;
}

//--------------------------------------------------------------------------------

Output::~Output()
{
}

//--------------------------------------------------------------------------------

void Output::OutputSolution(double time, double dt, int time_step, Variables &V, bool force_write)
{

  if(isTimeToWrite(time, dt, time_step, iod.output.frequency_dt, iod.output.frequency,
                   last_snapshot_time, force_write)) {
    if(iod.output.format == Output::VTP || iod.output.format == Output::VTP_and_XYZ)
      OutputSolutionVTP(time, time_step, V);
    else if(iod.output.format == Output::XYZ || iod.output.format == Output::VTP_and_XYZ)
      OutputSolutionXYZ(time, time_step, V);
  }

}

//--------------------------------------------------------------------------------

void Output::OutputSolutionXYZ(double time, int time_step, Variables &V) 

int iFrame, int iTimeStep, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H,
                     vector<double> &sigma_Pd, vector<double> &sigma_H, vector<double> &x, 
                     vector<double> &gamma, vector<int> &full_H)
{
  int MPI_rank = 0;
  MPI_Comm_rank(*comm, &MPI_rank);
  if(MPI_rank)
    return;

  // find bounds of simulation box
  double xmax, xmin, ymax, ymin, zmax, zmin;
  xmax = q_Pd[0][0]; xmin = q_Pd[0][0];
  ymax = q_Pd[0][1]; ymin = q_Pd[0][1];
  zmax = q_Pd[0][2]; zmin = q_Pd[0][2];

  for(int i=0; i<(int)q_Pd.size(); i++) {
    if(q_Pd[i][0]>xmax) xmax = q_Pd[i][0];
    if(q_Pd[i][0]<xmin) xmin = q_Pd[i][0];
    if(q_Pd[i][1]>ymax) ymax = q_Pd[i][1];
    if(q_Pd[i][1]<ymin) ymin = q_Pd[i][1];
    if(q_Pd[i][2]>zmax) zmax = q_Pd[i][2];
    if(q_Pd[i][2]<zmin) zmin = q_Pd[i][2];
  }

  // create solution file
  char full_fname[256];

  sprintf(full_fname, "%s.%d.dat", full_filename_base, iFrame);
  ofstream file(full_fname, ios::out);

  // write time step
  file << "ITEM: TIMESTEP" << endl;
  file << iTimeStep  << endl;

  // write number of atomic sites
  file << "ITEM: NUMBER OF ATOMS" << endl;
  file <<  q_Pd.size()+q_H.size() << endl;

  // write bounds of simulation box
  file << "ITEM: BOX BOUNDS ss ss ss" << endl;
  file << scientific << xmin << " " << scientific << xmax  << endl;
  file << scientific << ymin << " " << scientific << ymax  << endl;
  file << scientific << zmin << " " << scientific << zmax  << endl;

  // write properties on each site
  file << "ITEM: ATOMS id type x y z sigma H_fraction chemical_potential full_occupancy" << endl;

  // Pd sites
  for(int i=0; i<(int)q_Pd.size(); i++)
    file << i+1 << " " << 1 << " " << scientific << q_Pd[i][0] << " " << scientific << q_Pd[i][1] << " " << scientific << q_Pd[i][2] << " "
         << scientific << sigma_Pd[i] << " " << scientific << 1.0 << " " << scientific << 0.0 << " " << 1 << endl;
  // H sites
  for(int i=0; i<(int)q_H.size(); i++)
    file << i+1+q_Pd.size() << " " << 2 << " " << scientific << q_H[i][0] << " " << scientific << q_H[i][1] << " " << scientific << q_H[i][2] << " " 
         << scientific << sigma_H[i] << " " << scientific << x[i] << " " << scientific << gamma[i] << " " << full_H[i] << endl;

  return;
}

//--------------------------------------------------------------------------------

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string Output::getCurrentDateTime() 
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    if(strlen(tzname[1]) != 0)
      strcat(buf, tzname[1]); //daylight saving time
    else
      strcat(buf, tzname[0]); //standard time

    return buf;
}

//--------------------------------------------------------------------------------

