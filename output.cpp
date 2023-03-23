#include <time.h>
#include <output.h>
#include <input.h>
#include <Vector3D.h>
#include <version.h>
#include <cstring>
#include <unistd.h> //sleep
using namespace std;

//--------------------------------------------------------------------------------

Output::Output(MPI_Comm *comm_, Input* input, int argc, char* argv[]) : comm(comm_)
{

  printLogo(argc, argv);


  int MPI_rank = 0;
  MPI_Comm_rank(*comm, &MPI_rank);
  if(MPI_rank)
    return;
  
  sprintf(full_filename_base, "%s/%s", input->file.foldername, input->file.filename_base);
  sprintf(filename_base, "%s", input->file.filename_base);
  char f1[256];
  sprintf(f1, "%s_summary.txt", full_filename_base);

  summaryfile.open(f1, ios::out);
  summaryfile << "Revision: " << GIT_REV << " | " << "Branch: " << GIT_BRANCH << " | " << "Tag: " << GIT_TAG << endl;
  summaryfile << "Computation started at:" << endl;
  summaryfile << "  " << getCurrentDateTime() << endl;
  summaryfile << endl;
  summaryfile << "Inputs" << endl;
  summaryfile << "  T0 = " << input->file.T0 << endl;
  summaryfile << "  eps = " << input->file.eps << endl;
  summaryfile << "  rc = " << input->file.rc << endl;
  summaryfile << "  dt = " << input->file.dt << endl;
  summaryfile << "  t_final = " << input->file.t_final << endl;
  summaryfile << "  N = " << input->file.N << endl;
  summaryfile << "  xe = " << input->file.xe << endl;
  summaryfile << "  a_Pd = " << input->file.a_Pd << endl;
  summaryfile << "  sigma_Pd = " << input->file.sigma_Pd << endl;
  summaryfile << "  sigma_H = " << input->file.sigma_H << endl;
  summaryfile << "  output_frequency = " << input->file.output_frequency << endl;
  summaryfile << "  v = " << input->file.v << endl;
  summaryfile << "  Qm = " << input->file.Qm << endl;
  summaryfile << "  mubd = " << input->file.mubd << endl;
  summaryfile << "  t_subsurf = " << input->file.t_subsurf << endl;
}

//--------------------------------------------------------------------------------

Output::~Output()
{
  if(summaryfile.is_open()) summaryfile.close();
}

//--------------------------------------------------------------------------------

void Output::output_solution(int iFrame, int iTimeStep, vector<Vec3D> &q_Pd, vector<Vec3D> &q_H,
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

void Output::printLogo(int argc, char *argv[])
{
  int MPI_rank = 0, MPI_size = 0;
  MPI_Comm_size(*comm, &MPI_size);
  MPI_Comm_rank(*comm, &MPI_rank);

  if(!MPI_rank) {
    cout << endl;
    cout << "  ____                 _ _      _      _    ____   ____ " << endl;
    cout << " |  _ \\ __ _ _ __ __ _| | | ___| |    / \\  |___ \\ / ___|" << endl;
    cout << " | |_) / _` | '__/ _` | | |/ _ \\ |   / _ \\   __) | |    " << endl;
    cout << " |  __/ (_| | | | (_| | | |  __/ |  / ___ \\ / __/| |___ " << endl;
    cout << " |_|   \\__,_|_|  \\__,_|_|_|\\___|_| /_/   \\_\\_____|\\____|" << endl;
    cout << endl;
    cout << "Revision: " << GIT_REV << " | " << "Branch: " << GIT_BRANCH << " | " << "Tag: " << GIT_TAG << endl;
    cout << "Computation started at: " << getCurrentDateTime() << endl;
    cout << "Using " << MPI_size << " processor cores (including concurrent programs)." << endl;
    cout << "Command:";
    for(int i=0; i<argc; i++)
      cout << " " << argv[i];
    cout << endl;
    cout.flush();
  } else
    sleep(0.1);
}


