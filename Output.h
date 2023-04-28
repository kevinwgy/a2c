#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <IoData.h>
#include <vector>

struct Vec3D;
class Variables;

class Output {

  MPI_Comm &comm;
  IoData &iod;

  char filename_base[256];

  int iFrame;
  double last_snapshot_time;

public:

  Output(MPI_Comm &comm_, IoData &iod_);
  ~Output();
  
  void OutputSolution(double time, double dt, int time_step, Variables &V, bool force_write = false);

private:

  void OutputSolutionVTP(double time, int time_step, Variables &V);
  void OutputSolutionXYZ(double time, int time_step, Variables &V);

};
#endif
