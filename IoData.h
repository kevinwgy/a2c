#ifndef _IO_DATA_H_
#define _IO_DATA_H_
#include <cstdio>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>

/*********************************************************************
 * class IoData reads and processes the input data provided by the user
 *********************************************************************
*/
//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:
  map<int, DataType *> dataMap;

  void setup(const char *name, ClassAssigner *p) {
    SysMapObj<DataType> *smo = new SysMapObj<DataType>(name, &dataMap);
    if (p) p->addSmb(name, smo);
    else addSysSymbol(name, smo);
  }

  ~ObjectMap() {
    for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      delete it->second;
  }
};

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

struct SpaceData
{
  enum LatticeType {SC = 0, BCC = 1, FCC = 2, HCP = 3} lattice_type;

  enum SampleShape {CUBE = 0, SPHERE = 1, OCTAHEDRON = 2, OTHER = 3} sample_shape;

  int N; //!< sample size (NxNxN), can be overriden by "size"
  double size; //!< size (dimension: length) of the sample. Its specific meaning depends on sample shape

  double xe; //!< solute molar fraction

  double lattice_parameter; //!< default unit: Ang

  double sigma_0; //!< matrix atom freq. (Any^-2)
  double sigma_1; //1< solute atom frequency (Ang^-2)

  SpaceData();
  ~SpaceData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct PotentialData
{

  double eps; // tolerance
  double rc;  // cutoff distance (Ang)

  PotentialData();
  ~PotentialData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct DiffusionData
{
  double T0;  //!< temperature (K)

  //! absorption
  double v; //!< attempt frequency (s^-1)
  double Qm; //!< activation energy (eV)

  //! adsorption
  double vad; //!< attempt frequency (s^-1)
  double Qad; //!< activation energy (eV) 

  double xUb; //!< upper bound of hydrogen atomic fraction
  double xLb; //!< lower bound of hydrogen atomic fraction

  //! "boundary condition"
  double t_subsurf; //!< subsurface layer thickness
  double mubd; //!< chemical potential

  DiffusionData();
  ~DiffusionData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct TsData
{
  double dt;  //!< time-step (s)  TODO: if dt drops below 1e-8, we should change unit to reduce roundoff error
  double t_final;  //!< final time (unit: same as dt)

  TsData();
  ~TsData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct OutputData
{

  const char* foldername;
  const char* filename_base;

  int frequency;

  OutputData();
  ~OutputData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct RestartData
{
  //! continue setup
  int restart_file_num;

  //! number of sites in restart file
  int matrix_particle_number;
  int solute_particle_number;

  RestartData();
  ~RestartData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

class IoData 
{
  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  PotentialData potential;

  SpaceData space;

  DiffusionData diffusion;

  TsData ts;

  RestartData restart;

  OutputData output;
  
public:

  IoData() {}
  IoData(int, char**);
  ~Iodata() {}

  void readCmdLine(int, char**);
  void readCmdFile();
  void setupCmdFileVariables();
  
};

//------------------------------------------------------------------------------

#endif
