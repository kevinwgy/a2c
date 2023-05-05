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

struct StateVariable {
  //not used at the moment
  StateVariable() {}
  ~StateVariable() {}
  void setup(const char *, ClassAssigner * = 0) {}
};

//------------------------------------------------------------------------------

struct PlaneData {

  double cen_x, cen_y, cen_z, nx, ny, nz;

  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  PlaneData();
  ~PlaneData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct ParallelepipedData {

  double x0, y0, z0; //!< point 1

  double ax, ay, az; //!< axis 1 and its length
  double bx, by, bz; //!< axis 2 and its length
  double cx, cy, cz; //!< axis 3 and its length

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  ParallelepipedData();
  ~ParallelepipedData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SphereData {

  double cen_x, cen_y, cen_z, radius;

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  SphereData();
  ~SphereData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SpheroidData {

  double cen_x, cen_y, cen_z;
  double axis_x, axis_y, axis_z;
  double length, diameter;

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  SpheroidData();
  ~SpheroidData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct CylinderConeData {

  //! info about the cylinder
  double cen_x, cen_y, cen_z, nx, ny, nz, r, L;
  //! info about the cone (connected to the top of the cylinder)
  double opening_angle_degrees, cone_height;

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  CylinderConeData();
  ~CylinderConeData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct CylinderSphereData {

  double cen_x, cen_y, cen_z, nx, ny, nz, r, L;

  enum OnOff {Off = 0, On = 1};
  OnOff front_cap;
  OnOff back_cap;

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  CylinderSphereData();
  ~CylinderSphereData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct LatticeSiteData
{
  double la, lb, lc;
  int materialid;

  LatticeSiteData();
  ~LatticeSiteData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

struct LatticeData
{

  //! Lattice domain
  ObjectMap<PlaneData>          planeMap;
  ObjectMap<CylinderConeData>   cylinderconeMap;
  ObjectMap<CylinderSphereData> cylindersphereMap;
  ObjectMap<SphereData>         sphereMap;
  ObjectMap<ParallelepipedData> parallelepipedMap;
  ObjectMap<SpheroidData>       spheroidMap;
  const char *custom_geometry_specifier;

  double ax, ay, az, bx, by, bz, cx, cy, cz; //!< lattice vectors
  double ox, oy, oz; //!< origin of lattice, i.e. the (x,y,z) coords of site (0,0,0)
   
  enum Boolean {FALSE = 0, TRUE = 1} a_periodic, b_periodic, c_periodic; //!< whether a, b, c dirs have periodic b.c.

  ObjectMap<LatticeSiteData> siteMap;
  
  double dmin; //!< min spacing between sites on this lattice and those on other lattices

  LatticeData();
  ~LatticeData() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SpeciesData
{
  int id; //!< species id
  const char* name; //!< species name (optional)

  enum Boolean {FALSE = 0, TRUE = 1} diffusive; 

  double molar_fraction; //!< initial value (optional)
  
  double min_molar_fraction; //!< a smaller number to avoid division-by-zero

  SpeciesData();
  ~SpeciesData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

struct InteratomicPotentialData
{
  const char* filename;

  double eps; // tolerance
  double rc;  // cutoff distance (Ang)

  InteratomicPotentialData();
  ~InteratomicPotentialData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct MaterialData 
{

  const char* name;

  ObjectMap<SpeciesData> speciesMap;
    
  InteratomicPotentialData potential;

  MaterialData();
  ~MaterialData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

//! pecify initial conditions for material(s) and species(s) within a region
struct RegionalIcData
{
  ObjectMap<PlaneData>          planeMap;
  ObjectMap<CylinderConeData>   cylinderconeMap;
  ObjectMap<CylinderSphereData> cylindersphereMap;
  ObjectMap<SphereData>         sphereMap;
  ObjectMap<ParallelepipedData> parallelepipedMap;
  ObjectMap<SpheroidData>       spheroidMap;
  const char *custom_geometry_specifier;

  ObjectMap<SpeciesData> speciesMap;
    
  double atomic_frequency; //!< initial value (optional)
  double mean_displacement_x; //!< initial value (optional)
  double mean_displacement_y;
  double mean_displacement_z;
  double mean_momentum_x; //!< initial value (optional)
  double mean_momentum_y;
  double mean_momentum_z;

  RegionalIcData();
  ~RegionalIcData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

struct DiffusionData
{
  double T0;  //!< temperature (K)

  //! absorption
  double v; //!< attempt frequency (s^-1)
  double Qm; //!< activation energy (eV)

  //! adsorption
//  double vad; //!< attempt frequency (s^-1)
//  double Qad; //!< activation energy (eV) 

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
  enum FileFormat {VTP = 0, XYZ = 1, VTP_and_XYZ = 2} format;

  const char* prefix;
  const char* solution_filename_base;

  enum VerbosityLevel {LOW = 0, MEDIUM = 1, HIGH = 2} verbose;

  int frequency;
  double frequency_dt;

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
  int matrix_site_number;
  int solute_site_number;

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

  ObjectMap<LatticeData> latticeMap;

  ObjectMap<MaterialData> materialMap; //!< a ``material'' can/often has multiple species

  ObjectMap<RegionalIcData> icMap;

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
