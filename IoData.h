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

  int order; //!< set operation order (0, 1, ...)

  StateVariable initialConditions;

  PlaneData();
  ~PlaneData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SphereData {

  double cen_x, cen_y, cen_z, radius;

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  int order; //!< set operation order (0, 1, ...)

  StateVariable initialConditions;

  SphereData();
  ~SphereData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct ParallelepipedData {

  double x0, y0, z0; //!< point 1

  double ax, ay, az; //!< axis 1 and its length
  double bx, by, bz; //!< axis 2 and its length
  double cx, cy, cz; //!< axis 3 and its length

  int order; //!< set operation order (0, 1, ...)

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  StateVariable initialConditions;

  ParallelepipedData();
  ~ParallelepipedData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SpheroidData {

  double cen_x, cen_y, cen_z;
  double axis_x, axis_y, axis_z;
  double semi_length; //!< half length (along axis)
  double radius; //!< max. radius on the transverse plane

  enum InteriorOrExterior {INTERIOR = 0, EXTERIOR = 1} side;
  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  int order; //!< set operation order (0, 1, ...)

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

  int order; //!< set operation order (0, 1, ...)

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

  int order; //!< set operation order (0, 1, ...)

  StateVariable initialConditions;

  CylinderSphereData();
  ~CylinderSphereData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct CustomGeometryData {

  const char *specifier;

  enum Inclusion {OVERRIDE = 0, INTERSECTION = 1, UNION = 2} inclusion;

  int order; //!< set operation order (0, 1, ...)

  StateVariable initialConditions;

  CustomGeometryData();
  ~CustomGeometryData() {}

  void setup(const char *, ClassAssigner *);

};

//------------------------------------------------------------------------------
//! Glocal species ID & name
struct GlobalSpeciesData 
{
  const char* name;

  GlobalSpeciesData();
  ~GlobalSpeciesData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------
//! For specifying material composition & initial conditions
struct LocalSpeciesData
{
  int speciesID; //the global species ID

  enum Boolean {FALSE = 0, TRUE = 1} diffusive; 

  double molar_fraction; //!< initial value (optional)
  
  double min_molar_fraction; //!< a smaller number to avoid division-by-zero

  LocalSpeciesData();
  ~LocalSpeciesData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

struct LatticeSiteData
{
  double la, lb, lc;
  int materialid;

  //! For specifying the initial conditions for certain species
  ObjectMap<LocalSpeciesData> speciesMap; 
    
  //! Initial site condition
  double atomic_frequency; //!< initial value (optional)
  double mean_displacement_x; //!< initial value (optional)
  double mean_displacement_y;
  double mean_displacement_z;
  double mean_momentum_x; //!< initial value (optional)
  double mean_momentum_y;
  double mean_momentum_z;

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
  
  CustomGeometryData custom_geometry;
  
  double ax, ay, az, bx, by, bz, cx, cy, cz; //!< lattice vectors
  double ox, oy, oz; //!< origin of lattice, i.e. the (x,y,z) coords of site (0,0,0)
   
  enum Boolean {FALSE = 0, TRUE = 1} a_periodic, b_periodic, c_periodic; //!< whether a,b,c dirs are periodic

  ObjectMap<LatticeSiteData> siteMap;
  
  double dmin; //!< min spacing between sites on this lattice and any other lattice 

  LatticeData();
  ~LatticeData() {}

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
//! Each "material" has an interatomic potential
struct MaterialData 
{
  const char* name;

  //! The user specifieds the global species ID of the species involved in each material.
  //! Starting from species0, following the ascending order. 
  //! The function IoData::finalize is called immediately after reading the input file.
  //! It sets the array species[MAX_SPECIES] and nSpecies. Unused species have ID = -1
  const static int MAX_SPECIES = 12;
  int species0, species1, species2, species3, species4, species5, species6;
  int species7, species8, species9, species10, species11;
  int species[MAX_SPECIES];
  int nSpecies;

  InteratomicPotentialData potential;

  MaterialData();
  ~MaterialData() {}

  Assigner *getAssigner();
};

//------------------------------------------------------------------------------

//! specify initial conditions for material(s) and species(s) within a region
struct RegionalIcData
{
  //! Variables & structures for determining the region and the lattice sites where IC is applied
  ObjectMap<PlaneData>          planeMap;
  ObjectMap<CylinderConeData>   cylinderconeMap;
  ObjectMap<CylinderSphereData> cylindersphereMap;
  ObjectMap<SphereData>         sphereMap;
  ObjectMap<ParallelepipedData> parallelepipedMap;
  ObjectMap<SpheroidData>       spheroidMap;
  CustomGeometryData custom_geometry;
  int latticeID;
  int materialID;
  int siteID;
  
  //! For specifying the initial conditions for certain species
  ObjectMap<LocalSpeciesData> speciesMap; 
    
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

  //! "boundary condition"
  double t_subsurf; //!< subsurface layer thickness
  double mubd; //!< chemical potential

  DiffusionData();
  ~DiffusionData() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct ExplicitData {

  //!time-integration scheme used
  enum Type {FORWARD_EULER = 0, RUNGE_KUTTA_2 = 1, RUNGE_KUTTA_3 = 2} type;

  ExplicitData();
  ~ExplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TsData
{
  enum Type {EXPLICIT = 0, IMPLICIT = 1} type;
  int maxIts;
  double timestep;
  double cfl;
  double maxTime;

  ExplicitData expl;

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

  ObjectMap<GlobalSpeciesData> speciesMap; //!< global names & ids of species

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
  
  void finalize(); //!< should be called right after creating IoData & communicator
};

//------------------------------------------------------------------------------

#endif
