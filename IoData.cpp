#include <iostream>
#include <IoData.h>

//------------------------------------------------------------------------------

PlaneData::PlaneData()
{

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 0.0;

  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *PlaneData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 9, nullAssigner);

  new ClassDouble<PlaneData> (ca, "Point_x", this, &PlaneData::cen_x);
  new ClassDouble<PlaneData> (ca, "Point_y", this, &PlaneData::cen_y);
  new ClassDouble<PlaneData> (ca, "Point_z", this, &PlaneData::cen_z);
  new ClassDouble<PlaneData> (ca, "Normal_x", this, &PlaneData::nx);
  new ClassDouble<PlaneData> (ca, "Normal_y", this, &PlaneData::ny);
  new ClassDouble<PlaneData> (ca, "Normal_z", this, &PlaneData::nz);

  initialConditions.setup("InitialState", ca);

  new ClassToken<PlaneData> (ca, "Inclusion", this,
     reinterpret_cast<int PlaneData::*>(&PlaneData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<PlaneData>(ca, "OperationOrder", this, &PlaneData::order);

  return ca;
}

//------------------------------------------------------------------------------

SphereData::SphereData()
{
  
  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  radius = -1.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *SphereData::getAssigner()
{
  
  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);
  
  new ClassDouble<SphereData> (ca, "Center_x", this, &SphereData::cen_x);
  new ClassDouble<SphereData> (ca, "Center_y", this, &SphereData::cen_y);
  new ClassDouble<SphereData> (ca, "Center_z", this, &SphereData::cen_z);
  new ClassDouble<SphereData> (ca, "Radius", this, &SphereData::radius);
  
  new ClassToken<SphereData> (ca, "Side", this,
     reinterpret_cast<int SphereData::*>(&SphereData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<SphereData> (ca, "Inclusion", this,
     reinterpret_cast<int SphereData::*>(&SphereData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<SphereData>(ca, "OperationOrder", this, &SphereData::order);

  initialConditions.setup("InitialState", ca);
  
  return ca;
}

//------------------------------------------------------------------------------

ParallelepipedData::ParallelepipedData()
{

  x0 = y0 = z0 = 0.0;
  ax = ay = az = 0.0;
  bx = by = bz = 0.0;
  cx = cy = cz = 0.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *ParallelepipedData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 16, nullAssigner);

  new ClassDouble<ParallelepipedData> (ca, "X0", this, &ParallelepipedData::x0);
  new ClassDouble<ParallelepipedData> (ca, "Y0", this, &ParallelepipedData::y0);
  new ClassDouble<ParallelepipedData> (ca, "Z0", this, &ParallelepipedData::z0);

  new ClassDouble<ParallelepipedData> (ca, "Ax", this, &ParallelepipedData::ax);
  new ClassDouble<ParallelepipedData> (ca, "Ay", this, &ParallelepipedData::ay);
  new ClassDouble<ParallelepipedData> (ca, "Az", this, &ParallelepipedData::az);

  new ClassDouble<ParallelepipedData> (ca, "Bx", this, &ParallelepipedData::bx);
  new ClassDouble<ParallelepipedData> (ca, "By", this, &ParallelepipedData::by);
  new ClassDouble<ParallelepipedData> (ca, "Bz", this, &ParallelepipedData::bz);

  new ClassDouble<ParallelepipedData> (ca, "Cx", this, &ParallelepipedData::cx);
  new ClassDouble<ParallelepipedData> (ca, "Cy", this, &ParallelepipedData::cy);
  new ClassDouble<ParallelepipedData> (ca, "Cz", this, &ParallelepipedData::cz);

  new ClassToken<ParallelepipedData> (ca, "Side", this,
     reinterpret_cast<int ParallelepipedData::*>(&ParallelepipedData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<ParallelepipedData> (ca, "Inclusion", this,
     reinterpret_cast<int ParallelepipedData::*>(&ParallelepipedData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<ParallelepipedData>(ca, "OperationOrder", this, &ParallelepipedData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

SpheroidData::SpheroidData()
{

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;

  axis_x = 0.0;
  axis_y = 0.0;
  axis_z = 0.0;

  semi_length = 0.0;
  radius = 0.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *SpheroidData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 12, nullAssigner);

  new ClassDouble<SpheroidData> (ca, "Center_x", this, &SpheroidData::cen_x);
  new ClassDouble<SpheroidData> (ca, "Center_y", this, &SpheroidData::cen_y);
  new ClassDouble<SpheroidData> (ca, "Center_z", this, &SpheroidData::cen_z);
  new ClassDouble<SpheroidData> (ca, "Axis_x", this, &SpheroidData::axis_x);
  new ClassDouble<SpheroidData> (ca, "Axis_y", this, &SpheroidData::axis_y);
  new ClassDouble<SpheroidData> (ca, "Axis_z", this, &SpheroidData::axis_z);
  new ClassDouble<SpheroidData> (ca, "SemiLength", this, &SpheroidData::semi_length);
  new ClassDouble<SpheroidData> (ca, "Radius", this, &SpheroidData::radius);

  new ClassToken<SpheroidData> (ca, "Side", this,
     reinterpret_cast<int SpheroidData::*>(&SpheroidData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<SpheroidData> (ca, "Inclusion", this,
     reinterpret_cast<int SpheroidData::*>(&SpheroidData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<SpheroidData>(ca, "OperationOrder", this, &SpheroidData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

CylinderConeData::CylinderConeData() {

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 1.0;

  r = 1.0;
  L = 0.0;

  cone_height = 0.0;
  opening_angle_degrees = 45.0;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *CylinderConeData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 14, nullAssigner);

  new ClassDouble<CylinderConeData> (ca, "Axis_x", this, &CylinderConeData::nx);
  new ClassDouble<CylinderConeData> (ca, "Axis_y", this, &CylinderConeData::ny);
  new ClassDouble<CylinderConeData> (ca, "Axis_z", this, &CylinderConeData::nz);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_x", this, &CylinderConeData::cen_x);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_y", this, &CylinderConeData::cen_y);
  new ClassDouble<CylinderConeData> (ca, "BaseCenter_z", this, &CylinderConeData::cen_z);
  new ClassDouble<CylinderConeData> (ca, "CylinderRadius", this, &CylinderConeData::r);
  new ClassDouble<CylinderConeData> (ca, "CylinderHeight", this, &CylinderConeData::L);

  new ClassDouble<CylinderConeData> (ca, "ConeOpeningAngleInDegrees", this, &CylinderConeData::opening_angle_degrees);
  new ClassDouble<CylinderConeData> (ca, "ConeHeight", this, &CylinderConeData::cone_height);

  new ClassToken<CylinderConeData> (ca, "Side", this,
     reinterpret_cast<int CylinderConeData::*>(&CylinderConeData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<CylinderConeData> (ca, "Inclusion", this,
     reinterpret_cast<int CylinderConeData::*>(&CylinderConeData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<CylinderConeData>(ca, "OperationOrder", this, &CylinderConeData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//------------------------------------------------------------------------------

CylinderSphereData::CylinderSphereData() {

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 1.0;

  r = 1.0;
  L = 0.0;

  front_cap = Off;
  back_cap = Off;

  side = INTERIOR;
  inclusion = OVERRIDE;

  order = 0;
}

//------------------------------------------------------------------------------

Assigner *CylinderSphereData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 14, nullAssigner);

  new ClassDouble<CylinderSphereData> (ca, "Axis_x", this, &CylinderSphereData::nx);
  new ClassDouble<CylinderSphereData> (ca, "Axis_y", this, &CylinderSphereData::ny);
  new ClassDouble<CylinderSphereData> (ca, "Axis_z", this, &CylinderSphereData::nz);
  new ClassDouble<CylinderSphereData> (ca, "BaseCenter_x", this, &CylinderSphereData::cen_x);
  new ClassDouble<CylinderSphereData> (ca, "BaseCenter_y", this, &CylinderSphereData::cen_y);
  new ClassDouble<CylinderSphereData> (ca, "BaseCenter_z", this, &CylinderSphereData::cen_z);
  new ClassDouble<CylinderSphereData> (ca, "CylinderRadius", this, &CylinderSphereData::r);
  new ClassDouble<CylinderSphereData> (ca, "CylinderHeight", this, &CylinderSphereData::L);
  new ClassToken<CylinderSphereData> (ca, "FrontSphericalCap", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::front_cap), 2,
     "Off", 0, "On", 1);
  new ClassToken<CylinderSphereData> (ca, "BackSphericalCap", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::back_cap), 2,
     "Off", 0, "On", 1);

  new ClassToken<CylinderSphereData> (ca, "Side", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<CylinderSphereData> (ca, "Inclusion", this,
     reinterpret_cast<int CylinderSphereData::*>(&CylinderSphereData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<CylinderSphereData>(ca, "OperationOrder", this, &CylinderSphereData::order);

  initialConditions.setup("InitialState", ca);

  return ca;
}

//---------------------------------------------------------

CustomGeometryData::CustomGeometryData()
{
  specifier = "";

  inclusion = OVERRIDE;

  order = 0;
}

//---------------------------------------------------------

void CustomGeometryData::setup(const char *name)
{
  ClassAssigner *ca = new ClassAssigner(name, 4, 0);

  new ClassStr<CustomGeometryData> (ca, "GeometrySpecifier", this, &CustomGeometryData::specifier);

  new ClassToken<CustomGeometryData> (ca, "Inclusion", this,
     reinterpret_cast<int CustomGeometryData::*>(&CustomGeometryData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<CustomGeometryData>(ca, "OperationOrder", this, &CustomGeometryData::order);

  initialConditions.setup("InitialState", ca);
}

//---------------------------------------------------------

LatticeSiteData::LatticeSiteData()
{
  la = lb = lc = -1.0;
  materialid = -1;
}

//---------------------------------------------------------

Assigner *LatticeSiteData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);

  new ClassDouble<LatticeSiteData> (ca, "La", this, &LatticeSiteData::la);
  new ClassDouble<LatticeSiteData> (ca, "Lb", this, &LatticeSiteData::lb);
  new ClassDouble<LatticeSiteData> (ca, "Lc", this, &LatticeSiteData::lc);
  new ClassInt<LatticeSiteData> (ca, "MaterialID", this, &LatticeSiteData::materialid);
}

//---------------------------------------------------------

LatticeData::LatticeData()
{
  ax = 1.0; ay = 0.0; az = 0.0;
  bx = 0.0; by = 1.0; bz = 0.0;
  cx = 0.0; cy = 0.0; cz = 1.0;
  ox = oy = oz = 0.0;

  a_periodic = FALSE;
  b_periodic = FALSE;
  c_periodic = FALSE;

  dmin = 0.0; 
}

//---------------------------------------------------------

Assigner *LatticeData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner(name, 24, 0);

  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  parallelepiped.setup("Parallelepiped", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);

  custom_geometry.setup("CustomGeometry", ca);

  new ClassDouble<LatticeData> (ca, "OAx", this, &LatticeData::ax);
  new ClassDouble<LatticeData> (ca, "OAy", this, &LatticeData::ay);
  new ClassDouble<LatticeData> (ca, "OAz", this, &LatticeData::az);
  new ClassDouble<LatticeData> (ca, "OBx", this, &LatticeData::bx);
  new ClassDouble<LatticeData> (ca, "OBy", this, &LatticeData::by);
  new ClassDouble<LatticeData> (ca, "OBz", this, &LatticeData::bz);
  new ClassDouble<LatticeData> (ca, "OCx", this, &LatticeData::cx);
  new ClassDouble<LatticeData> (ca, "OCy", this, &LatticeData::cy);
  new ClassDouble<LatticeData> (ca, "OCz", this, &LatticeData::cz);

  new ClassDouble<LatticeData> (ca, "Ox", this, &LatticeData::ox);
  new ClassDouble<LatticeData> (ca, "Oy", this, &LatticeData::oy);
  new ClassDouble<LatticeData> (ca, "Oz", this, &LatticeData::oz);
  
  new ClassToken<LatticeData>(ca, "PeriodicOA", this,
                              reinterpret_cast<int LatticeData::*>(&LatticeData::a_periodic), 2,
                              "False", 0, "True", 1);
  new ClassToken<LatticeData>(ca, "PeriodicOB", this,
                              reinterpret_cast<int LatticeData::*>(&LatticeData::b_periodic), 2,
                              "False", 0, "True", 1);
  new ClassToken<LatticeData>(ca, "PeriodicOC", this,
                              reinterpret_cast<int LatticeData::*>(&LatticeData::c_periodic), 2,
                              "False", 0, "True", 1);

  siteMap.setup("Site", ca);

  new ClassDouble<LatticeData> (ca, "MinDistance", this, &LatticeData::dmin);
  
}

//---------------------------------------------------------

SpeciesData::SpeciesData()
{
  id = -1;
  name = "";
  diffusive = TRUE;
  molar_fraction = 1.0e-3;
  min_molar_fraction = 1.0e-4;
}

//---------------------------------------------------------

Assigner *SpeciesData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 5, nullAssigner);

  new ClassInt<SpeciesData> (ca, "ID", this, &SpeciesData::id);
  new ClassStr<SpeciesData> (ca, "SpeciesName", this, &SpeciesData::name);

  new ClassToken<SpeciesData>(ca, "Diffusive", this,
                               reinterpret_cast<int SpeciesData::*>(&SpeciesData::diffusion), 2,
                               "False", 0, "True", 1);

  new ClassDouble<SpeciesData> (ca, "MolarFraction", this, &SpeciesData::molar_fraction);

  new ClassDouble<SpeciesData> (ca, "MinimumMolarFraction", this, &SpeciesData::min_molar_fraction);

}

//---------------------------------------------------------

InteratomicPotentialData::InteratomicPotentialData() 
{
  filename = "";
  eps = 1.0e-5; 
  rc = 5.35; 
}

//---------------------------------------------------------

void InteratomicPotentialData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 3, 0);

  new ClassStr<InteratomicPotentialData> (ca, "File", this, &InteratomicPotentialData::filename);

  new ClassDouble<InteratomicPotentialData> (ca, "Tolerance", this, &InteratomicPotentialData::eps);

  new ClassDouble<InteratomicPotentialData> (ca, "CutOffDistance", this, &InteratomicPotentialData::rc);

}

//---------------------------------------------------------

MaterialData::MaterialData()
{
  name = "";
}

//---------------------------------------------------------

Assigner *MaterialData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner(name, 3, 0);

  new ClassStr<MaterialData> (ca, "MaterialName", this, &MaterialData::name);

  speciesMap.setup("Species", ca);

  potential.setup("InteratomicPotential");
}

//---------------------------------------------------------

RegionalIcData::RegionalIcData()
{
  atomic_frequency = 0.0; 
  mean_displacement_x = mean_displacement_y = mean_displacement_z = 0.0;
  mean_momentum_x     = mean_momentum_y     = mean_momentum_z     = 0.0;
}

//---------------------------------------------------------

Assigner *RegionalIcData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner(name, 15, 0);

  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  parallelepiped.setup("Parallelepiped", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);

  custom_geometry.setup("CustomGeometry", ca);

  new ClassDouble<SpeciesData> (ca, "AtomicFrequency", this, &SpeciesData::atomic_frequency);
  new ClassDouble<SpeciesData> (ca, "MeanDisplacementX", this, &SpeciesData::mean_displacement_x);
  new ClassDouble<SpeciesData> (ca, "MeanDisplacementY", this, &SpeciesData::mean_displacement_y);
  new ClassDouble<SpeciesData> (ca, "MeanDisplacementZ", this, &SpeciesData::mean_displacement_z);
  new ClassDouble<SpeciesData> (ca, "MeanMomentumX", this, &SpeciesData::mean_momentum_x);
  new ClassDouble<SpeciesData> (ca, "MeanMomentumY", this, &SpeciesData::mean_momentum_y);
  new ClassDouble<SpeciesData> (ca, "MeanMomentumZ", this, &SpeciesData::mean_momentum_z);

  speciesMap.setup("Species", ca);

}

//---------------------------------------------------------

DiffusionData::DiffusionData()
{
  T0 = 300; 
  v = 0.0;
  Qm = 0.0;
  t_subsurf = 0.0;
  mubd = -2.0; // dimensional (pressure)
}

//---------------------------------------------------------

void DiffusionData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 5, 0);

  new ClassDouble<DiffusionData> (ca, "Temperature", this, &DiffusionData::T0);

  new ClassDouble<DiffusionData> (ca, "AttemptFrequency", this, &DiffusionData::v);
  new ClassDouble<DiffusionData> (ca, "ActivationEnergy", this, &DiffusionData::Qm);

  new ClassDouble<DiffusionData> (ca, "SubsurfaceThickness", this, &DiffusionData::t_subsurf);
  new ClassDouble<DiffusionData> (ca, "SoluteChemicalPotential", this, &DiffusionData::mubd);

}

//---------------------------------------------------------

TsData::TsData()
{
  dt = 1.0e-6; 
  t_final = 0.01;
}

//---------------------------------------------------------

void TsData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 2, 0);

  new ClassDouble<TsData> (ca, "TimeStep", this, &TsData::dt);

  new ClassDouble<TsData> (ca, "MaxTime", this, &TsData::t_final);

}

//---------------------------------------------------------

OutputData::OutputData()
{
  format = VTP;
  prefix = "";
  solution_filename_base = "solution";
  frequency = 0;
  frequency_dt = -1.0;
  verbose = LOW;
}

//---------------------------------------------------------

void OutputData::setup(const char *name)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, 0);

  new ClassToken<OutputData>(ca, "Format", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::format), 4,
                               "VTP", 0, "XYZ", 1, "VTP_and_XYZ", 2, "XYZ_and_VTP", 2);

  new ClassStr<OutputData> (ca, "Prefix", this, &OutputData::prefix);

  new ClassStr<OutputData> (ca, "Solution", this, &OutputData::solution_filename_base);

  new ClassInt<OutputData> (ca, "Frequency", this, &OutputData::frequency);

  new ClassDouble<OutputData>(ca, "TimeInterval", this, &OutputData::frequency_dt);

  new ClassToken<OutputData>(ca, "VerboseScreenOutput", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::verbose), 3,
                               "Low", 0, "Medium", 1, "High", 2);
}

//---------------------------------------------------------

RestartData::RestartData()
{
  restart_file_num = 0;
  matrix_site_number = 0;
  solute_site_number = 0;
}

//---------------------------------------------------------

void RestartData::setup(const char *name)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, 0);

  new ClassInt<RestartData> (ca, "FileNum", this, &RestartData::restart_file_num); 

  new ClassInt<RestartData> (ca, "MatrixSiteNum", this, &RestartData::matrix_site_number);

  new ClassInt<RestartData> (ca, "SoluteSiteNum", this, &RestartData::solute_site_number); 

}

//-----------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  readCmdLine(argc, argv);
  readCmdFile();
}

//-----------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stderr,"ERROR: Input file not provided!\n");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//-----------------------------------------------------

void IoData::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stderr,"*** Error: could not open \'%s\'\n", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stderr,"*** Error: command file contained parsing errors\n");
    exit(error);
  }
  fclose(cmdFilePtr);
}

//-----------------------------------------------------

void IoData::setupCmdFileVariables()
{
  ClassAssigner *ca = new ClassAssigner(name, 7, 0);

  latticeMap.setup("Lattice", ca);

  materialMap.setup("Material", ca);

  icMap.setup("RegionalInitialCondition", ca);

  diffusion.seutp("Diffusion");

  ts.setup("Time");

  restart.setup("Restart");

  output.setup("Output");
}

//-----------------------------------------------------


//-----------------------------------------------------







