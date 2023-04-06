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

}

//------------------------------------------------------------------------------

Assigner *PlaneData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);

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
}

//------------------------------------------------------------------------------

Assigner *SphereData::getAssigner()
{
  
  ClassAssigner *ca = new ClassAssigner("normal", 7, nullAssigner);
  
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

}

//------------------------------------------------------------------------------

Assigner *ParallelepipedData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 15, nullAssigner);

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

  length = 0.0;
  diameter = 0.0;

  side = INTERIOR;
  inclusion = OVERRIDE;
}

//------------------------------------------------------------------------------

Assigner *SpheroidData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 11, nullAssigner);

  new ClassDouble<SpheroidData> (ca, "Center_x", this, &SpheroidData::cen_x);
  new ClassDouble<SpheroidData> (ca, "Center_y", this, &SpheroidData::cen_y);
  new ClassDouble<SpheroidData> (ca, "Center_z", this, &SpheroidData::cen_z);
  new ClassDouble<SpheroidData> (ca, "Axis_x", this, &SpheroidData::axis_x);
  new ClassDouble<SpheroidData> (ca, "Axis_y", this, &SpheroidData::axis_y);
  new ClassDouble<SpheroidData> (ca, "Axis_z", this, &SpheroidData::axis_z);
  new ClassDouble<SpheroidData> (ca, "Length", this, &SpheroidData::length);
  new ClassDouble<SpheroidData> (ca, "Diameter", this, &SpheroidData::diameter);

  new ClassToken<SpheroidData> (ca, "Side", this,
     reinterpret_cast<int SpheroidData::*>(&SpheroidData::side), 2,
     "Interior", 0, "Exterior", 1);

  new ClassToken<SpheroidData> (ca, "Inclusion", this,
     reinterpret_cast<int SpheroidData::*>(&SpheroidData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

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
}

//------------------------------------------------------------------------------

Assigner *CylinderConeData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 13, nullAssigner);

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
}

//------------------------------------------------------------------------------

Assigner *CylinderSphereData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 13, nullAssigner);

  new ClassDouble<CylinderSphereData> (ca, "Axis_x", this, &CylinderSphereData::nx);
  new ClassDouble<CylinderSphereData> (ca, "Axis_y", this, &CylinderSphereData::ny);
  new ClassDouble<CylinderSphereData> (ca, "Axis_z", this, &CylinderSphereData::nz);
  new ClassDouble<CylinderSphereData> (ca, "CylinderCenter_x", this, &CylinderSphereData::cen_x);
  new ClassDouble<CylinderSphereData> (ca, "CylinderCenter_y", this, &CylinderSphereData::cen_y);
  new ClassDouble<CylinderSphereData> (ca, "CylinderCenter_z", this, &CylinderSphereData::cen_z);
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

  initialConditions.setup("InitialState", ca);

  return ca;
}

//---------------------------------------------------------

LatticeSiteData::LatticeSiteData()
{
  la = lb = lc = -1.0;
  group_id = -1;
}

//---------------------------------------------------------

Assigner *LatticeSiteData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);

  new ClassDouble<LatticeSiteData> (ca, "La", this, &LatticeSiteData::la);
  new ClassDouble<LatticeSiteData> (ca, "Lb", this, &LatticeSiteData::lb);
  new ClassDouble<LatticeSiteData> (ca, "Lc", this, &LatticeSiteData::lc);
  new ClassInt<LatticeSiteData> (ca, "Group", this, &LatticeSiteData::group_id);
}

//---------------------------------------------------------

LatticeData::LatticeData()
{
  custom_geometry_specifier = "";
}

//---------------------------------------------------------

void LatticeData::setup(const char *name)
{

  ClassAssigner *ca = new ClassAssigner(name, 21, 0);

  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  parallelepiped.setup("Parallelepiped", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);
  new ClassStr<LatticeDomainData> (ca, "GeometrySpecifier", this, &LatticeDomainData::custom_geometry_specifier);

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
  
  siteMap.setup("Site", ca);

  new ClassDouble<LatticeData> (ca, "MinDistance", this, &LatticeData::dmin);
  
}

//---------------------------------------------------------

IcData::IcData()
{
  xe = 0.0; 
  sigma_0 = 1000.0;
  sigma_1 = 1000.0;
}

//---------------------------------------------------------

void SpaceData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 3, 0);

  new ClassDouble<SpaceData> (ca, "SoluteMolarFraction", this, &SpaceData::xe);
  new ClassDouble<SpaceData> (ca, "Sigma0", this, &SpaceData::sigma_0);
  new ClassDouble<SpaceData> (ca, "Sigma1", this, &SpaceData::sigma_1);
  
}

//---------------------------------------------------------

PotentialData::PotentialData() 
{
  filename = "";
  eps = 1.0e-5; 
  rc = 5.35; 
}

//---------------------------------------------------------

void PotentialData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 3, 0);

  new ClassStr<PotentialData> (ca, "File", this, &PotentialData::filename);

  new ClassDouble<PotentialData> (ca, "Tolerance", this, &PotentialData::eps);

  new ClassDouble<PotentialData> (ca, "CutOffDistance", this, &PotentialData::rc);

}

//---------------------------------------------------------

DiffusionData::DiffusionData()
{
  T0 = 300; 
  v = 0.0;
  Qm = 0.0;
  xUb = 1.0;
  xLb = 0.0;
  t_subsurf = 0.0;
  mubd = -2.0; // dimensional (pressure)
}

//---------------------------------------------------------

void DiffusionData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 7, 0);

  new ClassDouble<DiffusionData> (ca, "Temperature", this, &DiffusionData::T0);

  new ClassDouble<DiffusionData> (ca, "AttemptFrequency", this, &DiffusionData::v);
  new ClassDouble<DiffusionData> (ca, "ActivationEnergy", this, &DiffusionData::Qm);

  new ClassDouble<DiffusionData> (ca, "FractionUpperBound", this, &DiffusionData::xUb);
  new ClassDouble<DiffusionData> (ca, "FractionLowerBound", this, &DiffusionData::xLb);

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
  foldername = "";
  filename_base = "";
  frequency = 100;
}

//---------------------------------------------------------

void OutputData::setup(const char *name)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, 0);

  new ClassStr<OutputData> (ca, "Path", this, &OutputData::foldername);

  new ClassStr<OutputData> (ca, "FilePrefix", this, &OutputData::filename_base);

  new ClassInt<OutputData> (ca, "Frequency", this, &OutputData::frequency);

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

  potential.setup("Potential");

  ic.setup("InitialCondition");

  diffusion.seutp("Diffusion");

  ts.setup("Time");

  restart.setup("Restart");

  output.setup("Output");
}

//-----------------------------------------------------


//-----------------------------------------------------







