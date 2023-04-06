#include <iostream>
#include <IoData.h>

//---------------------------------------------------------

PlaneData::PlaneData()
{

  cen_x  = 0.0;
  cen_y  = 0.0;
  cen_z  = 0.0;
  nx     = 0.0;
  ny     = 0.0;
  nz     = 0.0;

}

//------------------------------------------------------------------------------

Assigner *PlaneData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 6, nullAssigner);

  new ClassDouble<PlaneData> (ca, "Point_x", this, &PlaneData::cen_x);
  new ClassDouble<PlaneData> (ca, "Point_y", this, &PlaneData::cen_y);
  new ClassDouble<PlaneData> (ca, "Point_z", this, &PlaneData::cen_z);
  new ClassDouble<PlaneData> (ca, "Normal_x", this, &PlaneData::nx);
  new ClassDouble<PlaneData> (ca, "Normal_y", this, &PlaneData::ny);
  new ClassDouble<PlaneData> (ca, "Normal_z", this, &PlaneData::nz);

  return ca;
}

//---------------------------------------------------------

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

}

//------------------------------------------------------------------------------

Assigner *SpheroidData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 8, nullAssigner);

  new ClassDouble<SpheroidData> (ca, "Center_x", this, &SpheroidData::cen_x);
  new ClassDouble<SpheroidData> (ca, "Center_y", this, &SpheroidData::cen_y);
  new ClassDouble<SpheroidData> (ca, "Center_z", this, &SpheroidData::cen_z);
  new ClassDouble<SpheroidData> (ca, "Axis_x", this, &SpheroidData::axis_x);
  new ClassDouble<SpheroidData> (ca, "Axis_y", this, &SpheroidData::axis_y);
  new ClassDouble<SpheroidData> (ca, "Axis_z", this, &SpheroidData::axis_z);
  new ClassDouble<SpheroidData> (ca, "Length", this, &SpheroidData::length);
  new ClassDouble<SpheroidData> (ca, "Diameter", this, &SpheroidData::diameter);

  return ca;
}

//---------------------------------------------------------

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
}

//------------------------------------------------------------------------------

Assigner *CylinderConeData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 10, nullAssigner);

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

  return ca;
}

//---------------------------------------------------------

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

}

//------------------------------------------------------------------------------

Assigner *CylinderSphereData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 10, nullAssigner);

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

  return ca;
}

//---------------------------------------------------------

LatticeDomainData::LatticeDomainData()
{
  custom_geometry_specifier = "";
}

//---------------------------------------------------------

void LatticeDomainData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 5, 0);

  planeMap.setup("Plane", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);

  new ClassStr<LatticeDomainData> (ca, "GeometrySpecifier", this, &LatticeDomainData::custom_geometry_specifier);

}

//---------------------------------------------------------

SpaceData::SpaceData() 
{
  lattice_type = FCC;
  interstice_type = NONE;
  sample_shape = CUBE;
  N = 6;
  size = -1.0;
  a = -1.0;
  b = -1.0;
  c = -1.0; 
  alpha = 90.0;
  beta  = 90.0;
  gamma = 90.0;
  
  xe = 0.0; 
  sigma_0 = 1000.0;
  sigma_1 = 1000.0;
}

//---------------------------------------------------------

void SpaceData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 15, 0);

  new ClassToken<SpaceData>(ca, "LatticeType", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::lattice_type), 9,
                               "SimpleCubic", 0, "SC", 0, 
                               "BodyCenteredCubic", 1, "BCC", 1,
                               "FaceCenteredCubic", 2, "FCC", 2,
                               "HexagonalClosePacked", 3, "HCP", 3,
                               "Custom", 4);

  new ClassToken<SpaceData>(ca, "IntersticeType", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::interstice_type), 5,
                               "None", 0, "Central", 1, "Octahedral", 2, "Tetrahedral", 3, 
                               "OctahedralAndTetrahedral", 4);

  new ClassToken<SpaceData>(ca, "SampleShape", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::sample_shape), 4,
                               "Cube", 0, "Sphere", 1, "Octahedron", 2, "Custom", 3);

  new ClassInt<SpaceData>(ca, "NumberOfUnitCells", this, &SpaceData::N);

  new ClassDouble<SpaceData> (ca, "SampleSize", this, &SpaceData::size);

  new ClassDouble<SpaceData> (ca, "SoluteMolarFraction", this, &SpaceData::xe);

  new ClassDouble<SpaceData> (ca, "LatticeConstant", this, &SpaceData::a);

  new ClassDouble<SpaceData> (ca, "LatticeConstant1", this, &SpaceData::a);

  new ClassDouble<SpaceData> (ca, "LatticeConstant2", this, &SpaceData::b);

  new ClassDouble<SpaceData> (ca, "LatticeConstant3", this, &SpaceData::c);

  new ClassDouble<SpaceData> (ca, "LatticeAngle23", this, &SpaceData::alpha);

  new ClassDouble<SpaceData> (ca, "LatticeAngle13", this, &SpaceData::beta);

  new ClassDouble<SpaceData> (ca, "LatticeAngle12", this, &SpaceData::gamma);

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
  potential.setup("InteratomicPotential");

  space.setup("Space");

  diffusion.seutp("Diffusion");

  ts.setup("Time");

  restart.setup("Restart");

  output.setup("Output");
}

//-----------------------------------------------------


//-----------------------------------------------------







