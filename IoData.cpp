#include <IoData.h>
#include <Utils.h>
#include <cassert>
#include <climits>
#include <cstring>

RootClassAssigner *nullAssigner = new RootClassAssigner;

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

  new ClassDouble<ParallelepipedData> (ca, "Px", this, &ParallelepipedData::x0);
  new ClassDouble<ParallelepipedData> (ca, "Py", this, &ParallelepipedData::y0);
  new ClassDouble<ParallelepipedData> (ca, "Pz", this, &ParallelepipedData::z0);

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

void CustomGeometryData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 4, father);

  new ClassStr<CustomGeometryData> (ca, "GeometrySpecifier", this, &CustomGeometryData::specifier);

  new ClassToken<CustomGeometryData> (ca, "Inclusion", this,
     reinterpret_cast<int CustomGeometryData::*>(&CustomGeometryData::inclusion), 3,
     "Override", 0, "Intersection", 1, "Union", 2);

  new ClassInt<CustomGeometryData>(ca, "OperationOrder", this, &CustomGeometryData::order);

  initialConditions.setup("InitialState", ca);
}

//---------------------------------------------------------

GlobalSpeciesSetData::GlobalSpeciesSetData()
{ }

//---------------------------------------------------------

void
GlobalSpeciesSetData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);
  speciesMap.setup("Species", ca);
}

//---------------------------------------------------------

GlobalSpeciesData::GlobalSpeciesData()
{
  name = "";
}

//---------------------------------------------------------

Assigner *GlobalSpeciesData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 1, nullAssigner);

  new ClassStr<GlobalSpeciesData> (ca, "SpeciesName", this, &GlobalSpeciesData::name);

  return ca;
}

//---------------------------------------------------------

LocalSpeciesData::LocalSpeciesData()
{
  //These default values are "invalid" on purpose --> meaning the user did not specify them
  speciesID = -1;
  molar_fraction = -1.0e-3;
  min_molar_fraction = -1.0e-4;

  diffusive = FALSE;
}

//---------------------------------------------------------

Assigner *LocalSpeciesData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 4, nullAssigner);

  new ClassInt<LocalSpeciesData> (ca, "SpeciesID", this, &LocalSpeciesData::speciesID);

  new ClassToken<LocalSpeciesData>(ca, "Diffusive", this,
                                   reinterpret_cast<int LocalSpeciesData::*>(&LocalSpeciesData::diffusive), 2,
                                   "False", 0, "True", 1);

  new ClassDouble<LocalSpeciesData> (ca, "MolarFraction", this, &LocalSpeciesData::molar_fraction);

  new ClassDouble<LocalSpeciesData> (ca, "MinimumMolarFraction", this, &LocalSpeciesData::min_molar_fraction);

  return ca;

}

//---------------------------------------------------------

LatticeSiteData::LatticeSiteData()
{
  la = lb = lc = -1.0;
  materialid = -1;

  tag = 0;

  atomic_frequency = -1.0; //using an invalid value as default
  mean_displacement_x = mean_displacement_y = mean_displacement_z = 0.0;
  mean_momentum_x     = mean_momentum_y     = mean_momentum_z     = 0.0;
}

//---------------------------------------------------------

Assigner *LatticeSiteData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 13, nullAssigner);

  new ClassDouble<LatticeSiteData> (ca, "La", this, &LatticeSiteData::la);
  new ClassDouble<LatticeSiteData> (ca, "Lb", this, &LatticeSiteData::lb);
  new ClassDouble<LatticeSiteData> (ca, "Lc", this, &LatticeSiteData::lc);
  new ClassInt<LatticeSiteData> (ca, "MaterialID", this, &LatticeSiteData::materialid);

  new ClassInt<LatticeSiteData> (ca, "Tag", this, &LatticeSiteData::tag);

  new ClassDouble<LatticeSiteData> (ca, "AtomicFrequency", this, &LatticeSiteData::atomic_frequency);
  new ClassDouble<LatticeSiteData> (ca, "MeanDisplacementX", this, &LatticeSiteData::mean_displacement_x);
  new ClassDouble<LatticeSiteData> (ca, "MeanDisplacementY", this, &LatticeSiteData::mean_displacement_y);
  new ClassDouble<LatticeSiteData> (ca, "MeanDisplacementZ", this, &LatticeSiteData::mean_displacement_z);
  new ClassDouble<LatticeSiteData> (ca, "MeanMomentumX", this, &LatticeSiteData::mean_momentum_x);
  new ClassDouble<LatticeSiteData> (ca, "MeanMomentumY", this, &LatticeSiteData::mean_momentum_y);
  new ClassDouble<LatticeSiteData> (ca, "MeanMomentumZ", this, &LatticeSiteData::mean_momentum_z);
  
  speciesMap.setup("LocalSpecies", ca);

  return ca;
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

  ClassAssigner *ca = new ClassAssigner("normal", 24, 0);

  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  parallelepipedMap.setup("Parallelepiped", ca);
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
  
  return ca;
}

//---------------------------------------------------------

InteratomicPotentialData::InteratomicPotentialData() 
{
  filename = "";
  eps = 1.0e-5; 
  rc = 5.35; 
}

//---------------------------------------------------------

void InteratomicPotentialData::setup(const char *name, ClassAssigner *father)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassStr<InteratomicPotentialData> (ca, "File", this, &InteratomicPotentialData::filename);

  new ClassDouble<InteratomicPotentialData> (ca, "Tolerance", this, &InteratomicPotentialData::eps);

  new ClassDouble<InteratomicPotentialData> (ca, "CutOffDistance", this, &InteratomicPotentialData::rc);

}

//---------------------------------------------------------

MaterialData::MaterialData()
{
  name = "";

  species0 = species1 = species2 = species3 = species4 = species5 = species6 = 
  species7 = species8 = species9 = species10 = species11 = -1; //invalid value as default

  min_molar_fraction = 0.0;
}

//---------------------------------------------------------

Assigner *MaterialData::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 15, 0);

  new ClassStr<MaterialData> (ca, "MaterialName", this, &MaterialData::name);

  new ClassInt<MaterialData> (ca, "Species0", this, &MaterialData::species0); 
  new ClassInt<MaterialData> (ca, "Species1", this, &MaterialData::species1); 
  new ClassInt<MaterialData> (ca, "Species2", this, &MaterialData::species2); 
  new ClassInt<MaterialData> (ca, "Species3", this, &MaterialData::species3); 
  new ClassInt<MaterialData> (ca, "Species4", this, &MaterialData::species4); 
  new ClassInt<MaterialData> (ca, "Species5", this, &MaterialData::species5); 
  new ClassInt<MaterialData> (ca, "Species6", this, &MaterialData::species6); 
  new ClassInt<MaterialData> (ca, "Species7", this, &MaterialData::species7); 
  new ClassInt<MaterialData> (ca, "Species8", this, &MaterialData::species8); 
  new ClassInt<MaterialData> (ca, "Species9", this, &MaterialData::species9); 
  new ClassInt<MaterialData> (ca, "Species10", this, &MaterialData::species10); 
  new ClassInt<MaterialData> (ca, "Species11", this, &MaterialData::species11); 

  new ClassDouble<MaterialData> (ca, "MinimumMolarFraction", this, &MaterialData::min_molar_fraction);

  potential.setup("InteratomicPotential");

  return ca;
}

//---------------------------------------------------------

RegionalIcData::RegionalIcData()
{
  // a negative value means the user did not specify it
  latticeID = -1;
  materialID = -1;
  siteID = -1;

  tag = 0;

  atomic_frequency = -1.0; //using an invalid value as default
  mean_displacement_x = mean_displacement_y = mean_displacement_z = 0.0;
  mean_momentum_x     = mean_momentum_y     = mean_momentum_z     = 0.0;
}

//---------------------------------------------------------

Assigner *RegionalIcData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 18, 0);

  planeMap.setup("Plane", ca);
  sphereMap.setup("Sphere", ca);
  parallelepipedMap.setup("Parallelepiped", ca);
  spheroidMap.setup("Spheroid", ca);
  cylinderconeMap.setup("CylinderAndCone", ca);
  cylindersphereMap.setup("CylinderWithSphericalCaps", ca);

  custom_geometry.setup("CustomGeometry", ca);

  new ClassInt<RegionalIcData> (ca, "LatticeID", this, &RegionalIcData::latticeID);
  new ClassInt<RegionalIcData> (ca, "MaterialID", this, &RegionalIcData::materialID);
  new ClassInt<RegionalIcData> (ca, "SiteID", this, &RegionalIcData::siteID);

  new ClassInt<RegionalIcData> (ca, "Tag", this, &RegionalIcData::tag);
  new ClassDouble<RegionalIcData> (ca, "AtomicFrequency", this, &RegionalIcData::atomic_frequency);
  new ClassDouble<RegionalIcData> (ca, "MeanDisplacementX", this, &RegionalIcData::mean_displacement_x);
  new ClassDouble<RegionalIcData> (ca, "MeanDisplacementY", this, &RegionalIcData::mean_displacement_y);
  new ClassDouble<RegionalIcData> (ca, "MeanDisplacementZ", this, &RegionalIcData::mean_displacement_z);
  new ClassDouble<RegionalIcData> (ca, "MeanMomentumX", this, &RegionalIcData::mean_momentum_x);
  new ClassDouble<RegionalIcData> (ca, "MeanMomentumY", this, &RegionalIcData::mean_momentum_y);
  new ClassDouble<RegionalIcData> (ca, "MeanMomentumZ", this, &RegionalIcData::mean_momentum_z);

  speciesMap.setup("LocalSpecies", ca);

  return ca;
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

void DiffusionData::setup(const char *name, ClassAssigner *father)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 5, father);

  new ClassDouble<DiffusionData> (ca, "Temperature", this, &DiffusionData::T0);

  new ClassDouble<DiffusionData> (ca, "AttemptFrequency", this, &DiffusionData::v);
  new ClassDouble<DiffusionData> (ca, "ActivationEnergy", this, &DiffusionData::Qm);

  new ClassDouble<DiffusionData> (ca, "SubsurfaceThickness", this, &DiffusionData::t_subsurf);
  new ClassDouble<DiffusionData> (ca, "SoluteChemicalPotential", this, &DiffusionData::mubd);

}

//---------------------------------------------------------


ExplicitData::ExplicitData()
{
  type = RUNGE_KUTTA_2;
}

//---------------------------------------------------------

void ExplicitData::setup(const char *name, ClassAssigner *father)
{

 ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<ExplicitData>
    (ca, "Type", this,
     reinterpret_cast<int ExplicitData::*>(&ExplicitData::type), 3,
     "ForwardEuler", 0, "RungeKutta2", 1, "RungeKutta3", 2);

}

//---------------------------------------------------------

TsData::TsData()
{
  type = EXPLICIT;
  maxIts = INT_MAX;
  timestep = -1.0;
  cfl = 0.5;
  maxTime = 1e10;
}

//---------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<TsData>(ca, "Type", this,
                         reinterpret_cast<int TsData::*>(&TsData::type), 2,
                         "Explicit", 0, "Implicit", 1);
  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "CFL", this, &TsData::cfl);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);

  expl.setup("Explicit", ca);
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
  precision = MEDIUM;
}

//---------------------------------------------------------

void OutputData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 7, father);

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

  new ClassToken<OutputData>(ca, "Precision", this,
                               reinterpret_cast<int OutputData::*>(&OutputData::precision), 3,
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

void RestartData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

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
 
  // read the first input argument that is not a "flag"
  for(int i=1; i<argc; i++) {
    if(!strcmp(argv[i], "-t"))
      continue;
    cmdFileName = argv[i];
    break;
  }
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
  ClassAssigner *ca = new ClassAssigner("normal", 8, 0);

  latticeMap.setup("Lattice", 0);

  materialMap.setup("Material", 0);

  speciesMap.setup("Species", 0);

  icMap.setup("RegionalInitialCondition", ca);

  diffusion.setup("Diffusion");

  ts.setup("Time");

  restart.setup("Restart");

  output.setup("Output");
}

//-----------------------------------------------------

void IoData::finalize()
{
  assert(MaterialData::MAX_SPECIES==12);
  for(auto&& mat : materialMap.dataMap) {
    mat.second->species[0] = mat.second->species0;
    mat.second->species[1] = mat.second->species1;
    mat.second->species[2] = mat.second->species2;
    mat.second->species[3] = mat.second->species3;
    mat.second->species[4] = mat.second->species4;
    mat.second->species[5] = mat.second->species5;
    mat.second->species[6] = mat.second->species6;
    mat.second->species[7] = mat.second->species7;
    mat.second->species[8] = mat.second->species8;
    mat.second->species[9] = mat.second->species9;
    mat.second->species[10] = mat.second->species10;
    mat.second->species[11] = mat.second->species11;

    // figure out nSpecies
    mat.second->nSpecies = 0;
    for(int i=0; i<MaterialData::MAX_SPECIES; i++) {
      if(mat.second->species[i]>=0)
        mat.second->nSpecies++;
      else
        break;
    }
    if(mat.second->nSpecies==0) {
      print_error("*** Error: Material[%d] does not have any species.\n", mat.first);
      exit_mpi();
    }
    for(int i=mat.second->nSpecies; i<MaterialData::MAX_SPECIES; i++) {
      if(mat.second->species[i]>=0) {
        print_error("*** Error: Found gap(s) in Material[%d] species indices.\n", mat.first);
        exit_mpi();
      }
    }
  }
}

//-----------------------------------------------------







