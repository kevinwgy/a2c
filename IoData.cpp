#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <IoData.h>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>

//---------------------------------------------------------

SpaceData::SpaceData() 
{
  lattice_type = FCC;
  interstice_type = NONE;
  sample_shape = CUBE;
  N = 6;
  size = -1.0;
  lattice_constant = 3.8;
  xe = 0.0; 
  sigma_0 = 1000.0;
  sigma_1 = 1000.0;
}

//---------------------------------------------------------

void SpaceData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 9, 0);

  new ClassToken<SpaceData>(ca, "LatticeType", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::lattice_type), 8,
                               "SimpleCubic", 0, "SC", 0, 
                               "BodyCenteredCubic", 1, "BCC", 1,
                               "FaceCenteredCubic", 2, "FCC", 2,
                               "HexagonalClosePacked", 3, "HCP", 3);

  new ClassToken<SpaceData>(ca, "IntersticeType", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::interstice_type), 4,
                               "None", 0, "Octahedral", 1, "Tetrahedral", 2, 
                               "OctahedralAndTetrahedral", 3);

  new ClassToken<SpaceData>(ca, "SampleShape", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::sample_shape), 4,
                               "Cube", 0, "Sphere", 1, "Octahedron", 2, "Other", 3);

  new ClassInt<SpaceData>(ca, "NumberOfUnitCells", this, &SpaceData::N);

  new ClassDouble<SpaceData> (ca, "SampleSize", this, &SpaceData::size);

  new ClassDouble<SpaceData> (ca, "SoluteMolarFraction", this, &SpaceData::xe);

  new ClassDouble<SpaceData> (ca, "LatticeConstant", this, &SpaceData::lattice_constant);

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







