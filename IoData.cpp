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
  sample_shape = cube;
  N = 6;
  size = -1.0;
  lattice_parameter = 3.8;
  sigma_0 = 1000.0;
  sigma_1 = 1000.0;
}

//---------------------------------------------------------

void SpaceData::setup(const char *name)
{
 
  ClassAssigner *ca = new ClassAssigner(name, 24, 0);

  new ClassToken<MeshData>(ca, "LatticeType", this,
                               reinterpret_cast<int SpaceData::*>(&SpaceData::lattice_type), 8,
                               "SimpleCubic", 0, "SC", 0, 
                               "BodyCenteredCubic", 1, "BCC", 1,
                               "FaceCenteredCubic", 2, "FCC", 2,
                               "Cylindrical", 2);
}

//---------------------------------------------------------

  // set defaults
  T0 = 300; 
  eps = 1.0e-5; 
  rc = 5.35; 
  N = 6; 
  xe = 0.0; 
  a_Pd = 3.8;
  sigma_Pd = 1000.0; 
  sigma_H = 1000.0; 
  
  dt = 1.0e-6; 
  t_final = 0.01;

  minimize_frequency = 10;
  output_frequency = 100;

  mubd = -2.0; // (Pa)

  foldername = "";
  filename_base = "";
 
  restart = 0;
  restart_file_num = 0;
}


InputFileData::~InputFileData()
{
//  delete[] foldername;
//  delete[] filename_base;
}


void InputFileData::setup(const char *name)
{
  ClassAssigner *ca = new ClassAssigner(name, 24, 0); //father);

  new ClassStr<InputFileData> (ca, "ResultFolder", this, &InputFileData::foldername);
//  sprintf(foldername,"%s%s", foldername, "/"); // append a slash
  new ClassStr<InputFileData> (ca, "ResultFilePrefix", this, &InputFileData::filename_base);

  new ClassDouble<InputFileData> (ca, "Temperature", this, &InputFileData::T0);
  new ClassDouble<InputFileData> (ca, "CutOffDistance", this, &InputFileData::rc);
  
  new ClassDouble<InputFileData> (ca, "Eps", this, &InputFileData::eps);
  new ClassDouble<InputFileData> (ca, "FractionUpperBound", this, &InputFileData::xUb);
  new ClassDouble<InputFileData> (ca, "FractionLowerBound", this, &InputFileData::xLb);

  new ClassInt<InputFileData> (ca, "SampleShapeType", this, &InputFileData::sample_shape);
  new ClassInt<InputFileData> (ca, "NumberOfUnitCells", this, &InputFileData::N);
  
  new ClassDouble<InputFileData> (ca, "InitialHFraction", this, &InputFileData::xe);
  new ClassDouble<InputFileData> (ca, "InitialPdLatticeConstant", this, &InputFileData::a_Pd);
  new ClassDouble<InputFileData> (ca, "InitialSigmaPd", this, &InputFileData::sigma_Pd);
  new ClassDouble<InputFileData> (ca, "InitialSigmaH", this, &InputFileData::sigma_H);

  new ClassDouble<InputFileData> (ca, "TimeStepSize", this, &InputFileData::dt);
  new ClassDouble<InputFileData> (ca, "FinalTime", this, &InputFileData::t_final);

  new ClassInt<InputFileData> (ca, "OutputFrequency", this, &InputFileData::output_frequency);

  new ClassDouble<InputFileData> (ca, "HydrogenChemicalPotential", this, &InputFileData::mubd);
  new ClassDouble<InputFileData> (ca, "SubsurfaceThickness", this, &InputFileData::t_subsurf);

  new ClassDouble<InputFileData> (ca, "AttemptFrequency", this, &InputFileData::v);
  new ClassDouble<InputFileData> (ca, "ActivationEnergy", this, &InputFileData::Qm);

  new ClassInt<InputFileData> (ca, "Restart", this, &InputFileData::restart);
  new ClassInt<InputFileData> (ca, "RestartFileNum", this, &InputFileData::restart_file_num); 
  new ClassInt<InputFileData> (ca, "PdSiteNum", this, &InputFileData::nPd_restart);
  new ClassInt<InputFileData> (ca, "HSiteNum", this, &InputFileData::nH_restart); 
}

//-----------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  readCmdLine(argc, argv);
  readCmdFile();
}

//-----------------------------------------------------

void Input::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stderr,"ERROR: Input file not provided!\n");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//-----------------------------------------------------

void Input::readCmdFile()
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

void Input::setupCmdFileVariables()
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







