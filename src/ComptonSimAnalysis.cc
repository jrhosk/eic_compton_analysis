#define compton_sim_analysis_cxx
// C++ Libraries

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>

// Root libraries

#include "TAxis.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"
#include "TMath.h"
#include "TApplication.h"
#include "TChain.h"

// Boost libraries

#include "boost/filesystem.hpp"
#include "boost/math/constants/constants.hpp"
#include <boost/math/special_functions/erf.hpp>

// Custom Libraries

#include "ComptonSimAnalysis.hh"
#include "FontColor.hh"
#include "PhysicalConstants.hh"

#ifdef compton_sim_analysis_cxx

ComptonSimAnalysis::ComptonSimAnalysis()
{
  fGraphicsShow     = false;
  fFileSet          = false;
  fComptonWeight    = false;
  fBackgroundWeight = false;
  fHaloWeight       = false;
  fApertureSize     = false;

  beam.beam_energy  = 5;         // GeV
  beam.laser_energy = 2.33e-9;   // GeV
  beam.polarization = -1;
  beam.sigma_ex     = 226.6e-6;  // meters
  beam.sigma_ey     = 99e-6;     // meters 
  beam.sigma_g      = 151.4e-6; 
  beam.strip_number = 200;

  compton.halo_amplitude = 7.2e-5;
  compton.halo_scale_x = 3.3;
  compton.halo_scale_y = 10;
  compton.gaussian_weight = 4.40643e-07;     // Should really be set each case.
  compton.aperture_size = 0.5;

  fFileOutput = "analysis.C";

}

ComptonSimAnalysis::~ComptonSimAnalysis()
{
  // Filler
}



void ComptonSimAnalysis::GetOptions(char **options){

  int i = 0;
  
  std::string flag;

 while(options[i] != NULL){
   flag = options[i];

   if(flag.compare("--file") == 0){
     std::string opt(options[i+1]);
     fFileLocation = opt;
     flag.clear();
     fFileSet = true;
     std::cout << red << "Loading root file:\t" 
    	       << fFileLocation 
    	       << white << std::endl;
   }
   if(flag.compare("--output-file") == 0){
     std::string opt(options[i+1]);
     fFileOutput = opt;
     flag.clear();
   }
    if(flag.compare("--graphics") == 0){
      flag.clear();
      fGraphicsShow = true;
    }
    if(flag.compare("--compton-weight") == 0){
      flag.clear();
      fComptonWeight = true;
    }
    if(flag.compare("--background-weight") == 0){
      flag.clear();
      fBackgroundWeight = true;
    }
    if(flag.compare("--halo-weight") == 0){
      flag.clear();
      fHaloWeight = true;
    }
    if(flag.compare("--polarization") == 0){
      flag.clear();
      fPolarization = atof(options[i+1]);
    }
    if(flag.compare("--aperture-size") == 0){
      flag.clear();
      fApertureSize = true;
      compton.aperture_size = atof(options[i+1]);
    }
    if(flag.compare("--number-strips") == 0){
      flag.clear();
      beam.strip_number = atof(options[i+1]);
    }
    if(flag.compare("--energy") == 0){
      flag.clear();
      beam.beam_energy = atof(options[i+1]);
      if(beam.beam_energy == 3) {
	beam.sigma_ex = 136e-6;
	beam.sigma_ey = 56e-6; 
      }
      if(beam.beam_energy == 5) {
	beam.sigma_ex = 226.6e-6;
	beam.sigma_ey = 99e-6; 
      }
      if(beam.beam_energy == 11) {
	beam.sigma_ex = 356e-6;
	beam.sigma_ey = 115e-6;
      }
    }

   if(flag.compare("--help") == 0){
     printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
     printf("Usage: ./analysis <options>\n");
     printf("         --graphics       \tGraphical output.\n");
     printf("         --file           \tDefines input root file.\n");
     printf("         --output-file    \tDefines output root file.\n");
     printf("         --compton-weight               \tWeighting flag.\n");
     printf("         --halo-weight                  \tWeighting flag.\n");
     printf("         --background-weight            \tWeighting flag.\n");
     printf("         --polarization  <polarization> \tBeam polarization.\n");
     printf("         --energy        <energy GeV>   \tBeam energy.\n");
     printf("         --aperture-size <size cm>      \tAperture size.\n");
     printf("         --number-strips <# strips>     \tNumber of detector strips.\n");
     printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
     exit(0);
   }
   i++;
 }
}

double ComptonSimAnalysis::CalculateLuminosity(double current)
{
  double wavelength = 532e-9;
  double crossing_angle = 2.58*(boost::math::constants::pi<double>()/180);
  double laser_power = 1e3;  // Joules/s

  double coeff = (1+std::cos(crossing_angle))/(boost::math::constants::root_two_pi<double>());
  double top   = (current*wavelength*laser_power);
  double bottom = (pchep::e_SI*pchep::hc*pchep::c_light*std::sin(crossing_angle)*std::sqrt(std::pow(beam.sigma_ex,2) + std::pow(beam.sigma_g,2)));
  
  compton.luminosity = (1e-4)*(coeff*top)/bottom;  // convert to cm^-2

  std::cout << blue << "Luminosity: " << compton.luminosity << white << std::endl;

  return(compton.luminosity);
}

bool ComptonSimAnalysis::OpenFile(TChain *fChain)
{
  return true;
}

double ComptonSimAnalysis::CalculateGaussianWeight()
{
  if(!fApertureSize){
    PrintError("Aperture size not defined. Using default aperture size (+- 0.5 cm). Please refer to --help for instruction on setting this value.");
    return 0;
  }
  // the factor of ten is a multplier for the halo and the 'e2' is to convert to cm
  compton.gaussian_weight = boost::math::erfc(compton.aperture_size/(std::sqrt(2)*10e2*beam.sigma_ey));
  std::cout << blue << "\n>>>>> Gaussian weight: " << compton.gaussian_weight << white << std::endl;

  return 0;
}

double ComptonSimAnalysis::CalculateHaloFraction()
{

  compton.halo_ratio = compton.halo_scale_x*compton.halo_scale_y*compton.halo_amplitude;
  return (compton.halo_ratio);
}

void ComptonSimAnalysis::InitGraphicsEngine(int Argc, char **Argv)
{
  std::cout << green << "<<<< Initialize Graphics Engine." << white << std::endl;
  app = new TApplication("App", &Argc, Argv);

}

void ComptonSimAnalysis::RunGraphicsEngine()
{
  std::cout << green << "<<<< Running Graphics Engine." << white << std::endl;
  app->Run();
}

double ComptonSimAnalysis::CalculateIntegratedCS()
{

  SetBeamEnergy(beam.beam_energy);
  SetLaserEnergy(beam.laser_energy);
  SetPolarization(beam.polarization);

  Initialize(); 
  compton.cross_section = 2*TMath::Pi()*this->GetFunction((char *)"cs")->Integral(0,1);

  return(2*TMath::Pi()*this->GetFunction((char *)"cs")->Integral(0,1));

}

void ComptonSimAnalysis::PrintError(const char *message)
{
  std::cout << red << message << "\n" << white << std::endl;
  return;
}

#endif

