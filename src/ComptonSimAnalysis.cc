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
#include "TH1.h"
#include "TGraph.h"
#include "TString.h"
#include "TMath.h"
#include "TApplication.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"

// Boost libraries

#include "boost/filesystem.hpp"
#include "boost/math/constants/constants.hpp"
#include <boost/math/special_functions/erf.hpp>

// Custom Libraries

#include "ComptonSimAnalysis.hh"
#include "MsgStream.hh"
#include "SysMsg.hh"
#include "PhysicalConstants.hh"

#ifdef compton_sim_analysis_cxx

ComptonSimAnalysis::ComptonSimAnalysis()
{
  fGraphicsShow      = false;
  fFileSet           = false;
  fFileLeftSet       = false;
  fFileRightSet      = false;
  fComptonWeight     = false;
  fBackgroundWeight  = false;
  fHaloWeight        = false;
  fApertureSize      = false;
  fAsymmetryAnalysis = false;

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
  compton.gaussian_weight = 4.40643e-07;
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
     Sys::SysError << "Loading root file:\t" 
		   << fFileLocation 
		   << Sys::endl;
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
      std::string opt(options[i+1]);
      if(opt.compare("left"))
	 fPolarization = -1;
      if(opt.compare("right"))
	 fPolarization = 1;
    }
    if(flag.compare("--aperture-size") == 0){
      flag.clear();
      fApertureSize = true;
      compton.aperture_size = atof(options[i+1]);
    }
    if(flag.compare("--asymmetry") == 0){
      flag.clear();
      fAsymmetryAnalysis = true;
    }
   if(flag.compare("--file-right") == 0){
     std::string opt(options[i+1]);
     fFileRight = opt;
     flag.clear();
     fFileRightSet = true;
     Sys::SysMsg << "Loading root file:\t" 
		 << fFileRight 
		 << Sys::endl;
   }
   if(flag.compare("--file-left") == 0){
     std::string opt(options[i+1]);
     fFileLeft = opt;
     flag.clear();
     fFileLeftSet = true;
     Sys::SysMsg << "Loading root file:\t" 
		 << fFileLeft 
		 << Sys::endl;
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

void ComptonSimAnalysis::AsymmetryAnalysis()
{
  if(!(fFileLeftSet && fFileRightSet)){
    PrintError("In order to access asymetry analysis, two root files must be provided. Please use --help for more info.");
    exit(1);
  }

  double cs_right = CalculateIntegratedCS(1);
  double cs_left  = CalculateIntegratedCS(-1);

  TFile *file_left = new TFile(fFileLeft.c_str());
  TTree *tree_left = (TTree *)file_left->Get("flux");

  TFile *file_right = new TFile(fFileRight.c_str());
  TTree *tree_right = (TTree *)file_right->Get("flux");

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 800);
  canvas->SetTitle("Compton Detector Asymmetry");

  int lentries = tree_left->GetEntries();
  int rentries = tree_right->GetEntries();

  if((lentries == 0) || (rentries == 0)){
    Sys::SysError << "No events found. Exiting." << Sys::endl;
    exit(1);
  }
  TH1D *left = new TH1D("left", "left", 251, 0, 251);
  TH1D *right = new TH1D("right", "right", 251, 0, 251);

  canvas->cd();

  std::cout << "Histograming data." << std::endl;

  tree_left->Draw("id>>left", "id < 251", "goff");
  left = (TH1D *)gDirectory->Get("left");
  left->Draw();

  tree_right->Draw("id>>right", "id < 251", "goff");
  right = (TH1D *)gDirectory->Get("right");
  right->Draw();

  std::vector <int> nLeft(251);
  std::vector <int> nRight(251);
  std::vector <double> Error(251);

  std::fill(nLeft.begin(), nLeft.end(), 0);  
  std::fill(nRight.begin(), nRight.end(), 0);  

  std::vector <int> *id_l = 0;
  std::vector <int> *id_r = 0;

  int ID = 0;

  tree_left->ResetBranchAddresses();
  tree_left->SetBranchStatus("*", 0);
  tree_left->SetBranchStatus("id", 1);

  tree_left->SetBranchAddress("id", &id_l);

  tree_right->ResetBranchAddresses();
  tree_right->SetBranchStatus("*", 0);
  tree_right->SetBranchStatus("id", 1);

  tree_right->SetBranchAddress("id", &id_r);

  if(!(tree_left->GetLeaf("id"))){
    Sys::SysError << "Error finding leaf" << Sys::endl;
    exit(1);
  }

  if(!(tree_right->GetLeaf("id"))){
    Sys::SysError << "Error finding leaf" << Sys::endl;
    exit(1);
  }

  std::cout << (int)id_l->size() << std::endl;

  for(int i = 0; i < lentries; i++){
    tree_left->GetEntry(i);
    for(int j = 0; j < (int)id_l->size(); j++){
      if(id_l->at(j) < 201){
   	ID = id_l->at(j);
  	nLeft[ID]++;
      }
    }
  }
  for(int i = 0; i < rentries; i++){
    tree_right->GetEntry(i);
    for(int j = 0; j < (int)id_r->size(); j++){
      if(id_r->at(j) < 201){
   	ID = id_r->at(j);
   	nRight[ID]++;
      }
    }
  }

  TF1 *fn_l = new TF1("fn_1", "[0]", 1, 250);
  fn_l->SetParameter(0,1);

  TF1 *fn_r = new TF1("fn_r", "[0]", 1, 250);
  fn_r->SetParameter(0,1);

  TH1D *yield_l = (TH1D *)left->Clone("yield_l");
  TH1D *yield_r = (TH1D *)right->Clone("yield_r");

  left->Multiply(fn_l, cs_left);
  right->Multiply(fn_r, cs_right);

  TH1D *diff = (TH1D *)left->Clone("diff");
  TH1D *sum  = (TH1D *)left->Clone("sum");

  diff->Add(right, -1);
  sum->Add(right, 1);

  diff->Divide(sum);
  diff->SetStats(0);

  std::fstream asymmetry;
  std::fstream yield;

  asymmetry.open("asymmetry.dat", std::fstream::out);
  yield.open("yield.dat", std::fstream::out);

  if(!yield.is_open()){ 
    std::cout << "Failure." << std::endl;
    exit(1);
  }

  std::cout << "Strip\t|\tAsymmetry\t|\tBin Error\t|\tCalc Error\t|\tnL\t|\tnR" << std::endl;

  for(Int_t i = 0; i < 251; i++){

    Error[i] = 1/(TMath::Sqrt(nLeft[i] + nRight[i]));
    if((nLeft[i] == 0) && (nRight[i]) == 0) Error[i] = 0;
    
    asymmetry << i << " " 
	      << diff->GetBinContent(i) << " "
	      <<  Error[i] << "\n";

    std::cout << i << " " 
       	      << diff->GetBinContent(i) << " "
     	      << Error[i] << std::endl;
    yield << i << " " 
     	  << yield_l->GetBinContent(i) + yield_r->GetBinContent(i) << " "
     	  << 1/TMath::Sqrt(yield_l->GetBinError(i) + yield_r->GetBinError(i)) << "\n";
  }

  asymmetry.close();
  yield.close();

  std::cout << "Histograming asymmetry." << std::endl;

  diff->Draw(); 
  diff->SetLineColor(9);
  diff->SetLineWidth(2);
  diff->SetTitle("Compton Detector Asymmetry");

  canvas->SaveAs("asymmetry_nowindow.png");
  canvas->SaveAs("asymmetry_nowindow.C");
  
  Sys::SysMsg << "Analysis done." << Sys::endl;


  exit(0);
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

  Sys::SysMsg << "Luminosity: " << compton.luminosity << Sys::endl;

  return(compton.luminosity);
}

bool ComptonSimAnalysis::OpenFile(TChain *fChain)
{
  return true;
}

double ComptonSimAnalysis::CalculateGaussianWeight()
{
  if(!fApertureSize){
    Sys::SysError << "Aperture size not defined. Using default aperture size (+- 0.5 cm). Please refer to --help for instruction on setting this value." << Sys::endl;
    return 0;
  }
  // the factor of ten is a multplier for the halo and the '1e2' is to convert to cm
  compton.gaussian_weight = boost::math::erfc(compton.aperture_size/(std::sqrt(2)*10e2*beam.sigma_ey));
  Sys::SysMsg << "\n>>>>> Gaussian weight: " << compton.gaussian_weight << Sys::endl;

  return 0;
}

double ComptonSimAnalysis::CalculateHaloFraction()
{

  compton.halo_ratio = compton.halo_scale_x*compton.halo_scale_y*compton.halo_amplitude;
  return (compton.halo_ratio);
}

void ComptonSimAnalysis::InitGraphicsEngine(int Argc, char **Argv)
{
  Sys::SysMsg << "<<<< Initialize Graphics Engine." << Sys::endl;
  app = new TApplication("App", &Argc, Argv);

}

void ComptonSimAnalysis::RunGraphicsEngine()
{
  Sys::SysMsg << "<<<< Running Graphics Engine." << Sys::endl;
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

double ComptonSimAnalysis::CalculateIntegratedCS(double polarization)
{

  SetBeamEnergy(beam.beam_energy);
  SetLaserEnergy(beam.laser_energy);
  SetPolarization(polarization);
  Initialize(); 

  return(2*TMath::Pi()*this->GetFunction((char *)"cs")->Integral(0,1));
}

void ComptonSimAnalysis::PrintError(const char *message)
{
  Sys::SysError << message << "\n" << Sys::endl;
  return;
}

#endif

