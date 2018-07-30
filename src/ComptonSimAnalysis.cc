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
#include "boost/math/special_functions/erf.hpp"

// Custom Libraries

#include "ComptonSimAnalysis.hh"
#include "MsgStream.hh"
#include "SysMsg.hh"
#include "PhysicalConstants.hh"

#ifdef compton_sim_analysis_cxx

ComptonSimAnalysis::ComptonSimAnalysis()
{
  fGraphicsShow       = false;
  fFileSet            = false;
  fFileLeftSet        = false;
  fFileRightSet       = false;
  fComptonWeight      = false;
  fBackgroundWeight   = false;
  fHaloWeight         = false;
  fApertureSize       = false;
  fAsymmetryAnalysis  = false;
  fResolutionAnalysis = false;
  fVetrocAnalysis     = false;

  beam.beam_energy  = 5;         // GeV
  beam.laser_energy = 2.33e-9;   // GeV
  beam.polarization = -1;
  beam.sigma_ex     = 226.6e-6;  // meters
  beam.sigma_ey     = 99e-6;     // meters 
  beam.sigma_g      = 151.4e-6; 
  beam.strip_number = 10300;

  compton.halo_amplitude = 7.2e-5;
  compton.halo_scale_x = 10.0;
  compton.halo_scale_y = 10.0;
  compton.gaussian_weight = 4.40643e-07;
  compton.aperture_size = 0.8;

  fFileOutput = "analysis";
  fEnergyCut = 0;                // MeV

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
      fFileOutput = "compton";
    }
    if(flag.compare("--background-weight") == 0){
      flag.clear();
      fBackgroundWeight = true;
      fFileOutput = "background";
    }
    if(flag.compare("--halo-weight") == 0){
      flag.clear();
      fHaloWeight = true;
      fFileOutput = "halo";
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
    if(flag.compare("--energy-cut") == 0){
      flag.clear();
      fEnergyCut = atof(options[i+1]);

      Sys::SysMsg << "Energy cut:\t" 
		  << options[i+1] << Sys::endl;
     
    }
    if(flag.compare("--strip-number") == 0){
      flag.clear();
      beam.strip_number = atof(options[i+1]);
      Sys::SysMsg << "Number of strips set: "
		  << beam.strip_number 
		  << Sys::endl;
    }
    if(flag.compare("--asymmetry") == 0){
      flag.clear();
      fAsymmetryAnalysis = true;
    }
    if(flag.compare("--resolution-analysis") == 0){
      flag.clear();
      fResolutionAnalysis = true;
    }
    if(flag.compare("--vftdc-analysis") == 0){
      flag.clear();
      fVetrocAnalysis = true;
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
    if(flag.compare("--halo-scale-x") == 0){
      flag.clear();
      compton.halo_scale_x = atof(options[i+1]);
    }
    if(flag.compare("--halo-scale-y") == 0){
      flag.clear();
      compton.halo_scale_y = atof(options[i+1]);
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
      if(beam.beam_energy == 10) {
	beam.sigma_ex = 434e-6;
	beam.sigma_ey = 199e-6;
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
    Sys::SysError << "In order to access asymetry analysis, two root files must be provided. Please use --help for more info." << Sys::endl;;
    exit(1);
  }

  double cs_right = CalculateIntegratedCS(1);
  double cs_left  = CalculateIntegratedCS(-1);

  Sys::SysCout << "Cross sections:\n" 
	    << "Left:\t" << cs_left
	    << "\nRight:\t" << cs_right << Sys::endl;
	     

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
  TH1D *left = new TH1D("left", "left", beam.strip_number, 1, beam.strip_number);
  TH1D *right = new TH1D("right", "right", beam.strip_number, 1, beam.strip_number);

  canvas->cd();

  Sys::SysCout << "Histograming data." << Sys::endl;

  std::vector <int> *id_l = 0;
  std::vector <int> *id_r = 0;

  std::vector <double> *totEdep_l = 0;
  std::vector <double> *totEdep_r = 0;

  int ID = 0;

  tree_left->ResetBranchAddresses();
  tree_left->SetBranchStatus("*", 0);
  tree_left->SetBranchStatus("id", 1);
  tree_left->SetBranchStatus("totEdep", 1);

  tree_left->SetBranchAddress("id", &id_l);
  tree_left->SetBranchAddress("totEdep", &totEdep_l);

  tree_right->ResetBranchAddresses();
  tree_right->SetBranchStatus("*", 0);
  tree_right->SetBranchStatus("id", 1);
  tree_right->SetBranchStatus("totEdep", 1);

  tree_right->SetBranchAddress("id", &id_r);
  tree_right->SetBranchAddress("totEdep", &totEdep_r);

  if(!(tree_left->GetLeaf("id"))){
    Sys::SysError << "Error finding leaf" << Sys::endl;
    exit(1);
  }

  if(!(tree_right->GetLeaf("id"))){
    Sys::SysError << "Error finding leaf" << Sys::endl;
    exit(1);
  }

  Sys::SysCout << (int)id_l->size() << Sys::endl;

  for(int i = 0; i < lentries; i++){
    tree_left->GetEntry(i);
    for(int j = 0; j < (int)id_l->size(); j++){
      if(id_l->at(j) < (beam.strip_number+1)){
   	ID = id_l->at(j);
	if(totEdep_l->at(j) > fEnergyCut)
	  left->Fill(ID);
      }
    }
  }
  for(int i = 0; i < rentries; i++){
    tree_right->GetEntry(i);
    for(int j = 0; j < (int)id_r->size(); j++){
      if(id_r->at(j) < (beam.strip_number+1)){
   	ID = id_r->at(j);
	if(totEdep_r->at(j) > fEnergyCut)
	  right->Fill(ID);
      }
    }
  }

  TF1 *fn_l = new TF1("fn_1", "[0]", 1, beam.strip_number+1);
  fn_l->SetParameter(0,1);

  TF1 *fn_r = new TF1("fn_r", "[0]", 1, beam.strip_number+1);
  fn_r->SetParameter(0,1);

  TH1D *raw_yield = (TH1D *)left->Clone("left");

  raw_yield->Add(right, 1);
  
  left->Multiply(fn_l, cs_left);
  right->Multiply(fn_r, cs_right);

  TH1D *asym = (TH1D *)left->Clone("asym");
  TH1D *yield  = (TH1D *)left->Clone("yield");

  asym->Add(right, -1);
  yield->Add(right, 1);

  asym->Divide(yield);
  asym->SetStats(0);

  std::fstream asym_out;
  std::fstream yield_out;

  asym_out.open("asymmetry.dat", std::fstream::out);
  yield_out.open("yield.dat", std::fstream::out);

  if(!yield_out.is_open()){ 
    Sys::SysError << "Failure." << Sys::endl;
    exit(1);
  }

  Sys::SysCout << "Strip\t|\tAsymmetry\t|\tBin Error\t|\tCalc Error\t|\tnL\t|\tnR" << Sys::endl;
  
  for(Int_t i = 1; i < (beam.strip_number+1); i++){
    
    if((raw_yield->GetBinContent(i)) == 0){
      raw_yield->SetBinError(i, 0.0);
      asym->SetBinError(i, 0.0);
    }
    else{
      raw_yield->SetBinError(i, 1/TMath::Sqrt(raw_yield->GetBinContent(i)));
      asym->SetBinError(i, 1/TMath::Sqrt(raw_yield->GetBinContent(i)));
    }
    
    asym_out << i << " " 
   	     << asym->GetBinContent(i) << " "
   	     << asym->GetBinError(i) << "\n";
    
    
    std::cout << i << " " 
   	      << asym->GetBinContent(i) << " "
       	      << asym->GetBinError(i) << std::endl;
    yield_out << i << " " 
   	      << raw_yield->GetBinContent(i) << " "
   	      << raw_yield->GetBinError(i) << "\n";
	      }
  asym_out.close();
  yield_out.close();
  
  Sys::SysMsg << "Histograming asymmetry." << Sys::endl;

  asym->Draw(); 
  asym->SetLineColor(9);
  asym->SetLineWidth(2);
  asym->SetTitle("Compton Detector Asymmetry");

  canvas->SaveAs("output/asymmetry_nowindow.png");
  canvas->SaveAs("output/asymmetry_nowindow.C");
  canvas->SaveAs("output/asymmetry_nowindow.root");
  delete canvas;
  
  Sys::SysMsg << "Analysis done." << Sys::endl;
  Sys::SysMsg << left->GetBinContent(2) << Sys::endl;

  if(fResolutionAnalysis){
    ScaleAsymmetry(left, right, raw_yield, 2);
    ScaleAsymmetry(left, right, raw_yield, 5);
    ScaleAsymmetry(left, right, raw_yield, 10);
    ScaleAsymmetry(left, right, raw_yield, 12);
    ScaleAsymmetry(left, right, raw_yield, 20);
  }

  exit(0);
}

void ComptonSimAnalysis::vfTDCAnalysis()
{

  Sys::SysMsg << "Processing vfTDC rootfile." << Sys::endl;

  if(!(fFileSet)){
    Sys::SysError << __FUNCTION__ << " In order to access vdTDC analysis, please provide vfTDC rootfile. Please use --help for more info." << Sys::endl;;
    exit(1);
  }

  TFile *file_vetroc = new TFile(fFileLocation.c_str());
  TTree *tree_vetroc = (TTree *)file_vetroc->Get("T");

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 800);
  canvas->SetTitle("Vetroc hits");

  int entries = tree_vetroc->GetEntries();
  Sys::SysMsg << "Events:." << entries << Sys::endl;

  if( entries == 0 ){
    Sys::SysError << __FUNCTION__  << " No events found. Exiting." << Sys::endl;
    exit(1);
  }
  TH1D *nhits_hist = new TH1D("nhits_hist", "nhits_hist", 192, 1, 192);
  TH1D *firsthit_hist = new TH1D("firsthit_hist", "firsthit_hist", 500, -500.0, 4500.0);

  canvas->cd();
  canvas->Divide(1,2);

  Sys::SysCout << "Histograming data." << Sys::endl;

  Double_t nhits[192];
  Double_t firsthit[192];

  int ID = 0;
  double TIME = 0;
  double max_time = 0;
  double min_time = 0;

  Sys::SysMsg << "Accessing rootfile." << Sys::endl;

  tree_vetroc->ResetBranchAddresses();
  tree_vetroc->SetBranchStatus("*", 0);
  tree_vetroc->SetBranchStatus("nhit", 1);
  tree_vetroc->SetBranchStatus("FirstHit", 1);

  tree_vetroc->SetBranchAddress("nhit", &nhits);
  tree_vetroc->SetBranchAddress("FirstHit", &firsthit);

  for(int i = 0; i < entries; i++){
    tree_vetroc->GetEntry(i);
    for(int j = (lvl_one_accept + 1); j < 192; j++){
      if( (nhits[j] > 0) && j == 148){
      // if( (nhits[j] > 0)){
   	ID = j;
	TIME = 0.001*(firsthit[j]-firsthit[lvl_one_accept]);
	if(TIME > max_time) max_time = TIME;
	if(TIME < min_time) min_time = TIME;

	if(TIME > 0){
	  nhits_hist->Fill(ID);
	  firsthit_hist->Fill(TIME);
	}
      }
    }
  }
  
  Sys::SysCout << "Histograming vfTDC hits." << Sys::endl;

  canvas->cd(1);
  nhits_hist->Draw(); 
  nhits_hist->SetLineColor(9);
  nhits_hist->SetLineWidth(2);
  nhits_hist->SetTitle("vfTDC Hits per Channel");

  canvas->cd(2);
  firsthit_hist->Draw(); 
  firsthit_hist->SetLineColor(9);
  firsthit_hist->SetLineWidth(2);
  firsthit_hist->SetTitle("vfTDC Hit Timing");
  firsthit_hist->GetXaxis()->SetRangeUser(min_time*0.9, max_time*1.1);

  canvas->SaveAs("output/vfTDC.png");
  canvas->SaveAs("output/vfTDC.C");

  delete canvas;
  
  Sys::SysMsg << "Analysis done." << Sys::endl;

  exit(0);
}

void ComptonSimAnalysis::ScaleAsymmetry(TH1D *l, TH1D* r, TH1D *y_raw, int multiplier = 1)
{
  if(!(fFileLeftSet && fFileRightSet)){
    Sys:: SysError << __FUNCTION__ << " In order to access asymetry analysis, two root files must be provided. Please use --help for more info." << Sys::endl;
    exit(1);
  }

  // Make a copy of the left and right histograms so that they are not edited directly

  TH1D *left  = (TH1D *)l->Clone("asym");
  TH1D *right  = (TH1D *)r->Clone("asym");
  TH1D *raw_yield = (TH1D *)y_raw->Clone("raw_yield");

  const int strips = (beam.strip_number)/multiplier;

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 800);
  canvas->SetTitle("Compton Detector Scaled Asymmetry");

  left->Rebin(multiplier);
  right->Rebin(multiplier);
  raw_yield->Rebin(multiplier);

  TH1D *asym = (TH1D *)left->Clone("asym");
  TH1D *yield  = (TH1D *)left->Clone("yield");

  asym->Add(right, -1);
  yield->Add(right, 1);

  asym->Divide(yield);
  asym->SetStats(0);

  std::fstream asym_out;
  std::fstream yield_out;

  asym_out.open(Form("asymmetry_m%d.dat", multiplier), std::fstream::out);
  yield_out.open(Form("yield_m%d.dat", multiplier), std::fstream::out);

  if(!yield_out.is_open()){ 
    std::cout << "Failure." << std::endl;
    exit(1);
  }

  Sys::SysCout << "Strip\t|\tAsymmetry\t|\tBin Error\t|\tCalc Error\t|\tnL\t|\tnR" << Sys::endl;

  for(Int_t i = 1; i < (strips+1); i++){


    if((raw_yield->GetBinContent(i)) == 0){
      raw_yield->SetBinError(i, 0.0);
      asym->SetBinError(i, 0.0);
    }
    else{
      raw_yield->SetBinError(i, 1/TMath::Sqrt(raw_yield->GetBinContent(i)));
      asym->SetBinError(i, 1/TMath::Sqrt(raw_yield->GetBinContent(i)));
    }

    asym_out << i << " " 
	     << asym->GetBinContent(i) << " " 
	     << asym->GetBinError(i) << "\n";

    std::cout << i << " " 
	      << raw_yield->GetBinContent(i) << " "
     	      << raw_yield->GetBinError(i) << " " 
     	      << asym->GetBinContent(i) << " "
      	      << asym->GetBinError(i) << std::endl;
    
    yield_out << i << " " 
	      << raw_yield->GetBinContent(i) << " "
	      << raw_yield->GetBinError(i) << "\n";
  }
  asym_out.close();
  yield_out.close();

  Sys::SysCout << "Histograming asymmetry." << Sys::endl;

  asym->Draw(); 
  asym->SetLineColor(9);
  asym->SetLineWidth(2);
  asym->SetTitle("Compton Detector Asymmetry");

  Sys::SysMsg << "Analysis scaling done." << Sys::endl;

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
    Sys::SysError << __FUNCTION__ << " Aperture size not defined. Using default aperture size (+- 0.5 cm). Please refer to --help for instruction on setting this value." << Sys::endl;
    return 0;
  }
  // the factor of ten is a multplier for the halo and the '1e2' is to convert to cm
  compton.gaussian_weight = boost::math::erfc(compton.aperture_size/(std::sqrt(2)*10e2*beam.sigma_ey));
  Sys::SysMsg << "\n>>>>> Gaussian weight: " << compton.gaussian_weight << Sys::endl;
  std::cout << compton.aperture_size << " " << beam.sigma_ey << std::endl;

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
  Generator::Initialize(beam.beam_energy, beam.laser_energy, beam.polarization); 

  compton.cross_section = 2*TMath::Pi()*this->GetFunction((char *)"cs")->Integral(0,1);

  return(2*TMath::Pi()*this->GetFunction((char *)"cs")->Integral(0,1));
}

double ComptonSimAnalysis::CalculateIntegratedCS(double polarization)
{
  Generator::Initialize(beam.beam_energy, beam.laser_energy, polarization); 

  return(2*TMath::Pi()*this->GetFunction((char *)"cs")->Integral(0,1));
}

#endif

