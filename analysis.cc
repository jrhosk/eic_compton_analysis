#include <iostream>
#include <string>
#include <stdio.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TSystem.h"
#include "TChainElement.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <vector>

// Simulation Analysis class
#include "ComptonSimAnalysis.hh"

// Generator
#include "Generator.hh"

// Addition QoL libraries
#include "PhysicalConstants.hh"
#include "MsgStream.hh"
#include "SysMsg.hh"

int main(int argc, char *argv[])
{

  ComptonSimAnalysis *simulation = new ComptonSimAnalysis();

  simulation->GetOptions(argv);

  if(simulation->fGraphicsShow) simulation->InitGraphicsEngine(argc, argv);
  if(simulation->fAsymmetryAnalysis) simulation->AsymmetryAnalysis();
  if(simulation->fVetrocAnalysis) simulation->vfTDCAnalysis();

  simulation->GenerateAsymmetry((char *)""); // The char * casting removes a deprecatred warning caused by difference between char * in C and C++

  simulation->CalculateIntegratedCS();
  Sys::SysCout << "Integrated cross section: " << simulation->compton.cross_section << Sys::endl;
  simulation->CalculateLuminosity(1);
  simulation->CalculateGaussianWeight();
  
  std::vector <int> *id    = 0;
  std::vector <int> *pid   = 0;
  std::vector <int> *otid  = 0;
  std::vector <int> *mtid  = 0;
  std::vector <int> *hitn  = 0;
  std::vector <int> *procID  = 0;

  std::vector <double> *totEdep  = 0;
  
  std::vector <int> BinContent(simulation->beam.strip_number);
  std::fill(BinContent.begin(), BinContent.end(), 0);  

  if(!(simulation->fFileSet)){
    Sys::SysError << "Must define rootfile first. Exiting." << Sys::endl;
    exit(1);
  }  

  TFile *file = new TFile(TString(simulation->fFileLocation));
  TTree *fChain = (TTree *)file->Get("flux");

  int entries = fChain->GetEntries();
  
  if(entries == 0){
    Sys::SysError << __FUNCTION__ << " No events found. Exiting" << Sys::endl;
    exit(1);
  }
  
  std::cout << "Number of events: " << entries << std::endl;
  
  fChain->ResetBranchAddresses();
  fChain->SetBranchStatus("*", 0);
  fChain->SetBranchStatus("id", 1);
  fChain->SetBranchStatus("pid", 1);
  fChain->SetBranchStatus("hitn", 1);
  fChain->SetBranchStatus("otid", 1);
  fChain->SetBranchStatus("mtid", 1);
  fChain->SetBranchStatus("totEdep", 1);
  fChain->SetBranchStatus("procID", 1);
  
  fChain->SetBranchAddress("id", &id);
  fChain->SetBranchAddress("pid", &pid);
  fChain->SetBranchAddress("hitn", &hitn);
  fChain->SetBranchAddress("otid", &otid);
  fChain->SetBranchAddress("mtid", &mtid);
  fChain->SetBranchAddress("totEdep", &totEdep);
  fChain->SetBranchAddress("procID", &procID);
  
  if(!(fChain->GetLeaf("id"))){
    Sys::SysError << "Error finding leaf" << Sys::endl;
    exit(1);
  }
    
  int ID = 0;
  
  for(int i = 0; i < entries; i++){
    fChain->GetEntry(i);
    for(int j = 0; j < (int)id->size(); j++){
      if(id->at(j) < (simulation->beam.strip_number)+1){

	if(totEdep->at(j) > simulation->fEnergyCut){
	  ID = id->at(j);
	  BinContent[ID]++;
	}

      }
    } 
  }
  
  double w = 1;
  if(simulation->fBackgroundWeight) w = (pchep::coulomb*simulation->pressure_conversion)/entries; 
  if(simulation->fComptonWeight)    w = (simulation->compton.luminosity*simulation->compton.cross_section)/entries; // units Hz corrected for cross_section/luminosity weighted
  if(simulation->fHaloWeight){
    simulation->CalculateHaloFraction();
    w = (simulation->compton.halo_ratio*pchep::coulomb*simulation->compton.gaussian_weight)/entries; // weighted for halo
    std::cout << "Halo: ratio:" << simulation->compton.halo_ratio << std::endl;
    std::cout << "\nGaussian weight:  " << simulation->compton.gaussian_weight << std::endl;

  }

  TH1D *hist = new TH1D("hist", "hist", simulation->beam.strip_number, 1, simulation->beam.strip_number);

  for(int i = 1; i < (simulation->beam.strip_number)+1; i++){
    int entry_count = BinContent[i]; 

    for(int j = 0; j < entry_count; j++){
      hist->Fill(i, w);  
    }
  }
  
  int max_bin =  hist->GetMaximumBin();
  double max_bin_content = hist->GetBinContent(max_bin);

  std::cout << "<<<<< Histogram maximum: " << hist->GetMaximumBin() << " " << hist->GetBinContent(max_bin) << Sys::endl;
  std::cout << "<<<<< Histogram has " << hist->GetEntries() << " entries.\n" << Sys::endl;
  std::cout << "<<<<< Integral:" << hist->Integral() << Sys::endl;
  
  std::fstream halo_out;

  halo_out.open("halo_out.dat", std::fstream::out | std::fstream::app);

  halo_out << simulation->compton.halo_scale_x << " "
	   << simulation->compton.halo_scale_y << " "
	   << (1e-3)*hist->Integral() << std::endl;

  TCanvas *canvas = new TCanvas("canvas","canvas",1200, 600);
  canvas->cd();
  
  hist->SetStats(0);
  hist->Draw("");
  hist->SetLineColor(4);
  hist->SetFillColor(9);
  hist->SetLineWidth(2);  
  hist->GetXaxis()->SetTitle("Strip");
  hist->GetYaxis()->SetTitle("Rate(Hz)");
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->SetTitle("Detector Rate");

  std::cout << "Entries:\t" << hist->GetEntries() << Sys::endl;

  canvas->SaveAs(Form("output/%s.C", (simulation->fFileOutput).c_str()));
  canvas->SaveAs(Form("output/%s.root", (simulation->fFileOutput).c_str()));
  
  if(simulation->fGraphicsShow) simulation->RunGraphicsEngine();
  
  return 0;
 }

