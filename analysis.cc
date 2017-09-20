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
#include "FontColor.hh"
#include "PhysicalConstants.hh"

int main(int argc, char *argv[])
{

  ComptonSimAnalysis *simulation = new ComptonSimAnalysis();

  simulation->GetOptions(argv);
  if(simulation->fGraphicsShow) simulation->InitGraphicsEngine(argc, argv);

  simulation->GenerateAsymmetry((char *)""); // The char * casting removes a deprecatred warning caused by difference between char * in C and C++

  simulation->CalculateIntegratedCS();
  std::cout << "Integrated cross section: " << simulation->compton.cross_section << std::endl;
  simulation->CalculateLuminosity(1);
  simulation->CalculateGaussianWeight();

  std::vector <int> *id    = 0;
  std::vector <int> *pid   = 0;
  std::vector <int> *otid  = 0;
  std::vector <int> *mtid  = 0;
  std::vector <int> *hitn  = 0;
  std::vector <int> *totEdep  = 0;
  
  std::vector <int> BinContent(200);
  std::fill(BinContent.begin(), BinContent.end(), 0);  

  if(!(simulation->fFileSet)){
    std::cout << "Must define rootfile first. Exiting." << std::endl;
    exit(1);
  }  

  TFile *file = new TFile(TString(simulation->fFileLocation));
  TTree *fChain = (TTree *)file->Get("flux");

  int entries = fChain->GetEntries();
  
  if(entries == 0){
    simulation->PrintError((const char*)"No events found. Exiting.");
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
  
  fChain->SetBranchAddress("id", &id);
  fChain->SetBranchAddress("pid", &pid);
  fChain->SetBranchAddress("hitn", &hitn);
  fChain->SetBranchAddress("otid", &otid);
  fChain->SetBranchAddress("mtid", &mtid);
  fChain->SetBranchAddress("totEdep", &totEdep);
  
  if(!(fChain->GetLeaf("id"))){
    std::cout << "Error finding leaf" << std::endl;
    exit(1);
  }
    
  Int_t ID = 0;
  
  for(Int_t i = 0; i < entries; i++){
    fChain->GetEntry(i);
    for(Int_t j = 0; j < (Int_t)id->size(); j++){
      if(id->at(j) < 201){
	ID = id->at(j);
	BinContent[ID]++;
      }
    }
  }
  double w = 1;
  if(simulation->fBackgroundWeight) w = (pchep::coulomb*simulation->pressure_conversion)/entries; 
  if(simulation->fComptonWeight)    w = (simulation->compton.luminosity*simulation->compton.cross_section)/entries; // units Hz corrected for cross_section/luminosity weighted
  if(simulation->fHaloWeight){
    simulation->CalculateHaloFraction();
    w = (simulation->compton.halo_ratio*pchep::coulomb*simulation->compton.gaussian_weight)/entries; // units Hz corrected for luminosity weighted for halo
  }

  TH1D *hist = new TH1D("hist", "hist", 200, 0, 200);

  for(Int_t i = 0; i < 201; i++){
    int entry_count = BinContent[i]; 
    std::cout << "Strip: " << i << " entry# " << entry_count << std::endl;
    for(Int_t j = 0; j < entry_count; j++){
      hist->Fill(i, w);  
    }
  }
  
    
  std::cout << "<<<<< Histogram has " << hist->GetEntries() << " entries.\n" << std::endl;
  std::cout << "<<<<< Integral:" << hist->Integral() << std::endl;
  
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

  canvas->SaveAs(TString(simulation->fFileOutput));
  
  if(simulation->fGraphicsShow) simulation->RunGraphicsEngine();
  
  return 0;
 }

