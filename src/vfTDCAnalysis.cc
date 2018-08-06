// C++ Libraries

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>

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

#include "vfTDCAnalysis.hh"
#include "ComptonSimAnalysis.hh"
#include "MsgStream.hh"
#include "SysMsg.hh"
#include "PhysicalConstants.hh"

vfTDCAnalysis::vfTDCAnalysis(){
  fFileSet = false;
}

vfTDCAnalysis::vfTDCAnalysis(std::string path){
  fFileSet = true;
  fFileLocation = path;
}

int vfTDCAnalysis::vfTDCAnalyze()
{

  std::cout << "Processing vfTDC rootfile." << std::endl;

  if(!(fFileSet)){
    Sys::SysError << __FUNCTION__ << " In order to access vdTDC analysis, please provide vfTDC rootfile. Please use --help for more info." << Sys::endl;;
    exit(1);
  }

  TFile *file_vetroc = new TFile(fFileLocation.c_str());
  TTree *tree_vetroc = (TTree *)file_vetroc->Get("T");

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 800);
  canvas->SetTitle("Vetroc hits");

  int entries = tree_vetroc->GetEntries();
  std::cout << "Events:." << entries << std::endl;

  if( entries == 0 ){
    Sys::SysError << __FUNCTION__  << " No events found. Exiting." << Sys::endl;
    exit(1);
  }
  TH1D *nhits_hist = new TH1D("nhits_hist", "nhits_hist", 192, 1, 192);
  TH1D *firsthit_hist = new TH1D("firsthit_hist", "firsthit_hist", 500, -500.0, 4500.0);

  canvas->cd();
  canvas->Divide(1,2);

  Sys::SysCout << "Histograming data." << Sys::endl;

  int ID = 0;
  int channel_key = 0;

  double nhits[192];
  double firsthit[192];
  double TIME = 0;
  double max_time = 0;
  double min_time = 0;

  std::vector <TH1D *> hist_chan;
  std::vector <TH1D *> hist_time;

  std::cout << "Accessing rootfile." << std::endl;

  tree_vetroc->ResetBranchAddresses();
  tree_vetroc->SetBranchStatus("*", 0);
  tree_vetroc->SetBranchStatus("nhit", 1);
  tree_vetroc->SetBranchStatus("FirstHit", 1);

  tree_vetroc->SetBranchAddress("nhit", &nhits);
  tree_vetroc->SetBranchAddress("FirstHit", &firsthit);

  for(int i = 0; i < entries; i++){
    tree_vetroc->GetEntry(i);
    for(int j = (lvl_one_accept + 1); j < 192; j++){
      if( (nhits[j] > 0)){
   	ID = j;
	TIME = 0.001*(firsthit[j]-firsthit[lvl_one_accept]);
	if(TIME > max_time) max_time = TIME;
	if(TIME < min_time) min_time = TIME;

	if(TIME > 0){
	  channel_key = GetChannelKey(j);
	  if( (channel_key) >= 0){
	    hist_chan[channel_key]->Fill(ID);
	    hist_time[channel_key]->Fill(TIME);
	  } else {
	    hist_chan.push_back(new TH1D(Form("hist_chan_%i", j), 
	     				 Form("hist_chan_%i", j), 
	     				 192, 1.0, 192.0) );
	    hist_time.push_back(new TH1D(Form("hist_time_%i", j), 
	     				 Form("hist_time_%i", j), 
	     				 500, -500.0, 4500.0) );
	    hist_chan.back()->Fill(ID);
	    hist_time.back()->Fill(TIME);
	  }

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

  for(int i = 0; i < (int)(hist_chan.size()); i++){

     canvas->cd(1);
     nhits_hist->Draw(); 
     nhits_hist->SetLineColor(9);
     nhits_hist->SetLineWidth(2);
     nhits_hist->SetTitle(Form("vfTDC Hits - Channel %i", i));

     hist_chan[i]->Draw("same");
     hist_chan[i]->SetLineColor(2);
     hist_chan[i]->SetLineWidth(2);
    
     canvas->cd(2);
     firsthit_hist->Draw(); 
     firsthit_hist->SetLineColor(9);
     firsthit_hist->SetLineWidth(2);
     firsthit_hist->SetTitle(Form("vfTDC Hit Timing - Channel %i", i));
     firsthit_hist->GetXaxis()->SetRangeUser(min_time*0.9, max_time*1.1);
     
     hist_time[i]->Draw("same");
     hist_time[i]->SetLineColor(2);
     hist_time[i]->SetLineWidth(2);
     
     canvas->SaveAs(Form("output/vfTDC_%i.png", channelList[i]));
     canvas->SaveAs(Form("output/vfTDC_%i.C", channelList[i]));
  }
  
  delete canvas;
  
  Sys::SysMsg << "Analysis done." << Sys::endl;

  exit(0);
}

int vfTDCAnalysis::GetChannelKey(int channel)
{
  for(int i = 0; i < (int)(channelList.size()); i++){
    if(channelList[i] == channel){
      return i;
    }
  }
  channelList.push_back(channel);
  return(-1);

}
