#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLT/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

bool    matchMuon   (HLTMuonCand, std::vector<HLTObjCand>, std::string);
void    dorate      (TFile*, TFile*);

void readNtuplesRate(){

  std::string thestr = "_XXX";

  TFile* inputfile0 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD20to30.root","READ");
  TFile* outfile0 = TFile::Open(Form("rates/rate_QCD20to30%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile0, outfile0);
  
  TFile* inputfile1 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD30to50.root","READ");
  TFile* outfile1 = TFile::Open(Form("rates/rate_QCD30to50%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile1, outfile1);

  TFile* inputfile2 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD50to80.root","READ");
  TFile* outfile2 = TFile::Open(Form("rates/rate_QCD50to80%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile2, outfile2);

  TFile* inputfile3 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD80to120.root","READ");
  TFile* outfile3 = TFile::Open(Form("rates/rate_QCD80to120%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile3, outfile3);

  TFile* inputfile4 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_WToMuNu.root","READ");
  TFile* outfile4 = TFile::Open(Form("rates/rate_WMuNu%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile4, outfile4);

  TFile* inputfile5 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_DYToLL.root","READ");
  TFile* outfile5 = TFile::Open(Form("rates/rate_DYLL%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile5, outfile5);
  
  return;

}

void dorate( TFile* inputfile,  TFile* outfile){
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree");
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }

  MuonEvent* ev = new MuonEvent();
  tree -> SetBranchAddress("event", &ev);

  TH1F* muon_size           = new TH1F("muon_size"           ,"muon_size"            ,   10,  0, 10 );
  TH1F* wp_passing          = new TH1F("wp_passing"          ,"wp_passing"           ,   60,  0,  60);

  int nentries = tree->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  bool passwp = false;
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {

    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    passwp          = false; 

//     if (!ev-> hlt.find("HLT_Mu20_v2")) continue;  // to evaluate rate wrt to reference trigger

    unsigned int nhltmuons = ev->hltmuons.size();
    muon_size       -> Fill( nhltmuons );

    for (int imu = 0; imu < nhltmuons; imu++){
      
      HLTMuonCand theL3 = ev -> hltmuons.at(imu);
      if ( matchMuon(ev ->hltmuons.at(imu), ev -> hlt.objects, "FILTER_TO_MATCH::REHLT") && passwp==false) {
        wp_passing -> Fill(3);
	    passwp = true;
	  }
    }
  }

  outfile            -> cd();
  
  muon_size          -> Write();
  wp_passing         -> Write();

  outfile            -> Close();  
  
  return;
}



bool matchMuon(HLTMuonCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1;
  float theDR = 100;
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) {
    if ( it->filterTag.compare(tagFilterName) == 0) {
      
      theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
      if (theDR < minDR){
        minDR = theDR;
        match = true;
      }
    }
  }
  
  return match;
}




