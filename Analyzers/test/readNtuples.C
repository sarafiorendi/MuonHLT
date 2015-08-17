#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
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

// const MuonEvent& GetEvent(){return *_event;};
double muonmass = 0.10565837;

bool selectTagMuon  (MuonCand );
bool selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);


void readNtuples(){

  TFile* inputfile = TFile::Open("ntuples/muonNtuples_2015B_MuonJsonv4_partial.root","READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open("provaOut.root","RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree");
  if (!tree) {
    std::cout << "tree not found " << std::endl;
    return;
  }

  TH1F* dimuon_mass         = new TH1F("dimuon_mass"         ,"dimuon_mass"          , 1500,  0, 150);
  TH1F* tagMuonPt           = new TH1F("tagMuonPt"           ,"tagMuonPt"            ,  150,  0, 150);
  TH1F* muonPt_den          = new TH1F("muonPt_den"          ,"muonPt_den"           ,  150,  0, 150);
  TH1F* muonPt_num          = new TH1F("muonPt_num"          ,"muonPt_num"           ,  150,  0, 150);
  TH1F* muonEta_den         = new TH1F("muonEta_den"         ,"muonEta_den"          ,  600, -3,   3);
  TH1F* muonEta_num         = new TH1F("muonEta_num"         ,"muonEta_num"          ,  600, -3,   3);

  MuonEvent* ev = new MuonEvent();
  TBranch* evBranch = tree->GetBranch("event");
  evBranch -> SetAddress(&ev);

  int nentries = tree->GetEntriesFast();
  std::cout << "Number of entries = " << nentries << std::endl;


  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
//     ClearEventVariables();
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
//     std::cout << " ---- analyzing event: " << ev -> eventNumber << std::endl;
    
    unsigned int nmuons = ev->muons.size();
    if (nmuons < 2) continue;

    unsigned int nhltmuons = ev->hltmuons.size();
    if (nhltmuons > 0) std::cout << "Number of hlt muons = " << nhltmuons << std::endl;
    
    if (!ev-> hlt.find("HLT_Mu20_v1")) continue;
    

    for (int imu = 0; imu < nmuons; imu++){
      
      // select the tag muon        
      if (! selectTagMuon(ev -> muons.at(imu))) continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hlt.objects, "hltL3fL1sMu16L1f0L2f10QL3Filtered20Q::HLT")) continue;
      tagMuonPt -> Fill(ev -> muons.at(imu).pt);
      
      for (int jmu = 0; jmu < nmuons; jmu++){
        // select the probe muon
        if (! selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
        muonPt_den  -> Fill( ev -> muons.at(jmu).pt );
        muonEta_den -> Fill( ev -> muons.at(jmu).eta);
        if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09::HLT")){
          muonPt_num  -> Fill( ev -> muons.at(jmu).pt);
          muonEta_num -> Fill( ev -> muons.at(jmu).eta);
        }
      }
      
    }

// HLT_IsoMu20_eta2p1_v1



  }
  



  outfile     -> cd();
  muonPt_den  -> Write();
  muonPt_num  -> Write();
  tagMuonPt   -> Write();
  
  muonEta_den -> Write();
  muonEta_num -> Write();
  
  dimuon_mass -> Write();
  
  muonPt_num  -> Draw();
  outfile    -> Close();  
  
  return;
}



bool matchMuon(MuonCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

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


bool selectTagMuon(MuonCand mu){
  
  if (!( mu.pt         > 25  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  
  
  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt         > 25  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  if (mu.charge * tagMu.charge > 0) return false;
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();
  dimuon_mass -> Fill(mumumass);
  if (! (mumumass > 86. && mumumass < 96. )) return false;
  
  return true;
}



