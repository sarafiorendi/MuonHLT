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

bool selectTagMuon  (MuonCand, TH1F* );
bool selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);


void readNtuples(){

  TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/redo/CMSSW_7_4_8/src/HLTrigger/Configuration/test/ntuples/muonNtuple_onZMuMuSkim.root","READ");
//   TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/redo/CMSSW_7_4_8/src/HLTrigger/Configuration/test/ntuples/muonNtupleMC_DYJetsToLL.root","READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open("efficiency_onZMuMuSkim_hltisoTag_offlineIsoTag_match.root","RECREATE");
//   TFile* outfile = TFile::Open("efficiency_onDYJetsToLL_MC.root","RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TTree *tree = (TTree*) inputfile -> Get("ntupleproducer/muonTree");
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }

  TH1F* dimuon_mass         = new TH1F("dimuon_mass"         ,"dimuon_mass"          , 1500,  0, 150);
  TH1F* tagiso              = new TH1F("tagiso"              ,"tagiso"               ,  100,  0, 1  );
  TH1F* tagMuonPt           = new TH1F("tagMuonPt"           ,"tagMuonPt"            ,  150,  0, 150);

  TH1F* muonPt_den          = new TH1F("muonPt_den"          ,"muonPt_den"           ,  150,  0, 150);
  TH1F* muonPt_num          = new TH1F("muonPt_num"          ,"muonPt_num"           ,  150,  0, 150);
  TH1F* muonPt_barrel_den   = new TH1F("muonPt_barrel_den"   ,"muonPt_barrel_den"    ,  150,  0, 150);
  TH1F* muonPt_barrel_num   = new TH1F("muonPt_barrel_num"   ,"muonPt_barrel_num"    ,  150,  0, 150);
  TH1F* muonPt_endcap_den   = new TH1F("muonPt_endcap_den"   ,"muonPt_endcap_den"    ,  150,  0, 150);
  TH1F* muonPt_endcap_num   = new TH1F("muonPt_endcap_num"   ,"muonPt_endcap_num"    ,  150,  0, 150);

  TH1F* muonEta_den         = new TH1F("muonEta_den"         ,"muonEta_den"          ,  600, -3,   3);
  TH1F* muonEta_num         = new TH1F("muonEta_num"         ,"muonEta_num"          ,  600, -3,   3);

  TH1F* muonPhi_den         = new TH1F("muonPhi_den"         ,"muonPhi_den"          ,   40, -3.2,   3.2);
  TH1F* muonPhi_num         = new TH1F("muonPhi_num"         ,"muonPhi_num"          ,   40, -3.2,   3.2);

  TH1F* nvtx_num            = new TH1F("nvtx_num"            ,"nvtx_num"             ,   60,  0,  60);
  TH1F* nvtx_den            = new TH1F("nvtx_den"            ,"nvtx_den"             ,   60,  0,  60);
  TH1F* nvtx_barrel_num     = new TH1F("nvtx_barrel_num"     ,"nvtx_barrel_num"      ,   60,  0,  60);
  TH1F* nvtx_barrel_den     = new TH1F("nvtx_barrel_den"     ,"nvtx_barrel_den"      ,   60,  0,  60);
  TH1F* nvtx_endcap_num     = new TH1F("nvtx_endcap_num"     ,"nvtx_endcap_num"      ,   60,  0,  60);
  TH1F* nvtx_endcap_den     = new TH1F("nvtx_endcap_den"     ,"nvtx_endcap_den"      ,   60,  0,  60);

  TH1F* muonIso_den         = new TH1F("muonIso_den"         ,"muonIso_den"          ,  100,  0,   1);
  TH1F* muonIso_num         = new TH1F("muonIso_num"         ,"muonIso_num"          ,  100,  0,   1);
  TH1F* muonIso_barrel_den  = new TH1F("muonIso_barrel_den"  ,"muonIso_barrel_den"   ,  100,  0,   1);
  TH1F* muonIso_barrel_num  = new TH1F("muonIso_barrel_num"  ,"muonIso_barrel_num"   ,  100,  0,   1);
  TH1F* muonIso_endcap_den  = new TH1F("muonIso_endcap_den"  ,"muonIso_endcap_den"   ,  100,  0,   1);
  TH1F* muonIso_endcap_num  = new TH1F("muonIso_endcap_num"  ,"muonIso_endcap_num"   ,  100,  0,   1);
  TH1F* muonIso04_den       = new TH1F("muonIso04_den"       ,"muonIso04_den"        ,  100,  0,   1);
  TH1F* muonIso04_num       = new TH1F("muonIso04_num"       ,"muonIso04_num"        ,  100,  0,   1);
  
  double offlineiso   = 100;
  double offlineiso04 = 100;

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
//     if (nhltmuons > 0) std::cout << "Number of hlt muons = " << nhltmuons << std::endl;
    
    if (!ev-> hlt.find("HLT_IsoMu20_v2")) continue;
    

    for (int imu = 0; imu < nmuons; imu++){
      
      // select the tag muon        
      if (! selectTagMuon(ev -> muons.at(imu), tagiso)) continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hlt.objects, "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09::REHLT")) continue;
      tagMuonPt -> Fill(ev -> muons.at(imu).pt);
      
      for (int jmu = 0; jmu < nmuons; jmu++){
        // select the probe muon
        if (! selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;

        offlineiso = ev -> muons.at(jmu).chargedDep_dR03 + std::max(0.,
                     ev -> muons.at(jmu).photonDep_dR03 + ev -> muons.at(jmu).neutralDep_dR03 - 0.5*ev -> muons.at(jmu).puPt_dR03);
        offlineiso = offlineiso / ev -> muons.at(jmu).pt;

        offlineiso04 = ev -> muons.at(jmu).chargedDep_dR04 + std::max(0.,
                       ev -> muons.at(jmu).photonDep_dR04 + ev -> muons.at(jmu).neutralDep_dR04 - 0.5*ev -> muons.at(jmu).puPt_dR04);
        offlineiso04 = offlineiso04 / ev -> muons.at(jmu).pt;
        
//         if (!(offlineiso04 < 0.12)) continue;
        if (!(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, "hltL3fL1sMu16L1f0L2f10QL3Filtered20Q::REHLT"))) continue;

        muonPt_den    -> Fill( ev -> muons.at(jmu).pt );
        muonEta_den   -> Fill( ev -> muons.at(jmu).eta);
        muonPhi_den   -> Fill( ev -> muons.at(jmu).phi);
        muonIso_den   -> Fill( offlineiso   );
        muonIso04_den -> Fill( offlineiso04 );
        nvtx_den      -> Fill( ev -> nVtx   );
        if (fabs(ev -> muons.at(jmu).eta) < 1.479){
          muonPt_barrel_den  -> Fill( ev -> muons.at(jmu).pt );
          nvtx_barrel_den    -> Fill( ev -> nVtx );
          muonIso_barrel_den -> Fill( offlineiso );
        }
        else{
          muonPt_endcap_den  -> Fill( ev -> muons.at(jmu).pt );
          nvtx_endcap_den    -> Fill( ev -> nVtx );
          muonIso_endcap_den -> Fill( offlineiso );
        }
        if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09::REHLT")){
          muonPt_num    -> Fill( ev -> muons.at(jmu).pt);
          muonEta_num   -> Fill( ev -> muons.at(jmu).eta);
          muonPhi_num   -> Fill( ev -> muons.at(jmu).phi);
          muonIso_num   -> Fill( offlineiso   );
          muonIso04_num -> Fill( offlineiso04 );
          nvtx_num      -> Fill( ev -> nVtx   );
          if (fabs(ev -> muons.at(jmu).eta) < 1.479){
            muonPt_barrel_num  -> Fill( ev -> muons.at(jmu).pt );
            nvtx_barrel_num    -> Fill( ev -> nVtx );
            muonIso_barrel_num -> Fill( offlineiso );
          }
          else{
            muonPt_endcap_num  -> Fill( ev -> muons.at(jmu).pt );
            nvtx_endcap_num    -> Fill( ev -> nVtx );
            muonIso_endcap_num -> Fill( offlineiso );
          }
        }
      }
      
    }

// HLT_IsoMu20_eta2p1_v1



  }
  



  outfile            -> cd();
  muonPt_den         -> Write();
  muonPt_num         -> Write();
  muonPt_barrel_den  -> Write();
  muonPt_barrel_num  -> Write();
  muonPt_endcap_den  -> Write();
  muonPt_endcap_num  -> Write();
  tagMuonPt          -> Write();
  
  muonEta_den        -> Write();
  muonEta_num        -> Write();
  muonPhi_den        -> Write();
  muonPhi_num        -> Write();
  
  nvtx_num           -> Write();
  nvtx_den           -> Write();
  nvtx_barrel_num    -> Write();
  nvtx_barrel_den    -> Write();
  nvtx_endcap_num    -> Write();
  nvtx_endcap_den    -> Write();
  
  muonIso_num        -> Write();
  muonIso_den        -> Write();
  muonIso_barrel_num -> Write();
  muonIso_barrel_den -> Write();
  muonIso_endcap_num -> Write();
  muonIso_endcap_den -> Write();
  muonIso04_num      -> Write();
  muonIso04_den      -> Write();

  dimuon_mass        -> Write();
  tagiso             -> Write();
  
  muonPt_num         -> Draw();
  outfile            -> Close();  
  
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


bool selectTagMuon(MuonCand mu, TH1F* tagiso){
  
  if (!( mu.pt         > 20  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0.,
                       mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  tagiso -> Fill(offlineiso04);
  if (!(offlineiso04   < 0.12)) return false;
  
  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt         > 20  )) return false; 
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



