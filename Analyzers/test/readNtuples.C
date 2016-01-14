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

double muonmass = 0.10565837;

bool selectTagMuon  (MuonCand, TH1F* );
bool selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);

std::string L3filter_18   = "hltL3fL1sMu16L1f0L2f10QL3Filtered18Q::HLT";
std::string isofilter_18  = "hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09::HLT";
std::string isofilter_o18 = "hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3OldCaloIsotrkIsoFiltered0p09::HLT";
std::string L3filter, isofilter;

double offlineIsoCut = 0.15;

void readNtuples(){

  std::string hltname = "HLT_IsoMu18_v2";
  
  if (hltname.find("HLT_IsoMu18_v2") == 0 ){
     L3filter  = L3filter_18;
     isofilter = isofilter_18;
  }
  else if (hltname.find("HLT_OldIsoMu18_v1") == 0 ){
     L3filter  = L3filter_18;
     isofilter = isofilter_o18;
  }
  else{
    std::cout << "no filters corresponding to given HLT path" << std::endl;
    return;
  }


  TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/CMSSW_7_4_15/src/MuonHLT/Analyzers/test/muonNtuple_run258158.root","READ");
//   TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/redo/CMSSW_7_4_8/src/HLTrigger/Configuration/test/ntuples/muonNtupleMC_DYJetsToLL.root","READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open(Form("efficiency_fullpath_onZMuMuSkim_run258158_%s_offlineIso0p15.root", hltname.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree");
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
  TH1F* muonPt_eta0_den     = new TH1F("muonPt_eta0_den"     ,"muonPt_eta0_den"      ,  150,  0, 150);
  TH1F* muonPt_eta0_num     = new TH1F("muonPt_eta0_num"     ,"muonPt_eta0_num"      ,  150,  0, 150);
  TH1F* muonPt_eta1_den     = new TH1F("muonPt_eta1_den"     ,"muonPt_eta1_den"      ,  150,  0, 150);
  TH1F* muonPt_eta1_num     = new TH1F("muonPt_eta1_num"     ,"muonPt_eta1_num"      ,  150,  0, 150);
  TH1F* muonPt_eta2_den     = new TH1F("muonPt_eta2_den"     ,"muonPt_eta2_den"      ,  150,  0, 150);
  TH1F* muonPt_eta2_num     = new TH1F("muonPt_eta2_num"     ,"muonPt_eta2_num"      ,  150,  0, 150);
 
  TH1F* muonEta_den         = new TH1F("muonEta_den"         ,"muonEta_den"          ,  600, -3,   3);
  TH1F* muonEta_num         = new TH1F("muonEta_num"         ,"muonEta_num"          ,  600, -3,   3);

  TH1F* muonPhi_den         = new TH1F("muonPhi_den"         ,"muonPhi_den"          ,   40, -3.2,   3.2);
  TH1F* muonPhi_num         = new TH1F("muonPhi_num"         ,"muonPhi_num"          ,   40, -3.2,   3.2);

  TH1F* nvtx_event          = new TH1F("nvtx_event"          ,"nvtx_event"           ,   60,  0,  60);
  TH1F* nvtx_num            = new TH1F("nvtx_num"            ,"nvtx_num"             ,   60,  0,  60);
  TH1F* nvtx_den            = new TH1F("nvtx_den"            ,"nvtx_den"             ,   60,  0,  60);
  TH1F* nvtx_barrel_num     = new TH1F("nvtx_barrel_num"     ,"nvtx_barrel_num"      ,   60,  0,  60);
  TH1F* nvtx_barrel_den     = new TH1F("nvtx_barrel_den"     ,"nvtx_barrel_den"      ,   60,  0,  60);
  TH1F* nvtx_endcap_num     = new TH1F("nvtx_endcap_num"     ,"nvtx_endcap_num"      ,   60,  0,  60);
  TH1F* nvtx_endcap_den     = new TH1F("nvtx_endcap_den"     ,"nvtx_endcap_den"      ,   60,  0,  60);
  TH1F* nvtx_eta0_num       = new TH1F("nvtx_eta0_num"       ,"nvtx_eta0_num"        ,   60,  0,  60);
  TH1F* nvtx_eta0_den       = new TH1F("nvtx_eta0_den"       ,"nvtx_eta0_den"        ,   60,  0,  60);
  TH1F* nvtx_eta1_num       = new TH1F("nvtx_eta1_num"       ,"nvtx_eta1_num"        ,   60,  0,  60);
  TH1F* nvtx_eta1_den       = new TH1F("nvtx_eta1_den"       ,"nvtx_eta1_den"        ,   60,  0,  60);
  TH1F* nvtx_eta2_num       = new TH1F("nvtx_eta2_num"       ,"nvtx_eta2_num"        ,   60,  0,  60);
  TH1F* nvtx_eta2_den       = new TH1F("nvtx_eta2_den"       ,"nvtx_eta2_den"        ,   60,  0,  60);

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
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    
    unsigned int nmuons = ev->muons.size();
    if (nmuons < 2) continue;
    
    unsigned int nhltmuons = ev->hltmuons.size();
//     if (nhltmuons > 0) std::cout << "Number of hlt muons = " << nhltmuons << std::endl;
    
    if (!ev-> hlt.find(hltname)) continue;
    nvtx_event      -> Fill( ev -> nVtx   );
    

    for (int imu = 0; imu < nmuons; imu++){
      
      // select the tag muon        
      if (! selectTagMuon(ev -> muons.at(imu), tagiso))                   continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hlt.objects, isofilter)) continue;
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
        
        if (!(offlineiso04 < offlineIsoCut)) continue;
//         if (!(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, L3filter))) continue;

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

        if (fabs(ev -> muons.at(jmu).eta) < 0.9){
          muonPt_eta0_den  -> Fill( ev -> muons.at(jmu).pt );
          nvtx_eta0_den    -> Fill( ev -> nVtx );
        }
        else if (fabs(ev -> muons.at(jmu).eta) > 0.9 && fabs(ev -> muons.at(jmu).eta) < 1.2){
          muonPt_eta1_den  -> Fill( ev -> muons.at(jmu).pt );
          nvtx_eta1_den    -> Fill( ev -> nVtx );
        }
        else{
          muonPt_eta2_den  -> Fill( ev -> muons.at(jmu).pt );
          nvtx_eta2_den    -> Fill( ev -> nVtx );
        }


        if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, isofilter)){
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

		  if (fabs(ev -> muons.at(jmu).eta) < 0.9){
			muonPt_eta0_num  -> Fill( ev -> muons.at(jmu).pt );
			nvtx_eta0_num    -> Fill( ev -> nVtx );
		  }
		  else if (fabs(ev -> muons.at(jmu).eta) > 0.9 && fabs(ev -> muons.at(jmu).eta) < 1.2){
			muonPt_eta1_num  -> Fill( ev -> muons.at(jmu).pt );
			nvtx_eta1_num    -> Fill( ev -> nVtx );
		  }
		  else{
			muonPt_eta2_num  -> Fill( ev -> muons.at(jmu).pt );
			nvtx_eta2_num    -> Fill( ev -> nVtx );
		  }


        }
      }
      
    }

  }
  
  outfile            -> cd();
  muonPt_den         -> Write();
  muonPt_num         -> Write();
  muonPt_barrel_den  -> Write();
  muonPt_barrel_num  -> Write();
  muonPt_endcap_den  -> Write();
  muonPt_endcap_num  -> Write();
  muonPt_eta0_den    -> Write();
  muonPt_eta0_num    -> Write();
  muonPt_eta1_den    -> Write();
  muonPt_eta1_num    -> Write();
  muonPt_eta2_den    -> Write();
  muonPt_eta2_num    -> Write();
  tagMuonPt          -> Write();
  
  muonEta_den        -> Write();
  muonEta_num        -> Write();
  muonPhi_den        -> Write();
  muonPhi_num        -> Write();
  
  nvtx_event         -> Write();
  nvtx_num           -> Write();
  nvtx_den           -> Write();
  nvtx_barrel_num    -> Write();
  nvtx_barrel_den    -> Write();
  nvtx_endcap_num    -> Write();
  nvtx_endcap_den    -> Write();
  nvtx_eta0_num      -> Write();
  nvtx_eta0_den      -> Write();
  nvtx_eta1_num      -> Write();
  nvtx_eta1_den      -> Write();
  nvtx_eta2_num      -> Write();
  nvtx_eta2_den      -> Write();
  
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
  
  if (!( mu.pt         > 18  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0.,
                       mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  tagiso -> Fill(offlineiso04);
  if (!(offlineiso04   < offlineIsoCut)) return false;
  
  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt         > 18  )) return false; 
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



