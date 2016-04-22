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

double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double eta_bins[14] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};
double iso_bins[12] = { 0  , 0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1   };

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


  TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/CMSSW_8_0_0/src/MuonHLT/Analyzers/test/muonNtuple_run258158.root","READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open(Form("efficiency_fullpath_onZMuMuSkim_run258158_%s_offlineIso0p15.root", hltname.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree");
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }
    
  TH1F* dimuon_mass             = new TH1F("dimuon_mass"            ,"dimuon_mass"      , 1500,  0,  150 );
  TH1F* tagiso                  = new TH1F("tagiso"                 ,"tagiso"           ,  100,  0,  1   );
  TH1F* tagMuonPt               = new TH1F("tagMuonPt"              ,"tagMuonPt"        ,  150,  0,  150 );
  TH1F* nvtx_event              = new TH1F("nvtx_event"             ,"nvtx_event"       ,   60,  0,   60 );
 
  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"           ,   11,  pt_bins );
  TEfficiency* muonPt_barrel    = new TEfficiency("muonPt_barrel"   ,"muonPt_barrel"    ,   11,  pt_bins );
  TEfficiency* muonPt_endcap    = new TEfficiency("muonPt_endcap"   ,"muonPt_endcap"    ,   11,  pt_bins );
  TEfficiency* muonPt_eta0      = new TEfficiency("muonPt_eta0"     ,"muonPt_eta0"      ,   11,  pt_bins );
  TEfficiency* muonPt_eta1      = new TEfficiency("muonPt_eta1"     ,"muonPt_eta1"      ,   11,  pt_bins );
  TEfficiency* muonPt_eta2      = new TEfficiency("muonPt_eta2"     ,"muonPt_eta2"      ,   11,  pt_bins );
  
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"          ,   13, eta_bins );
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"          ,   40, -3.2, 3.2);
 
  TEfficiency* nvtx             = new TEfficiency("nvtx"            ,"nvtx"             ,   60,    0,  60);
  TEfficiency* nvtx_barrel      = new TEfficiency("nvtx_barrel"     ,"nvtx_barrel"      ,   60,    0,  60);
  TEfficiency* nvtx_endcap      = new TEfficiency("nvtx_endcap"     ,"nvtx_endcap"      ,   60,    0,  60);
  TEfficiency* nvtx_eta0        = new TEfficiency("nvtx_eta0"       ,"nvtx_eta0"        ,   60,    0,  60);
  TEfficiency* nvtx_eta1        = new TEfficiency("nvtx_eta1"       ,"nvtx_eta1"        ,   60,    0,  60);
  TEfficiency* nvtx_eta2        = new TEfficiency("nvtx_eta2"       ,"nvtx_eta2"        ,   60,    0,  60);
   
  TEfficiency* muonIso          = new TEfficiency("muonIso"         ,"muonIso"          ,  100,    0,   1);
  TEfficiency* muonIso_barrel   = new TEfficiency("muonIso_barrel"  ,"muonIso_barrel"   ,  100,    0,   1);
  TEfficiency* muonIso_endcap   = new TEfficiency("muonIso_endcap"  ,"muonIso_endcap"   ,  100,    0,   1);
  
  double offlineiso04 = 100;

  MuonEvent* ev      = new MuonEvent();
  TBranch*  evBranch = tree->GetBranch("event");
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
      if (! selectTagMuon(ev -> muons.at(imu), tagiso))                      continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hlt.objects, isofilter))    continue;
      tagMuonPt -> Fill(ev -> muons.at(imu).pt);
      
      for (int jmu = 0; jmu < nmuons; jmu++){

        bool pass   = false;

        // select the probe muon
        if (! selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;

        offlineiso04 = ev -> muons.at(jmu).chargedDep_dR04 + std::max(0.,
                       ev -> muons.at(jmu).photonDep_dR04 + ev -> muons.at(jmu).neutralDep_dR04 - 0.5*ev -> muons.at(jmu).puPt_dR04);
        offlineiso04 = offlineiso04 / ev -> muons.at(jmu).pt;
        
        if (!(offlineiso04 < offlineIsoCut)) continue;
        if (!(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, L3filter))) continue;

        // match probe muon to the interesting filter and fill numerator histograms
        if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, isofilter))       pass = true; 

        muonPt       -> Fill( pass, ev -> muons.at(jmu).pt );
        muonEta      -> Fill( pass, ev -> muons.at(jmu).eta);
        muonPhi      -> Fill( pass, ev -> muons.at(jmu).phi);
        muonIso      -> Fill( pass, offlineiso04           );
        nvtx         -> Fill( pass, ev -> nVtx             );
        
        if (fabs(ev -> muons.at(jmu).eta) < 1.479){
          muonPt_barrel     -> Fill( pass, ev -> muons.at(jmu).pt );
          nvtx_barrel       -> Fill( pass, ev -> nVtx );
          muonIso_barrel    -> Fill( pass, offlineiso04 );
        }
        else{
          muonPt_endcap     -> Fill( pass, ev -> muons.at(jmu).pt );
          nvtx_endcap       -> Fill( pass, ev -> nVtx             );
          muonIso_endcap    -> Fill( pass, offlineiso04           );
        }

        if (fabs(ev -> muons.at(jmu).eta) < 0.9){
          muonPt_eta0     -> Fill( pass, ev -> muons.at(jmu).pt );
          nvtx_eta0       -> Fill( pass, ev -> nVtx );
        }
        else if (fabs(ev -> muons.at(jmu).eta) > 0.9 && fabs(ev -> muons.at(jmu).eta) < 1.2){
          muonPt_eta1     -> Fill( pass, ev -> muons.at(jmu).pt );
          nvtx_eta1       -> Fill( pass, ev -> nVtx );
        }
        else{
          muonPt_eta2     -> Fill( pass, ev -> muons.at(jmu).pt );
          nvtx_eta2       -> Fill( pass, ev -> nVtx );
        }
      }
    }
  }
  
  outfile           -> cd();
  muonPt            -> Write();
  muonPt_barrel     -> Write();
  muonPt_endcap     -> Write();
  muonPt_eta0       -> Write();
  muonPt_eta1       -> Write();
  muonPt_eta2       -> Write();
  tagMuonPt         -> Write();
  
  muonEta           -> Write();
  muonPhi           -> Write();
  
  nvtx_event        -> Write();
  nvtx              -> Write();
  nvtx_barrel       -> Write();
  nvtx_endcap       -> Write();
  nvtx_eta0         -> Write();
  nvtx_eta1         -> Write();
  nvtx_eta2         -> Write();
  
  muonIso           -> Write();
  muonIso_barrel    -> Write();
  muonIso_endcap    -> Write();

  dimuon_mass       -> Write();
  tagiso            -> Write();
  
  outfile           -> Close();  
  
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



