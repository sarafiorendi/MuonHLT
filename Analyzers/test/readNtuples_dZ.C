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
HLTMuonCand  matchL3(MuonCand, std::vector<HLTMuonCand>);

std::string L3filter_17   = "hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17::HLT";
std::string L3filter_8    = "hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8::HLT";
std::string isofilter     = "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4::HLT";

double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double eta_bins[14] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4};

double offlineIsoCut = 1000;

void readNtuples_dZ(){

  std::string hltname = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4";

  TFile* inputfile = TFile::Open("root://eoscms//eos/cms/store/group/phys_muon/fiorendi/13TeV/2016D/SingleMuon/crab_promptReco_v2_goldenJson_271036-277148_asMenu/160803_163959/0000/muonNtuple_983.root","READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open(Form("testDimuon.root", hltname.c_str()),"RECREATE");
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

  TH1F* dZ                      = new TH1F("dZ"                     ,"dZ"               , 2000,-10,   10 );
  TH1F* z_mu17                  = new TH1F("z_mu17"                 ,"z_mu17"           , 1000,-10,   10 );
  TH1F* z_mu8                   = new TH1F("z_mu8"                  ,"z_mu8"            , 1000,-10,   10 );
 
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
      
      MuonCand Mu17Cand = ev -> muons.at(imu);
      // select the tag muon        
      if (! selectTagMuon(Mu17Cand, tagiso))                      continue;
//       if (! matchMuon(Mu17Cand, ev -> hlt.objects, L3filter_17))  continue;
//       if (! matchMuon(Mu17Cand, ev -> hlt.objects, isofilter))    continue;
      tagMuonPt -> Fill(Mu17Cand.pt);
      
      for (int jmu = 0; jmu < nmuons; jmu++){

        bool pass   = false;

        MuonCand Mu8Cand = ev -> muons.at(jmu);
        // select the probe muon
        if (! selectProbeMuon(Mu8Cand, Mu17Cand, dimuon_mass)) continue;

        offlineiso04 = Mu8Cand.chargedDep_dR04 + std::max(0.,
                       Mu8Cand.photonDep_dR04 + Mu8Cand.neutralDep_dR04 - 0.5*Mu8Cand.puPt_dR04);
        offlineiso04 = offlineiso04 / Mu8Cand.pt;
        
//         if (!(offlineiso04 < offlineIsoCut)) continue;
//         if (!(matchMuon(Mu8Cand, ev -> hlt.objects, L3filter_8))) continue;
//         if (!(matchMuon(Mu8Cand, ev -> hlt.objects, isofilter ))) continue;
        
        HLTMuonCand theL3_17 = matchL3(Mu17Cand, ev -> hltmuons);
        HLTMuonCand theL3_8  = matchL3(Mu8Cand , ev -> hltmuons);
//         if (matchMuon(Mu8Cand, ev -> hlt.objects, "hltL3fL1sMu16L1f0L2f10QL3Filtered20Q::REHLT") 
//             && theL3.pt != -1000
//             && theL3.pt > cut_pt
//            ) 
//           foundL3 = true;

        z_mu17 -> Fill( theL3_17.vz);
        z_mu8  -> Fill( theL3_8 .vz);

        // match probe muon to the interesting filter and fill numerator histograms
//         if (matchMuon(Mu8Cand, ev -> hlt.objects, isofilter))       pass = true; 

//         muonPt       -> Fill( pass, Mu8Cand.pt );
//         muonEta      -> Fill( pass, Mu8Cand.eta);
//         muonPhi      -> Fill( pass, Mu8Cand.phi);
//         muonIso      -> Fill( pass, offlineiso04           );
//         nvtx         -> Fill( pass, ev -> nVtx             );
//         
//         if (fabs(Mu8Cand.eta) < 1.479){
//           muonPt_barrel     -> Fill( pass, Mu8Cand.pt );
//           nvtx_barrel       -> Fill( pass, ev -> nVtx );
//           muonIso_barrel    -> Fill( pass, offlineiso04 );
//         }
//         else{
//           muonPt_endcap     -> Fill( pass, Mu8Cand.pt );
//           nvtx_endcap       -> Fill( pass, ev -> nVtx             );
//           muonIso_endcap    -> Fill( pass, offlineiso04           );
//         }
// 
      }
    }
  }
  
  outfile           -> cd();
  
  z_mu17            -> Write();        
  z_mu8             -> Write();        

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
  
//   muonIso           -> Write();
//   muonIso_barrel    -> Write();
//   muonIso_endcap    -> Write();

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
//     std::cout << it->filterTag << std::endl;
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
  if (!(offlineiso04   < offlineIsoCut)) return false;
  
  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt         >  8  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
//   if (mu.charge * tagMu.charge > 0) return false;
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();
  dimuon_mass -> Fill(mumumass);
//   if (! (mumumass > 86. && mumumass < 96. )) return false;
  
  return true;
}


HLTMuonCand matchL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

  bool match = false;
  int nL3 = L3cands.size();

  float minDR = 0.1;
  float theDR = 100;
  HLTMuonCand theL3;
  theL3.pt        = -1000;
  theL3.eta       = -1000;
  theL3.phi       = -1000;
  theL3.trkpt     = -1000;
  theL3.ecalDep   = -1000;
  theL3.hcalDep   = -1000;
  theL3.trkDep    = -1000;
  theL3.vz        = -1000;
  theL3.n_pix_hit = -1000;
  theL3.n_trk_hit = -1000;
  theL3.n_pix_lay = -1000;
  theL3.n_trk_lay = -1000;
  
  for ( std::vector<HLTMuonCand>::const_iterator it = L3cands.begin(); it != L3cands.end(); ++it ) {
      
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
      theL3 = *it;
    }
  }
  
  return theL3;
}

