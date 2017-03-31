#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TEfficiency.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLT/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

double    muonmass  = 0.10565837;
std::string wpstring;  

// define pt threshold
float cut_pt          = 24  ; 
float cut_pt_offline  = 27  ; 
// define effective areas

float cut_trk         = 0;

std::string L3filter    = "hltL3fL1sMu22L1f0L2f10QL3Filtered24Q::REHLT";
std::string isofilter   = "hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09::HLT";

bool         selectTagMuon  (MuonCand );
bool         selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool         matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
HLTMuonCand  matchL3        (MuonCand, std::vector<HLTMuonCand>);
bool         matchGENMuon   (MuonCand, std::vector<GenParticleCand>);

double pt_bins[17]  = { 0, 15, 18, 20, 22, 24, 26, 28, 30, 35, 40, 50, 60, 80, 120, 300};// 1000} ;
double eta_bins[14] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4}    ;


void readNtuplesMC(){

  cut_trk         = 0.09    ;

//   TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/CMSSW_9_0_0_pre6/src/HLTrigger/Configuration/test/muonNtuple_samePtsas2016.root","READ");
//   std::cout << "input file: " << inputfile -> GetName() << std::endl;
  TChain* tree = new TChain("muonNtuples/muonTree");
  for (int i = 0; i < 50; i++){
    tree -> Add(Form("root://eoscms//eos/cms/store/group/phys_muon/fiorendi/13TeV/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TkIsolation_2017_81X_split/170327_192245/0000/muonNtuple_%d.root/muonNtuples/muonTree", i));
  }

  TFile* outfile = TFile::Open("efficiency_onDYToLLTSG.root","RECREATE");
//   TFile* outfile = TFile::Open("efficiency_onDYToLLRelVal_samePtAs2016.root","RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
//   TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree");
//   if (!tree) {
//     std::cout << " *** tree not found *** " << std::endl;
//     return;
//   }
    
  TH1F* dimuon_mass             = new TH1F("dimuon_mass"            ,"dimuon_mass"          , 1500,  0, 150);
  TH1F* tagMuonPt               = new TH1F("tagMuonPt"              ,"tagMuonPt"            ,  150,  0, 150);
  TH1F* nvtx_event              = new TH1F("nvtx_event"             ,"nvtx_event"           ,   60,  0,  60);
  TH1F* den_pt                  = new TH1F("den_pt"                 ,"den_pt"               ,  150,  0, 150);
  TH1F* num_pt                  = new TH1F("num_pt"                 ,"num_pt"               ,  150,  0, 150);


  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"               ,  15, pt_bins);//18
  TEfficiency* muonPt_barrel    = new TEfficiency("muonPt_barrel"   ,"muonPt_barrel"        ,  15, pt_bins);
  TEfficiency* muonPt_endcap    = new TEfficiency("muonPt_endcap"   ,"muonPt_endcap"        ,  15, pt_bins);
     
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"              ,  13, eta_bins);
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"              ,  20, -3.2,3.2);
    
  TEfficiency* nvtx             = new TEfficiency("nvtx"            , "nvtx"                ,   30,  0,  60);
  TEfficiency* nvtx_barrel      = new TEfficiency("nvtx_barrel"     , "nvtx_barrel"         ,   30,  0,  60);
  TEfficiency* nvtx_endcap      = new TEfficiency("nvtx_endcap"     , "nvtx_endcap"         ,   30,  0,  60);

  TEfficiency* muonEta2016      = new TEfficiency("muonEta2016"     ,"muonEta2016"          ,  13, eta_bins);
  TEfficiency* muonPhi2016      = new TEfficiency("muonPhi2016"     ,"muonPhi2016"          ,  20, -3.2,3.2);
  TEfficiency* nvtx2016         = new TEfficiency("nvtx2016"        , "nvtx2016"            ,  30,  0,  60 );
  TEfficiency* muonPt2016       = new TEfficiency("muonPt2016"      ,"muonP2016t"           ,  15,  pt_bins);//18
  
  TH1F* trkIso_barrel           = new TH1F("trkIso_barrel"          ,"trkIso_barrel"        ,   200,  0,  0.2);
  TH1F* trkIso_endcap           = new TH1F("trkIso_endcap"          ,"trkIso_endcap"        ,   200,  0,  0.2);
  TH1F* trkIso2016_barrel       = new TH1F("trkIso2016_barrel"      ,"trkIso2016_barrel"    ,   200,  0,  0.2);
  TH1F* trkIso2016_endcap       = new TH1F("trkIso2016_endcap"      ,"trkIso2016_endcap"    ,   200,  0,  0.2);

  TH2F* trkDeps           = new TH2F("trkDeps"          ,"trkDeps"        ,   200,  0,  20, 200, 0,20);

  double offlineiso04 = 100;

  MuonEvent* ev = new MuonEvent();
  tree -> SetBranchStatus("*",1);
  tree -> SetBranchAddress( "event", &ev);
//   TBranch* evBranch = tree->GetBranch("event");
//   evBranch -> SetAddress(&ev);

  int nentries = tree->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  float w = 1;
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    
    unsigned int nmuons    = ev->muons.size();
    unsigned int nhltmuons = ev->hltmuons.size();
    if (nmuons < 2) continue;
    
    if (!ev-> hltTag.find("HLT_IsoMu24_v4")) continue;

    nvtx_event -> Fill( ev -> nVtx );
//     if ( ev -> trueNI < 28 || ev -> trueNI > 32 ) continue;

    for (int imu = 0; imu < nmuons; imu++){
      
      // select the tag muon        
      if (! selectTagMuon(ev -> muons.at(imu)))                          continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hltTag.objects, isofilter)) continue;
      tagMuonPt -> Fill(ev -> muons.at(imu).pt, w);
      
      for (int jmu = 0; jmu < nmuons; jmu++){
        MuonCand probeMu = ev -> muons.at(jmu);
        // select the probe muon
        if (! selectProbeMuon(probeMu, ev -> muons.at(imu), dimuon_mass)) continue;
//         if (! matchGENMuon   (probeMu, ev -> genParticles))               continue;
        bool foundL3         = false;
        bool pass            = false;
        bool passPt          = false;
        bool pass2016        = false;
        bool pass2016Pt      = false;
        bool isBarrel        = false;

        offlineiso04 = probeMu.chargedDep_dR04 + std::max(0.,
                       probeMu.photonDep_dR04 + probeMu.neutralDep_dR04 - 0.5*probeMu.puPt_dR04);
        offlineiso04 = offlineiso04 / probeMu.pt;

        if (!(offlineiso04 < 0.15)) continue;
        if (! matchMuon(probeMu, ev -> hlt.objects, "hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22::REHLT")) continue;

        // to eval eficiency of the path in the menu 
//         if (matchMuon(probeMu, ev -> hlt.objects, "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09::REHLT")) pass = true;

		if (fabs(probeMu.eta) < 1.479)  isBarrel = true;

        HLTMuonCand theL3 = matchL3(probeMu, ev -> hltmuons);
//         if (matchMuon(probeMu, ev -> hlt.objects, "hltL3fL1sMu22L1f0L2f10QL3Filtered24QL3pfhcalIsoRhoFilteredHB0p21HE0p22::REHLT") 
//      
        if ( theL3.pt != -1000
              && theL3.pt > cut_pt
           ) 
        foundL3 = true;
        
        if (!foundL3) continue;  
        
        if (probeMu.pt < 0 || probeMu.pt > 1000) continue;
        den_pt -> Fill( probeMu.pt);  

        double hltTrkIso      = (theL3.trkDep) / theL3.pt;
        double hltTrkIso2016  = (theL3.trkDep2016) / theL3.pt;
        
        if ( hltTrkIso < cut_trk && probeMu.pt > cut_pt_offline  ) pass = true;
        if ( hltTrkIso < cut_trk )                                 passPt = true;

        if ( hltTrkIso2016 < cut_trk && probeMu.pt > cut_pt_offline  ) pass2016 = true;
        if ( hltTrkIso2016 < cut_trk )                                 pass2016Pt = true;

        muonPt      -> Fill( passPt    , probeMu.pt  );
        muonPt2016  -> Fill( pass2016Pt, probeMu.pt  );
        if (probeMu.pt > cut_pt_offline ){
          muonEta       -> Fill( pass,      probeMu.eta );
          muonEta2016   -> Fill( pass2016,  probeMu.eta );
          muonPhi       -> Fill( pass,      probeMu.phi );
          muonPhi2016   -> Fill( pass2016,  probeMu.phi );
          nvtx          -> Fill( pass,      ev -> nVtx  );
          nvtx2016      -> Fill( pass2016,  ev -> nVtx  );
          
          trkDeps  -> Fill( theL3.trkDep, theL3.trkDep2016);
        }
        if (isBarrel){
          muonPt_barrel    -> Fill( passPt, probeMu.pt  );
          if (probeMu.pt > cut_pt_offline ) {
            nvtx_barrel        -> Fill( pass  , ev -> nVtx   );
            trkIso_barrel      -> Fill        ( hltTrkIso    );
            trkIso2016_barrel  -> Fill        ( hltTrkIso2016);
          }
        }
        
        else{
          muonPt_endcap    -> Fill( passPt, probeMu.pt  );
          if (probeMu.pt > cut_pt_offline ){
            nvtx_endcap       -> Fill( pass, ev -> nVtx);
            trkIso_endcap     -> Fill( hltTrkIso       );
            trkIso2016_endcap -> Fill( hltTrkIso2016   );
          }
        }
        
        
      }
      
    }



  }
  



  outfile         -> cd();
  muonPt          -> Write();
  muonPt_barrel   -> Write();
  muonPt_endcap   -> Write();
  tagMuonPt       -> Write();
  trkDeps         -> Write();
  
  den_pt          -> Write();
  num_pt          -> Write();
  
  muonEta         -> Write();
  muonPhi         -> Write();

  muonPt2016      -> Write();
  muonEta2016     -> Write();
  muonPhi2016     -> Write();
  nvtx2016        -> Write();
  
  nvtx_event      -> Write();
  nvtx            -> Write();
  nvtx_barrel     -> Write();
  nvtx_endcap     -> Write();
  
  trkIso_barrel     -> Write();
  trkIso_endcap     -> Write();
  trkIso2016_barrel -> Write();
  trkIso2016_endcap -> Write();
  

  dimuon_mass     -> Write();

  outfile         -> Close();  
  
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
  
  if (!( mu.pt         > 26  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0.,
                       mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  if (!(offlineiso04   < 0.15 )) return false;

  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
//   if (!( mu.pt         > cut_pt_offline  )) return false; 
  if (!( fabs(mu.eta)  < 2.4             )) return false; 
  if (!( mu.isTight    == 1              )) return false; 
  
  if (mu.charge * tagMu.charge > 0) return false;
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();
  dimuon_mass -> Fill(mumumass);
  if (! (mumumass > 86. && mumumass < 96. )) return false;
  
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
  theL3.ecalDep05 = -1000;
  theL3.hcalDep05 = -1000;
  theL3.ecalDep1  = -1000;
  theL3.hcalDep1  = -1000;
  
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



bool matchGENMuon(MuonCand mu, std::vector<GenParticleCand> gen){

  bool match = false;

  float minDR = 0.5;
  float theDR = 100;
  for ( std::vector<GenParticleCand>::const_iterator it = gen.begin(); it != gen.end(); ++it ) {
// 	if (it -> status != 1 ) continue;
	theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
	if (theDR < minDR){
	  minDR = theDR;
	  match = true;
    }
  }
  return match;
}


