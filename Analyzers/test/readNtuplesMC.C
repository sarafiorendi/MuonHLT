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
const int nhist_iso = 1;

// old wp
// EB 0.11 EE 0.08 HB 0.21 HE 0.22 TRK 0.09
// new loose wp
// EB 0.08 EE 0.06 HB 0.13 HE 0.13 TRK 0.08

// define pt threshold
float cut_pt          = 24; 
float cut_pt_offline  = 25; 
// define effective areas
// 2015 menu
// float a_ecal_barrel   = 0.153;
// float a_hcal_barrel   = 0.060;
// float a_ecal_endcap   = 0.072;
// float a_hcal_endcap   = 0.107;
// 2016 menu
float a_ecal_barrel   = 0.153;
float a_hcal_barrel   = 0.074;
float a_ecal_endcap   = 0.071;
float a_hcal_endcap   = 0.100;

// define isolation thresholds
float cut_ecal_barrel = 0.11 ;
float cut_hcal_barrel = 0.13 ;
float cut_ecal_endcap = 0.08 ; 
float cut_hcal_endcap = 0.13 ; 
float cut_trk         = 0.08 ;

bool         selectTagMuon  (MuonCand );
bool         selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool         matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
HLTMuonCand  matchL3        (MuonCand, std::vector<HLTMuonCand>);
bool         matchGENMuon   (MuonCand, std::vector<GenParticleCand>);

double pt_bins[12]  = { 20, 24, 27, 30, 35, 40, 45, 50, 60, 70 , 90, 150} ;
double eta_bins[14] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4}    ;
double iso_bins[12] = {0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1}                   ;


void readNtuplesMC(){

  TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_DYToLL.root","READ");
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open("efficiency_onDYToLL_MC_allnPU_pt24_oldEcal_newHcal.root","RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree");
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }

//   TFile* wfile = TFile::Open("weightsPU_25ns.root","READ");
//   if (!wfile) std::cout << "file not found!" << std::endl;
//   TH1F*  whist = (TH1F*) wfile -> Get("data");
//   if (!wfile) std::cout << "hist not found!" << std::endl;

  TH1F* dimuon_mass         = new TH1F("dimuon_mass"         ,"dimuon_mass"          , 1500,  0, 150);
  TH1F* tagMuonPt           = new TH1F("tagMuonPt"           ,"tagMuonPt"            ,  150,  0, 150);
  TH1F* nvtx_event          = new TH1F("nvtx_event"          ,"nvtx_event"           ,   60,  0,  60);
  TH1F* den_pt              = new TH1F("den_pt"              ,"den_pt"               ,  150,  0, 150);
  TH1F* num_pt              = new TH1F("num_pt"              ,"num_pt"               ,  150,  0, 150);


  TEfficiency* muonPt          = new TEfficiency("muonPt"          ,"muonPt"           ,  11, pt_bins);//18
  TEfficiency* muonPt_barrel   = new TEfficiency("muonPt_barrel"   ,"muonPt_barrel"    ,  11, pt_bins);
  TEfficiency* muonPt_endcap   = new TEfficiency("muonPt_endcap"   ,"muonPt_endcap"    ,  11, pt_bins);

  TEfficiency* muonEta         = new TEfficiency("muonEta"         ,"muonEta"          ,  13, eta_bins);
  TEfficiency* muonPhi         = new TEfficiency("muonPhi"         ,"muonPhi"          ,  20, -3.2,3.2);

  TEfficiency* nvtx            = new TEfficiency("nvtx"           ,"nvtx"              ,   30,  0,  60);
  TEfficiency* nvtx_barrel     = new TEfficiency("nvtx_barrel"    ,"nvtx_barrel"       ,   30,  0,  60);
  TEfficiency* nvtx_endcap     = new TEfficiency("nvtx_endcap"    ,"nvtx_endcap"       ,   30,  0,  60);

  TEfficiency* muonIso         = new TEfficiency("muonIso"        ,"muonIso"           ,   11, iso_bins);
  TEfficiency* muonIso_barrel  = new TEfficiency("muonIso_barrel" ,"muonIso_barrel"    ,   11, iso_bins);
  TEfficiency* muonIso_endcap  = new TEfficiency("muonIso_endcap" ,"muonIso_endcap"    ,   11, iso_bins);
  TEfficiency* muonIso04       = new TEfficiency("muonIso04"      ,"muonIso04"         ,   11, iso_bins);

  TH1F* ecalIso_barrel      = new TH1F("ecalIso_barrel"      ,"ecalIso_barrel"       ,   50,  0,  0.5);
  TH1F* hcalIso_barrel      = new TH1F("hcalIso_barrel"      ,"hcalIso_barrel"       ,   50,  0,  0.5);
  TH1F* ecalIso_endcap      = new TH1F("ecalIso_endcap"      ,"ecalIso_endcap"       ,   50,  0,  0.5);
  TH1F* hcalIso_endcap      = new TH1F("hcalIso_endcap"      ,"hcalIso_endcap"       ,   50,  0,  0.5);

  TEfficiency *h_ecal_iso_barrel    [nhist_iso];
  TEfficiency *h_hcal_iso_barrel    [nhist_iso];
  TEfficiency *h_trk_iso_barrel     [nhist_iso];

  TEfficiency *h_ecal_iso_endcap    [nhist_iso];
  TEfficiency *h_hcal_iso_endcap    [nhist_iso];
  TEfficiency *h_trk_iso_endcap     [nhist_iso];

  for (int i=0; i< nhist_iso; i++) {
    h_ecal_iso_barrel[i] = new TEfficiency(Form("HMuonPt_barrel_ecal_iso%d",i), "", 1000, 0, 1000);
    h_hcal_iso_barrel[i] = new TEfficiency(Form("HMuonPt_barrel_hcal_iso%d",i), "", 1000, 0, 1000);
    h_trk_iso_barrel[i]  = new TEfficiency(Form("HMuonPt_barrel_trk_iso%d" ,i), "", 1000, 0, 1000);
    h_ecal_iso_endcap[i] = new TEfficiency(Form("HMuonPt_endcap_ecal_iso%d",i), "", 1000, 0, 1000);
    h_hcal_iso_endcap[i] = new TEfficiency(Form("HMuonPt_endcap_hcal_iso%d",i), "", 1000, 0, 1000);
    h_trk_iso_endcap[i]  = new TEfficiency(Form("HMuonPt_endcap_trk_iso%d" ,i), "", 1000, 0, 1000);
  }

    
  double offlineiso   = 100;
  double offlineiso04 = 100;
  float a_ecal, a_hcal;
  float cut_ecal, cut_hcal;

  MuonEvent* ev = new MuonEvent();
  TBranch* evBranch = tree->GetBranch("event");
  evBranch -> SetAddress(&ev);

  int nentries = tree->GetEntriesFast();
  std::cout << "Number of entries = " << nentries << std::endl;

  float w = 1;
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
//     std::cout << " ---- analyzing event: " << ev -> eventNumber << std::endl;
    
    unsigned int nmuons    = ev->muons.size();
    unsigned int nhltmuons = ev->hltmuons.size();
    if (nmuons < 2) continue;
    
    if (!ev-> hlt.find("HLT_Mu20_v2")) continue;
// 	int wbin = whist->FindBin(ev -> nVtx); 
// 	w        = whist->GetBinContent(wbin); 
    nvtx_event -> Fill( ev -> nVtx );
//     if ( ev -> trueNI < 28 || ev -> trueNI > 32 ) continue;

    for (int imu = 0; imu < nmuons; imu++){
      
      // select the tag muon        
      if (! selectTagMuon(ev -> muons.at(imu)))                               continue;
      if (! matchGENMuon (ev -> muons.at(imu), ev -> genParticles))           continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hlt.objects, "hltL3fL1sMu16L1f0L2f10QL3Filtered20Q::REHLT")) continue;
      tagMuonPt -> Fill(ev -> muons.at(imu).pt, w);
      
      for (int jmu = 0; jmu < nmuons; jmu++){
        // select the probe muon
        if (! selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
        if (! matchGENMuon   (ev -> muons.at(jmu), ev -> genParticles))               continue;
        bool foundL3  = false;
        bool pass     = false;
        bool passECal = false;
        bool passHCal = false;
        bool passTrk  = false;
        bool isBarrel = false;
        bool passECalNominal   = false;
        bool passHCalNominal   = false;

        offlineiso = ev -> muons.at(jmu).chargedDep_dR03 + std::max(0.,
                     ev -> muons.at(jmu).photonDep_dR03 + ev -> muons.at(jmu).neutralDep_dR03 - 0.5*ev -> muons.at(jmu).puPt_dR03);
        offlineiso = offlineiso / ev -> muons.at(jmu).pt;

        offlineiso04 = ev -> muons.at(jmu).chargedDep_dR04 + std::max(0.,
                       ev -> muons.at(jmu).photonDep_dR04 + ev -> muons.at(jmu).neutralDep_dR04 - 0.5*ev -> muons.at(jmu).puPt_dR04);
        offlineiso04 = offlineiso04 / ev -> muons.at(jmu).pt;

        if (!(offlineiso04 < 0.15)) continue;

        // to eval eficiency of the path in the menu 
//         if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09::REHLT")) pass = true;

		if (fabs(ev -> muons.at(jmu).eta) < 1.479){
	      isBarrel = true;
		  a_ecal   = a_ecal_barrel;
		  a_hcal   = a_hcal_barrel;
		  cut_hcal = cut_ecal_barrel;
		  cut_ecal = cut_hcal_barrel;
		}
		else{
		  a_ecal   = a_ecal_endcap;
		  a_hcal   = a_hcal_endcap;
		  cut_ecal = cut_ecal_endcap;
		  cut_hcal = cut_hcal_endcap;
		}


        HLTMuonCand theL3 = matchL3(ev -> muons.at(jmu), ev -> hltmuons);
        if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, "hltL3fL1sMu16L1f0L2f10QL3Filtered20Q::REHLT") 
            && theL3.pt != -1000
            && theL3.pt > cut_pt
           ) 
          foundL3 = true;
        
        if (!foundL3) continue;  
        den_pt -> Fill( ev -> muons.at(jmu).pt);  

        double hltEcalIso = max(0., theL3.ecalDep05   - a_ecal * ev -> hlt.rho ) / theL3.pt;
        double hltHcalIso = max(0., theL3.hcalDep1    - a_hcal * ev -> hlt.rho ) / theL3.pt;
        double hltTrkIso  = (theL3.trkDep) / theL3.pt;
        
        if (hltEcalIso < cut_ecal && hltHcalIso < cut_hcal && hltTrkIso < cut_trk && foundL3) pass = true;
        if (hltEcalIso < cut_ecal && foundL3) passECalNominal = true;
        if (hltHcalIso < cut_hcal && foundL3) passHCalNominal = true;
 
        muonPt    -> FillWeighted( pass, w, ev -> muons.at(jmu).pt  );
        muonEta   -> FillWeighted( pass, w, ev -> muons.at(jmu).eta );
        muonPhi   -> FillWeighted( pass, w, ev -> muons.at(jmu).phi );
        muonIso   -> FillWeighted( pass, w, offlineiso              );
        muonIso04 -> FillWeighted( pass, w, offlineiso04            );
        nvtx      -> FillWeighted( pass, w, ev -> nVtx              );

        if (isBarrel){
          muonPt_barrel    -> FillWeighted( pass, w, ev -> muons.at(jmu).pt  );
          muonIso_barrel   -> FillWeighted( pass, w, offlineiso              );
          nvtx_barrel      -> FillWeighted( pass, w, ev -> nVtx              );
          ecalIso_barrel   -> Fill( hltEcalIso );
          hcalIso_barrel   -> Fill( hltHcalIso );
        }
        else{
          muonPt_endcap    -> FillWeighted( pass, w, ev -> muons.at(jmu).pt  );
          muonIso_endcap   -> FillWeighted( pass, w, offlineiso              );
          nvtx_endcap      -> FillWeighted( pass, w, ev -> nVtx              );
          ecalIso_endcap   -> Fill( hltEcalIso );
          hcalIso_endcap   -> Fill( hltHcalIso );
        }
        
        if (hltEcalIso < cut_ecal && 
            hltHcalIso < cut_hcal &&
            hltTrkIso  < cut_trk 
            ) num_pt -> Fill( ev -> muons.at(jmu).pt);  
        
	    for (int ih=0; ih < nhist_iso; ih++) {
          if (hltEcalIso < ih*0.01 && foundL3) passECal = true;
          if (hltHcalIso < ih*0.01 && foundL3) passHCal = true;
          if (hltTrkIso  < ih*0.01 && foundL3) passTrk  = true;
          if (isBarrel){
            h_ecal_iso_barrel[ih] -> FillWeighted( passECal, w, ev -> muons.at(jmu).pt) ;
            if (passECalNominal) h_hcal_iso_barrel[ih] -> FillWeighted( passHCal, w, ev -> muons.at(jmu).pt) ;
            if (passECalNominal && passHCalNominal) h_trk_iso_barrel[ih] -> FillWeighted( passTrk , w, ev -> muons.at(jmu).pt) ;
          }
          else{
            h_ecal_iso_endcap[ih] -> FillWeighted( passECal                    , w, ev -> muons.at(jmu).pt) ;
            if (passECalNominal) h_hcal_iso_endcap[ih] -> FillWeighted( passHCal, w, ev -> muons.at(jmu).pt) ;
            if (passECalNominal && passHCalNominal) h_trk_iso_endcap[ih] -> FillWeighted( passTrk , w, ev -> muons.at(jmu).pt) ;
          }
          passECal = false;
          passHCal = false;
          passTrk  = false;
        }
        
      }
      
    }



  }
  



  outfile         -> cd();
  muonPt          -> Write();
  muonPt_barrel   -> Write();
  muonPt_endcap   -> Write();
  tagMuonPt       -> Write();
  
  den_pt          -> Write();
  num_pt          -> Write();
  
  muonEta         -> Write();
  muonPhi         -> Write();
  
  nvtx_event      -> Write();
  nvtx            -> Write();
  nvtx_barrel     -> Write();
  nvtx_endcap     -> Write();
  
  muonIso         -> Write();
  muonIso_barrel  -> Write();
  muonIso_endcap  -> Write();
  muonIso04       -> Write();

  ecalIso_barrel  -> Write();
  hcalIso_barrel  -> Write();
  ecalIso_endcap  -> Write();
  hcalIso_endcap  -> Write();

  dimuon_mass     -> Write();

  for (int i=0; i< nhist_iso; i++) {
    h_ecal_iso_barrel[i] -> Write();
    h_hcal_iso_barrel[i] -> Write();
    h_trk_iso_barrel[i]  -> Write();
    h_ecal_iso_endcap[i] -> Write();
    h_hcal_iso_endcap[i] -> Write();
    h_trk_iso_endcap[i]  -> Write();
  }


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
  
  if (!( mu.pt         > 22  )) return false; 
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
  
  if (!( mu.pt         > cut_pt_offline  )) return false; 
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


