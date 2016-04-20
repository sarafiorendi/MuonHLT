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


bool  matchMuon  (HLTMuonCand, std::vector<HLTObjCand>, std::string);
void  dorate     (TFile*, TFile*);

// old wp
// EB 0.11 EE 0.08 HB 0.21 HE 0.22 TRK 0.09
// new loose wp
// EB 0.08 EE 0.06 HB 0.13 HE 0.13 TRK 0.08

// define pt threshold
float cut_pt  = 24; 
// define effective areas
// 2016 menu - ECal Tune
float a_ecal_barrel   = 0.135;   
float a_ecal_endcap   = 0.080;
float a_hcal_barrel   = 0.110;
float a_hcal_endcap   = 0.163;

// define isolation thresholds
float cut_ecal_barrel = 0;
float cut_hcal_barrel = 0;
float cut_ecal_endcap = 0;
float cut_hcal_endcap = 0;
float cut_trk         = 0;

// ********************************************************
int   which = 0; // 0 = 2015, 1 = loose 2016
// ********************************************************


void readNtuplesRate(){

  if (which == 0){
	cut_ecal_barrel = 0.11    ;
	cut_hcal_barrel = 0.21    ;
	cut_ecal_endcap = 0.08    ;
	cut_hcal_endcap = 0.22    ;
	cut_trk         = 0.09    ;
  }  
  else if (which == 1){
	cut_ecal_barrel = 0.08   ;
	cut_hcal_barrel = 0.13   ;
	cut_ecal_endcap = 0.06   ;
	cut_hcal_endcap = 0.14   ;
	cut_trk         = 0.08   ;
  }  

  std::string thestr = "_pu30_pt24_wptest";

  TFile* inputfile0 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD20to30.root","READ");
  TFile* outfile0 = TFile::Open(Form("rates/abs_rate_QCD20to30%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile0, outfile0);
  
  TFile* inputfile1 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD30to50.root","READ");
  TFile* outfile1 = TFile::Open(Form("rates/abs_rate_QCD30to50%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile1, outfile1);

  TFile* inputfile2 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD50to80.root","READ");
  TFile* outfile2 = TFile::Open(Form("rates/abs_rate_QCD50to80%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile2, outfile2);

  TFile* inputfile3 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_QCD80to120.root","READ");
  TFile* outfile3 = TFile::Open(Form("rates/abs_rate_QCD80to120%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile3, outfile3);

  TFile* inputfile4 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_WToMuNu.root","READ");
  TFile* outfile4 = TFile::Open(Form("rates/abs_rate_WMuNu%s.root", thestr.c_str()),"RECREATE");
  dorate(inputfile4, outfile4);

  TFile* inputfile5 = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/76/CMSSW_7_6_3/src/HLTrigger/Configuration/test/muonNtuple_DYToLL.root","READ");
  TFile* outfile5 = TFile::Open(Form("rates/abs_rate_DYLL%s.root", thestr.c_str()),"RECREATE");
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

  TH1F* dimuon_mass         = new TH1F("dimuon_mass"         ,"dimuon_mass"          , 1500,  0, 150);
  TH1F* muon_size           = new TH1F("muon_size"           ,"muon_size"            ,   10,  0, 10 );
  TH1F* tagiso              = new TH1F("tagiso"              ,"tagiso"               ,  100,  0, 1  );
  TH1F* tagMuonPt           = new TH1F("tagMuonPt"           ,"tagMuonPt"            ,  150,  0, 150);
  TH1F* n_passing_events    = new TH1F("n_passing_events"    ,"n_passing_events"     ,    5,  0,   5);
  TH1F* nvtx_event          = new TH1F("nvtx_event"          ,"nvtx_event"           ,   60,  0,  60);
  TH1F* nvtx_barrel         = new TH1F("nvtx_barrel"         ,"nvtx_barrel"          ,   60,  0,  60);
  TH1F* nvtx_endcap         = new TH1F("nvtx_endcap"         ,"nvtx_endcap"          ,   60,  0,  60);
  TH1F* nvtx_barrel_hcal    = new TH1F("nvtx_barrel_hcal"    ,"nvtx_barrel_hcal"     ,   60,  0,  60);
  TH1F* nvtx_endcap_hcal    = new TH1F("nvtx_endcap_hcal"    ,"nvtx_endcap_hcal"     ,   60,  0,  60);
  TH1F* nvtx_barrel_trk     = new TH1F("nvtx_barrel_trk"     ,"nvtx_barrel_trk"      ,   60,  0,  60);
  TH1F* nvtx_endcap_trk     = new TH1F("nvtx_endcap_trk"     ,"nvtx_endcap_trk"      ,   60,  0,  60);
  TH1F* wp_passing          = new TH1F("wp_passing"          ,"wp_passing"           ,   60,  0,  60);

  TH1F* ecalIso_barrel      = new TH1F("ecalIso_barrel"      ,"ecalIso_barrel"       ,   50,  0,  0.5);
  TH1F* hcalIso_barrel      = new TH1F("hcalIso_barrel"      ,"hcalIso_barrel"       ,   50,  0,  0.5);
  TH1F* ecalIso_endcap      = new TH1F("ecalIso_endcap"      ,"ecalIso_endcap"       ,   50,  0,  0.5);
  TH1F* hcalIso_endcap      = new TH1F("hcalIso_endcap"      ,"hcalIso_endcap"       ,   50,  0,  0.5);

  TH1F* ecalDep     = new TH1F("ecalDep"      ,"ecalDep"       ,   500,  0,  50);
  TH1F* hcalDep     = new TH1F("hcalDep"      ,"hcalDep"       ,   500,  0,  50);
  TH1F* trkPt       = new TH1F("trkPt"        ,"trkPt"         ,   500,  0, 100);
  TH1F* rho         = new TH1F("rho"          ,"rho"           ,   500,  0, 100);

  const int nhist_iso = 50;
  TH1F *h_ecal_iso_barrel    [nhist_iso];
  TH1F *h_hcal_iso_barrel    [nhist_iso];
  TH1F *h_trk_iso_barrel     [nhist_iso];

  TH1F *h_ecal_iso_endcap    [nhist_iso];
  TH1F *h_hcal_iso_endcap    [nhist_iso];
  TH1F *h_trk_iso_endcap     [nhist_iso];

  for (int i=0; i< nhist_iso; i++) {
    h_ecal_iso_barrel[i] = new TH1F(Form("HMuonPt_barrel_ecal_iso%d",i), "", 60, 0, 60);
    h_hcal_iso_barrel[i] = new TH1F(Form("HMuonPt_barrel_hcal_iso%d",i), "", 60, 0, 60);
    h_trk_iso_barrel[i]  = new TH1F(Form("HMuonPt_barrel_trk_iso%d" ,i), "", 60, 0, 60);
    h_ecal_iso_endcap[i] = new TH1F(Form("HMuonPt_endcap_ecal_iso%d",i), "", 60, 0, 60);
    h_hcal_iso_endcap[i] = new TH1F(Form("HMuonPt_endcap_hcal_iso%d",i), "", 60, 0, 60);
    h_trk_iso_endcap[i]  = new TH1F(Form("HMuonPt_endcap_trk_iso%d" ,i), "", 60, 0, 60);
  }

  
  float a_ecal  , a_hcal  ;
  float cut_ecal, cut_hcal;
  bool passEcalIso_b[nhist_iso],  passHcalIso_b[nhist_iso], passTrkIso_b[nhist_iso] ;
  bool passEcalIso_e[nhist_iso],  passHcalIso_e[nhist_iso], passTrkIso_e[nhist_iso] ;

  MuonEvent* ev = new MuonEvent();

  int count = 0;
  tree -> SetBranchAddress("event", &ev);

  bool isBarrel = false;
  int nentries = tree->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  bool passbarrel      = false;
  bool passendcap      = false;
  bool passbarrel_hcal = false;
  bool passendcap_hcal = false;
  bool passbarrel_trk  = false;
  bool passendcap_trk  = false;
  bool passwp          = false;
  bool passref         = false;

  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {

    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    
    passbarrel      = false;
    passendcap      = false;
    passbarrel_hcal = false;
    passendcap_hcal = false;
    passbarrel_trk  = false;
    passendcap_trk  = false;
    passwp          = false; 
    for (int ii=0; ii < nhist_iso ; ii++) {  
      passEcalIso_b    [ii] = false;  
      passHcalIso_b    [ii] = false;  
      passTrkIso_b     [ii] = false;  
      passEcalIso_e    [ii] = false;  
      passHcalIso_e    [ii] = false;  
      passTrkIso_e     [ii] = false;  
    }    

//     if (!ev-> hlt.find("HLT_Mu20_v2")) continue;
    if ( ev -> trueNI < 28 || ev -> trueNI > 32 ) continue;
    count++;
    
    passref = false;

    unsigned int nhltmuons = ev->hltmuons.size();
    muon_size       -> Fill( nhltmuons );
    for (int imu = 0; imu < nhltmuons; imu++){
      
      isBarrel = false;
      if (! matchMuon(ev ->hltmuons.at(imu), ev -> hlt.objects, "hltL3fL1sMu16L1f0L2f10QL3Filtered20Q::REHLT")) continue;

      HLTMuonCand theL3 = ev -> hltmuons.at(imu);
      tagMuonPt -> Fill(theL3.pt);
  	  if (!passref) {
  	    n_passing_events      -> Fill( 2 );
        passref = true;
      }

      if (! (theL3.pt > cut_pt) )  continue;
	  nvtx_event      -> Fill( 3         );

      if (fabs(theL3.eta) < 1.479){
		a_ecal   = a_ecal_barrel;
		a_hcal   = a_hcal_barrel;
		cut_ecal = cut_ecal_barrel;
		cut_hcal = cut_hcal_barrel;
      }
      else{
		a_ecal   = a_ecal_endcap;
		a_hcal   = a_hcal_endcap;
		cut_ecal = cut_ecal_endcap;
		cut_hcal = cut_hcal_endcap;
      }

	  double hltEcalIso = max(0., (theL3.ecalDep05 - a_ecal * ev -> hlt.rho ) )/ theL3.pt;
	  double hltHcalIso = max(0., (theL3.hcalDep1  - a_hcal * ev -> hlt.rho ) )/ theL3.pt;
	  double hltTrkIso  = (theL3.trkDep) / theL3.pt;

      if (fabs(theL3.eta) < 1.479){
        ecalIso_barrel   -> Fill( hltEcalIso );
        hcalIso_barrel   -> Fill( hltHcalIso );
        isBarrel = true;
        if (passbarrel == false ) nvtx_barrel -> Fill( 3  );
        passbarrel = true;
        if (passbarrel_hcal == false && hltEcalIso < cut_ecal) nvtx_barrel_hcal -> Fill( 3  );
        passbarrel_hcal = true;
        if (passbarrel_trk  == false && hltEcalIso < cut_ecal && hltHcalIso < cut_hcal) nvtx_barrel_trk -> Fill( 3  );
        passbarrel_trk  = true;
      }
      else{
        ecalIso_endcap   -> Fill( hltEcalIso );
        hcalIso_endcap   -> Fill( hltHcalIso );
        if (passendcap == false ) nvtx_endcap -> Fill( 3  );
        passendcap = true;
        if (passendcap_hcal == false && hltEcalIso < cut_ecal) nvtx_endcap_hcal -> Fill( 3  );
        passendcap_hcal = true;
        if (passendcap_trk  == false && hltEcalIso < cut_ecal && hltHcalIso < cut_hcal) nvtx_endcap_trk -> Fill( 3  );
        passendcap_trk  = true;
      }


      ecalDep -> Fill(theL3.ecalDep05);
      hcalDep -> Fill(theL3.hcalDep1 );
      trkPt   -> Fill(theL3.trkpt    );
      rho     -> Fill(ev -> hlt.rho  );

      if (hltEcalIso <= cut_ecal && hltHcalIso < cut_hcal && hltTrkIso < cut_trk && passwp==false && theL3.pt > cut_pt ) {
        wp_passing -> Fill(3);
	    passwp = true;
	  }
	  for (int ih=0; ih < nhist_iso; ih++) {
		if (isBarrel){
		  if (hltEcalIso < ih*0.01 && passEcalIso_b[ih]==false ) {
			h_ecal_iso_barrel[ih] -> Fill(3);
			passEcalIso_b[ih] = true;
		  }
		  if (hltEcalIso < cut_ecal && hltHcalIso < ih*0.01 && passHcalIso_b[ih]==false ) {
			h_hcal_iso_barrel[ih] -> Fill(3);
			passHcalIso_b[ih] = true;
		  }
		  if (hltEcalIso < cut_ecal && hltHcalIso < cut_hcal && hltTrkIso  < ih*0.01 && passTrkIso_b[ih]==false ) {
			h_trk_iso_barrel[ih] -> Fill(3);
			passTrkIso_b[ih] = true;
		  }
		}
		else{
		  if (hltEcalIso < ih*0.01 && passEcalIso_e[ih]==false ) {
			h_ecal_iso_endcap[ih] -> Fill(3);
			passEcalIso_e[ih] = true;
		  }
		  if (hltEcalIso < cut_ecal && hltHcalIso < ih*0.01 && passHcalIso_e[ih]==false ) {
			h_hcal_iso_endcap[ih] -> Fill(3);
			passHcalIso_e[ih] = true;
		  }
		  if (hltEcalIso < cut_ecal && hltHcalIso < cut_hcal && hltTrkIso  < ih*0.01 && passTrkIso_e[ih]==false ) {
			h_trk_iso_endcap[ih] -> Fill(3);
			passTrkIso_e[ih] = true;
		  }
		}
	  }
   
    }




  }
  std::cout << "events in PU range: " << count << std::endl;
  


  outfile            -> cd();
  
  n_passing_events   -> Write();
  nvtx_event         -> Write();
  nvtx_barrel        -> Write();
  nvtx_endcap        -> Write();
  nvtx_barrel_hcal   -> Write();
  nvtx_endcap_hcal   -> Write();
  nvtx_barrel_trk    -> Write();
  nvtx_endcap_trk    -> Write();

  muon_size          -> Write();
  dimuon_mass        -> Write();
  tagiso             -> Write();
  wp_passing         -> Write();
  tagMuonPt          -> Write();

  ecalIso_barrel     -> Write();
  hcalIso_barrel     -> Write();
  ecalIso_endcap     -> Write();
  hcalIso_endcap     -> Write();
  
  ecalDep -> Write();
  hcalDep -> Write();
  trkPt   -> Write();
  rho     -> Write();
 
  
  for (int i=0; i< nhist_iso; i++) {
    h_ecal_iso_barrel[i] -> Write();
    h_hcal_iso_barrel[i] -> Write();
    h_trk_iso_barrel[i]  -> Write();
    h_ecal_iso_endcap[i] -> Write();
    h_hcal_iso_endcap[i] -> Write();
    h_trk_iso_endcap[i]  -> Write();
  }

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




