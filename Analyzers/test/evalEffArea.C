#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <string.h>

#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLT/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

std::string isofilter   = "hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09::REHLT";


bool                selectTagMuon  (MuonCand, TH1F* );
bool                selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool                matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
HLTMuonCand         matchL3        (MuonCand, std::vector<HLTMuonCand>);
bool                matchGENMuon   (MuonCand, std::vector<GenParticleCand>);
void                dofit          ();
pair<float, float>  doReallyFit    (TH1F* , TH1F* , TH1F*, int , int, std::string, std::string );
void                calcEffArea    (pair<float, float>, pair<float, float>, std::string );

const int    nbins    = 80;
const int    minbin   =  0;
const int    maxbin   = 80;
const double muonmass =  0.10565837;

TH1F * tagiso                  = new TH1F("tagiso"                   ,"tagiso"                        ,  100,      0,    1  );

TH1F * HMeanRhoVsNVtx          = new TH1F("HMeanRhoVsNVtx"           ,"MeanRhoVsNVtx"                 ,nbins, minbin, maxbin);
TH1F * HMeanRhoECALVsNVtx      = new TH1F("HMeanRhoECALVsNVtx"       ,"MeanRhoECALVsNVtx"             ,nbins, minbin, maxbin);
TH1F * HMeanRhoHCALVsNVtx      = new TH1F("HMeanRhoHCALVsNVtx"       ,"MeanRhoHCALVsNVtx"             ,nbins, minbin, maxbin);
TH1F * HRhoVsNVtx              = new TH1F("HRhoVsNVtx"               ,"RhoVsNVtx"                     ,nbins, minbin, maxbin);
TH1F * HRhoECALVsNVtx          = new TH1F("HRhoECALVsNVtx"           ,"RhoECALVsNVtx"                 ,nbins, minbin, maxbin);
TH1F * HRhoHCALVsNVtx          = new TH1F("HRhoHCALVsNVtx"           ,"RhoHCALVsNVtx"                 ,nbins, minbin, maxbin);


TH1F * HNVtx                   = new TH1F("HNVtx"                    ,"NVtx"                          ,nbins, minbin, maxbin);
 
TH1F * HMeanECalIsoVsNVtx      = new TH1F("HMeanECalIsoVsNVtx"       ,"MeanECalIsoVsNVtx"             ,nbins, minbin, maxbin);
TH1F * HMeanHCalIsoVsNVtx      = new TH1F("HMeanHCalIsoVsNVtx"       ,"MeanHCalIsoVsNVtx"             ,nbins, minbin, maxbin);
TH1F * HECalIsoVsNVtx          = new TH1F("HECalIsoVsNVtx"           ,"HECalIsoVsNVtx"                ,nbins, minbin, maxbin);
TH1F * HHCalIsoVsNVtx          = new TH1F("HHCalIsoVsNVtx"           ,"HHCalIsoVsNVtx"                ,nbins, minbin, maxbin);
TH1F * HNVtxForIso             = new TH1F("HNVtxForIso"              ,"NVtx for mean iso"             ,nbins, minbin, maxbin);
 
TH1F * HMeanECalIsoVsNVtx_eta0 = new TH1F("HMeanECalIsoVsNVtx_eta0"  ,"MeanECalIsoVsNVtx_eta0"        ,nbins, minbin, maxbin);
TH1F * HECalIsoVsNVtx_eta0     = new TH1F("HECalIsoVsNVtx_eta0"      ,"ECalIsoVsNVtx_eta0"            ,nbins, minbin, maxbin);
TH1F * HMeanHCalIsoVsNVtx_eta0 = new TH1F("HMeanHCalIsoVsNVtx_eta0"  ,"MeanHCalIsoVsNVtx_eta0"        ,nbins, minbin, maxbin);
TH1F * HHCalIsoVsNVtx_eta0     = new TH1F("HHCalIsoVsNVtx_eta0"      ,"HCalIsoVsNVtx_eta0"            ,nbins, minbin, maxbin);
TH1F * HNVtxForIso_eta0        = new TH1F("HNVtxForIso_eta0"         ,"NVtx for mean iso_eta0"        ,nbins, minbin, maxbin);
 
TH1F * HMeanECalIsoVsNVtx_eta1 = new TH1F("HMeanECalIsoVsNVtx_eta1"  ,"MeanECalIsoVsNVtx_eta1"        ,nbins, minbin, maxbin);
TH1F * HECalIsoVsNVtx_eta1     = new TH1F("HECalIsoVsNVtx_eta1"      ,"ECalIsoVsNVtx_eta1"            ,nbins, minbin, maxbin);
TH1F * HMeanHCalIsoVsNVtx_eta1 = new TH1F("HMeanHCalIsoVsNVtx_eta1"  ,"MeanHCalIsoVsNVtx_eta1"        ,nbins, minbin, maxbin);
TH1F * HHCalIsoVsNVtx_eta1     = new TH1F("HHCalIsoVsNVtx_eta1"      ,"HCalIsoVsNVtx_eta1"            ,nbins, minbin, maxbin);
TH1F * HNVtxForIso_eta1        = new TH1F("HNVtxForIso_eta1"         ,"NVtx for mean iso_eta1"        ,nbins, minbin, maxbin);

TH2F * HNVtxECalDep_barrel     = new TH2F("HNVtxECalDep_barrel"      ,"NVtx vs ECal dep bar"          ,nbins, minbin, maxbin, 200, -10, 10);
TH2F * HNVtxECalDep_endcap     = new TH2F("HNVtxECalDep_endcap"      ,"NVtx vs ECal dep end"          ,nbins, minbin, maxbin, 200, -10, 10);
TH2F * HNVtxHCalDep_barrel     = new TH2F("HNVtxHCalDep_barrel"      ,"NVtx vs HCal dep bar"          ,nbins, minbin, maxbin, 200, -10, 10);
TH2F * HNVtxHCalDep_endcap     = new TH2F("HNVtxHCalDep_endcap"      ,"NVtx vs HCal dep end"          ,nbins, minbin, maxbin, 200, -10, 10);

TH2F * HRhoECalDep_barrel      = new TH2F("HRhoECalDep_barrel"       ,"Rho vs ECal dep bar"           ,nbins, minbin, maxbin, 200, -10, 10);
TH2F * HRhoECalDep_endcap      = new TH2F("HRhoECalDep_endcap"       ,"Rho vs ECal dep bar"           ,nbins, minbin, maxbin, 200, -10, 10);
TH2F * HRhoHCalDep_barrel      = new TH2F("HRhoHCalDep_barrel"       ,"Rho vs HCal dep bar"           ,nbins, minbin, maxbin, 200, -10, 10);
TH2F * HRhoHCalDep_endcap      = new TH2F("HRhoHCalDep_endcap"       ,"Rho vs HCal dep bar"           ,nbins, minbin, maxbin, 200, -10, 10);

TH2F * HRhoECalIso_barrel      = new TH2F("HRhoECalIso_barrel"       ,"Rho vs ECal isolation"         ,nbins, minbin, maxbin, 1000, 0, 10);
TH2F * HRhoECalIso_endcap      = new TH2F("HRhoECalIso_endcap"       ,"Rho vs ECal isolation"         ,nbins, minbin, maxbin, 1000, 0, 10);

TProfile* hprof                = new TProfile("hprof","Profile of ECalDep versus nvtx", nbins, minbin, maxbin, -10, 10);
TH1F* tagMuonPt                = new TH1F("tagMuonPt"              ,"tagMuonPt"            ,  150,  0, 150);

TH1F * HECalDep         = new TH1F("HECalDep"           ,"HECalDep"                ,200, -1, 49);
TH1F * HHCalDep         = new TH1F("HHCalDep"           ,"HHCalDep"                ,200, -1, 49);


void evalEffArea(){


  TChain* tree = new TChain("muonNtuples/muonTree");
  for (int i = 1; i < 10; i++){
    tree -> Add(Form("root://eoscms//eos/cms/store/group/phys_muon/fiorendi/13TeV/MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_isolation2017_eff_83Xsample_fillOffline/170507_161411/0000/muonNtuple_%d.root/muonNtuples/muonTree", i));
  }

//   TFile* inputfile = TFile::Open("/afs/cern.ch/work/f/fiorendi/private/MuonHLTRegMuVtx/CMSSW_8_0_0/src/HLTrigger/Configuration/test/muonNtuple_DYToLL_ECalTune.root","READ");
//   std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TFile* outfile = TFile::Open("effArea_83X_pt24.root","RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  TH1F* dimuon_mass              = new TH1F("dimuon_mass"              ,"dimuon_mass"                   , 1500,  0, 150);
         

  MuonEvent* ev = new MuonEvent();
  tree -> SetBranchStatus("*",1);
  tree -> SetBranchAddress( "event", &ev);

  int nentries = tree->GetEntries();
//   nentries = 500000;
  std::cout << "Number of entries = " << nentries << std::endl;


  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    
    unsigned int nmuons = ev->muons.size();
    if (nmuons < 2) continue;

    unsigned int nhltmuons = ev->hltmuons.size();
    std::cout << nhltmuons << std::endl;
    
//     for ( std::vector<std::string>::const_iterator it = ev-> hlt.triggers.begin(); it != ev-> hlt.triggers.end(); ++it ) {
//       std::cout << *it << std::endl;
// 	}
    
    if (!ev-> hlt.find("HLT_IsoMu24_v4")) continue;

    for (int imu = 0; imu < nmuons; imu++){
      
      // select the tag muon        
      if (! selectTagMuon(ev -> muons.at(imu), tagiso)) continue;
      if (! matchMuon(ev -> muons.at(imu), ev -> hlt.objects, isofilter )) continue;
//       if (! matchGENMuon (ev -> muons.at(imu), ev -> genParticles))           continue;
      tagMuonPt -> Fill(ev -> muons.at(imu).pt);
      
      for (int jmu = 0; jmu < nmuons; jmu++){
        // select the probe muon
        if (! selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
//         if (! matchGENMuon   (ev -> muons.at(jmu), ev -> genParticles))               continue;
        HRhoVsNVtx     -> Fill(ev -> nVtx, ev -> hlt.rho);
        HRhoECALVsNVtx -> Fill(ev -> nVtx, ev -> hlt.rho_ecal);
        HRhoHCALVsNVtx -> Fill(ev -> nVtx, ev -> hlt.rho_hcal);
        HNVtx      -> Fill(ev -> nVtx);
        
//         muonPt_den  -> Fill( ev -> muons.at(jmu).pt );

        if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, "hltL3fL1sMu22L1f0L2f10QL3Filtered24Q::REHLT")){
          HNVtxForIso-> Fill(ev -> nVtx);
          HLTMuonCand theL3 = matchL3(ev -> muons.at(jmu), ev -> hltmuons);
          if (theL3.pt == -1000) {std::cout << "no L3 found" << std::endl; continue;}
          if (theL3.pt < 24    ) continue;
          
          HECalDep -> Fill(theL3.ecalDep);
          HHCalDep -> Fill(theL3.hcalDep);
		  HECalIsoVsNVtx -> Fill(ev -> nVtx, theL3.ecalDep);
		  HHCalIsoVsNVtx -> Fill(ev -> nVtx, theL3.hcalDep);

		  if (fabs(ev -> muons.at(jmu).eta) <= 1.479) {
		    HECalIsoVsNVtx_eta0 -> Fill( ev -> nVtx,    theL3.ecalDep );
		    HHCalIsoVsNVtx_eta0 -> Fill( ev -> nVtx,    theL3.hcalDep  );
		    HNVtxForIso_eta0    -> Fill( ev -> nVtx                     );
		    HNVtxECalDep_barrel -> Fill( ev -> nVtx,    theL3.ecalDep );
		    HNVtxHCalDep_barrel -> Fill( ev -> nVtx,    theL3.hcalDep  );
		    hprof               -> Fill( ev -> nVtx,    theL3.ecalDep );
		    HRhoECalDep_barrel  -> Fill( ev -> hlt.rho_ecal, theL3.ecalDep );
		    HRhoECalIso_barrel  -> Fill( ev -> hlt.rho_ecal, theL3.ecalDep/theL3.pt );
		    HRhoHCalDep_barrel  -> Fill( ev -> hlt.rho, theL3.hcalDep  );
		  }
		  else { 
			HECalIsoVsNVtx_eta1 -> Fill( ev -> nVtx, theL3.ecalDep    );
			HHCalIsoVsNVtx_eta1 -> Fill( ev -> nVtx, theL3.hcalDep     );
			HNVtxForIso_eta1    -> Fill( ev -> nVtx                     );
		    HNVtxECalDep_endcap -> Fill( ev -> nVtx, theL3.ecalDep    );
		    HNVtxHCalDep_endcap -> Fill( ev -> nVtx, theL3.hcalDep     );
		    HRhoECalDep_endcap  -> Fill( ev -> hlt.rho, theL3.ecalDep );
		    HRhoHCalDep_endcap  -> Fill( ev -> hlt.rho, theL3.hcalDep  );
		    HRhoECalIso_endcap  -> Fill( ev -> hlt.rho, theL3.ecalDep/theL3.pt );
		  }

        }
      } // end loop on second muon
    } // end loop on first muon
  } // end loop on events
  

  dofit();

  
  outfile                 -> cd();
  
  tagMuonPt               -> Write();
  tagiso                  -> Write();
  dimuon_mass             -> Write();
  HECalDep                -> Write();
  HHCalDep                -> Write();
  
  hprof                   -> Write();
  HNVtxECalDep_barrel     -> Write();
  HNVtxECalDep_endcap     -> Write();
  HNVtxHCalDep_barrel     -> Write();
  HNVtxHCalDep_endcap     -> Write();
  
  HRhoECalDep_barrel      -> Write();
  HRhoECalDep_endcap      -> Write();
  HRhoHCalDep_barrel      -> Write();
  HRhoHCalDep_endcap      -> Write();

  HRhoECalIso_barrel      -> Write();
  HRhoECalIso_endcap      -> Write();
    
  HMeanRhoVsNVtx          -> Write();
  HMeanRhoECALVsNVtx      -> Write();
  HMeanRhoHCALVsNVtx      -> Write();
  HRhoVsNVtx              -> Write();
  HRhoECALVsNVtx          -> Write();
  HRhoHCALVsNVtx          -> Write();
  HNVtx                   -> Write();

  HMeanECalIsoVsNVtx      -> Write();
  HECalIsoVsNVtx          -> Write();
  HMeanHCalIsoVsNVtx      -> Write();
  HHCalIsoVsNVtx          -> Write();
  HNVtxForIso             -> Write();

  HMeanECalIsoVsNVtx_eta0 -> Write();
  HECalIsoVsNVtx_eta0     -> Write();
  HMeanHCalIsoVsNVtx_eta0 -> Write();
  HHCalIsoVsNVtx_eta0     -> Write();
  HNVtxForIso_eta0        -> Write();

  HMeanECalIsoVsNVtx_eta1 -> Write();
  HECalIsoVsNVtx_eta1     -> Write();
  HMeanHCalIsoVsNVtx_eta1 -> Write();
  HHCalIsoVsNVtx_eta1     -> Write();
  HNVtxForIso_eta1        -> Write();

  outfile                 -> Close();    
//   inputfile               -> Close();    
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
  
  if (!( mu.pt         > 26  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0.,
                       mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  tagiso -> Fill(offlineiso04);
  if (!(offlineiso04   < 0.15)) return false;
  
  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt         > 26  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  if (mu.charge * tagMu.charge > 0) return false;
  
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
  theL3.pt = -1000;
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



void dofit(){

  int minF = 20 ;
  int maxF = 60 ;
  
  gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);

  pair<float, float> a_pair        = doReallyFit( HRhoVsNVtx    , HNVtx, HMeanRhoVsNVtx    , minF, maxF, "NVtx", "<#rho>" );
  pair<float, float> a_ecal_pair   = doReallyFit( HRhoECALVsNVtx, HNVtx, HMeanRhoECALVsNVtx, minF, maxF, "NVtx", "<#rho>" );
  pair<float, float> a_hcal_pair   = doReallyFit( HRhoHCALVsNVtx, HNVtx, HMeanRhoHCALVsNVtx, minF, maxF, "NVtx", "<#rho>" );
  
  pair<float, float> k_ecal        = doReallyFit( HECalIsoVsNVtx     , HNVtxForIso     , HMeanECalIsoVsNVtx     , minF, maxF, "NVtx", "<ecal iso>" );
  pair<float, float> k_ecal_barrel = doReallyFit( HECalIsoVsNVtx_eta0, HNVtxForIso_eta0, HMeanECalIsoVsNVtx_eta0, minF, maxF, "NVtx", "<ecal iso>" );
  pair<float, float> k_ecal_endcap = doReallyFit( HECalIsoVsNVtx_eta1, HNVtxForIso_eta1, HMeanECalIsoVsNVtx_eta1, minF, maxF, "NVtx", "<ecal iso>" );

  pair<float, float> k_hcal        = doReallyFit( HHCalIsoVsNVtx     , HNVtxForIso     , HMeanHCalIsoVsNVtx     , minF, maxF, "NVtx", "<hcal iso>" );
  pair<float, float> k_hcal_barrel = doReallyFit( HHCalIsoVsNVtx_eta0, HNVtxForIso_eta0, HMeanHCalIsoVsNVtx_eta0, minF, maxF, "NVtx", "<hcal iso>" );
  pair<float, float> k_hcal_endcap = doReallyFit( HHCalIsoVsNVtx_eta1, HNVtxForIso_eta1, HMeanHCalIsoVsNVtx_eta1, minF, maxF, "NVtx", "<hcal iso>" );

  HMeanECalIsoVsNVtx -> Draw();
  
  calcEffArea(a_pair, k_ecal       , "ecal: whole range: " );
  calcEffArea(a_pair, k_ecal_barrel, "ecal: barrel     : " );
  calcEffArea(a_pair, k_ecal_endcap, "ecal: endcap     : " );
  calcEffArea(a_pair, k_hcal       , "hcal: whole range: " );
  calcEffArea(a_pair, k_hcal_barrel, "hcal: barrel     : " );
  calcEffArea(a_pair, k_hcal_endcap, "hcal: endcap     : " );

  calcEffArea(a_ecal_pair, k_ecal       , "ecal rho_ecal: whole range: " );
  calcEffArea(a_ecal_pair, k_ecal_barrel, "ecal rho_ecal: barrel     : " );
  calcEffArea(a_ecal_pair, k_ecal_endcap, "ecal rho_ecal: endcap     : " );
  calcEffArea(a_hcal_pair, k_hcal       , "hcal rho_hcal: whole range: " );
  calcEffArea(a_hcal_pair, k_hcal_barrel, "hcal rho_hcal: barrel     : " );
  calcEffArea(a_hcal_pair, k_hcal_endcap, "hcal rho_hcal: endcap     : " );


}



void calcEffArea(pair<float, float> a, pair<float, float> k, std::string str){
  
  float value  = k.first / a.first;
  float error  = pow(a.second/a.first,2) + pow(k.second/k.first,2);  
  error = sqrt(error);
  error = error*value;
  
  std::cout << str.c_str() << std::fixed << std::setprecision(5) << value <<" +/- " << error << std::endl;
  
}



pair<float, float> doReallyFit(TH1F* h_num, TH1F* h_den, TH1F* h_ratio, int minF, int maxF, std::string xtitle, std::string ytitle){

  pair<float, float> result;

  TF1 *f_pol = new TF1("f_pol", "pol1", minF, maxF);

  h_ratio -> Divide( h_num, h_den );
  h_ratio -> Fit   ( f_pol, "R"   );
  result.first   = f_pol -> GetParameter(1);
  result.second  = f_pol -> GetParError(1);

  h_ratio -> SetTitle("");
  h_ratio -> GetXaxis() -> SetTitle(xtitle.c_str());
  h_ratio -> GetYaxis() -> SetTitle(ytitle.c_str());
  h_ratio -> SetLineColor(kBlack);

  return result;
  
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



