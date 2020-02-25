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


bool matchToL1      (HLTMuonCand, L1MuonCand );
bool matchToL1      (HLTJetCand,  L1MuonCand );
void dorate         (TChain*, TFile* , TFile*);

// define pt threshold
const int n_lumi_columns = 1;

float cut_pt_l1[]   = {22,    18.,   14.,  12.,  10., 9,     8.,   7.   , 6.};//10.5, 14, 16.5, 18.5, 19.5, 20, 20.5, 21., 22}; 
float cut_eta_l1[]  = {2.5,   1.5,   1.5,  1.5,  1.5, 1.5,  1.5,   1.5 , 1.5};//10.5, 14, 16.5, 18.5, 19.5, 20, 20.5, 21., 22}; 
const int n_pt_l1  = 1;

float cut_pt_hlt[40];
const int nhist_pt  = 40;

float cut_dxy_s[40];
const int nhist_dxy  = 40;

int qualcut   = 11;
int mylumicol = 7;


void readNtuples_ROCPt_ZB(){

  for (int i=0; i < nhist_dxy; i++){
    cut_dxy_s[i] = 0.2*i - 0.2;
  }
  for (int i=0; i < nhist_pt; i++){
    cut_pt_hlt[i] = 0.5*i;
  }

  gROOT->ProcessLine("gErrorIgnoreLevel = kFatal;");

  TChain* tree0 = new TChain("muonNtuples/muonTree");
  TFile* tfile0  = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2018/EphemeralZeroBias1/crab_singlemu_data2018_rate_fill6944_EphZeroBias1_fromDAQ_emuL1_ignoreFilters/skim_hasL1Mu.root", "r");
  tree0 -> Add(tfile0->GetName());

//   TFile* tfile1  = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2017/EphemeralZeroBias2/crab_singlemu_nofilteriso_jet_rate_Run305178_305183_EphZeroBias2/skim_hasL1Mu.root", "r");
//   tree0 -> Add(tfile1->GetName());
  
  TFile* outfile0 = TFile::Open("plotPt_ZB_HLTMu.root","RECREATE");
  dorate(tree0, outfile0, tfile0);

  std::cout << "have you added all the nentries files??? " << std::endl;
  TH1F* nentries_nofilter_trueni_0 = (TH1F*)tfile0 -> Get("nentries_nofilter_trueni");
//   TH1F* nentries_nofilter_trueni_1 = (TH1F*)tfile1 -> Get("nentries_nofilter_trueni");
//   nentries_nofilter_trueni_0 -> Add(nentries_nofilter_trueni_1);
  outfile0 ->cd();
  nentries_nofilter_trueni_0 -> SetName("nentries_nofilter_trueni");
  nentries_nofilter_trueni_0 -> Write();

  outfile0            -> Close();  
  return;

}


void dorate( TChain* tree,  TFile* outfile, TFile* infile){

  TH1F* L3muon_size         = new TH1F("L3muon_size"         ,"L3muon_size"          ,   10,   0,  10 );
  TH1F* L1muon_size         = new TH1F("L1muon_size"         ,"L1muon_size"          ,   10,   0,  10 );
  TH1F* lumis               = new TH1F("lumis"               ,"lumis"                ,   200,  0,  200);

  TH1F *pt_hlt        [nhist_pt];
  TH1F *dxy_err_hlt   [nhist_dxy];
  TH1F *pt_dxy_hlt    [nhist_dxy*nhist_pt];
  for (int i=0; i < nhist_pt; i++) {
    pt_hlt[i]            = new TH1F(Form("pt_hlt_%d"  , i),     "", 600,  0,  60);
  }  
  for (int i=0; i < nhist_dxy; i++) {
    dxy_err_hlt[i]       = new TH1F(Form("dxy_err_hlt_%d"  , i),     "", 1000,  0,  500);
  }
  int k = 0;
  for (int i=0; i < nhist_pt; i++) {
    for (int j=0; j < nhist_dxy; j++) { 
      pt_dxy_hlt[k]       = new TH1F(Form("pt_dxy_hlt_%d"       , k),     "", 1000,  0,  500);
      k++;
    }
  }


  bool passpt[nhist_pt];
  bool passdxy[nhist_dxy];
  bool passcomb[nhist_dxy*nhist_pt];

  MuonEvent* ev = new MuonEvent();
  tree -> SetBranchStatus("*",1);
  tree -> SetBranchAddress( "event", &ev);

  int nentries = tree->GetEntries();
  std::cout << "Number of entries = "   << nentries  << std::endl;
      
  int count = 0;
  int thecut = 0;
  
  for (int ilumi = 0; ilumi < n_lumi_columns; ilumi++)
  {
    ilumi = mylumicol;
    for (Int_t eventNo=0; eventNo < nentries; eventNo++)
    {
  
      Int_t IgetEvent   = tree   -> GetEvent(eventNo);
      count++;
   std::cout << "ok 127" << std::endl;
      
      // reset all booleans
      for (int ibool = 0; ibool < nhist_pt ; ibool++){
        passpt[ibool]= false;
      }
      for (int ibool = 0; ibool < nhist_dxy ; ibool++){
        passdxy[ibool] = false;
      }
      int k = 0;
      for (int ibool = 0; ibool < nhist_pt ; ibool++){
        for (int jbool = 0; jbool < nhist_dxy ; jbool++){
          passcomb[k] = false;
          k++;
        }
      }
      
     unsigned int nhltjets  = ev->Jets.size();
      unsigned int nhltmuons = ev->iterL3muons.size();
      unsigned int nL1muons  = ev->L1muons.size();
      L3muon_size     -> Fill( nhltmuons );
      L1muon_size     -> Fill( nL1muons );
      lumis           -> Fill(ev -> luminosityBlockNumber);
   
      // build my HLT path, eg. L1 22, Jet 10
      for (int il1=0; il1 < nL1muons; il1++)
      {
        L1MuonCand theL1 = ev -> L1muons.at(il1);
        if ( (theL1.pt >= cut_pt_l1[ilumi] && fabs(theL1.eta) < cut_eta_l1[ilumi] && theL1.quality > qualcut) ||
             (theL1.pt >= cut_pt_l1[0]     && fabs(theL1.eta) < cut_eta_l1[0]     && theL1.quality > qualcut))
          {
          for (int il3=0; il3 < nhltmuons; il3++){
            HLTMuonCand theL3 = ev -> iterL3muons.at(il3);
            if (matchToL1(theL3,theL1)){      
              for (int i_ptcut = 0; i_ptcut < nhist_pt; i_ptcut++){  
                if (theL3.pt >  cut_pt_hlt[i_ptcut] && passpt[i_ptcut]==false){ 
                  pt_hlt[i_ptcut] -> Fill(theL3.pt);
                  passpt[i_ptcut] = true;
                }
              }
                                    
               for (int i_dxycut = 0; i_dxycut < nhist_dxy; i_dxycut++){  
                if (fabs(theL3.dxy/theL3.dxy_error) >  cut_dxy_s[i_dxycut] && passdxy[i_dxycut]==false){ 
                  dxy_err_hlt[i_dxycut] -> Fill(fabs(theL3.dxy/theL3.dxy_error));
                  passdxy[i_dxycut] = true;
                }
              }
              int k = 0;
              for (int i_ptcut = 0; i_ptcut < nhist_pt; i_ptcut++){  
                for (int i_dxycut = 0; i_dxycut < nhist_dxy; i_dxycut++){  
                  if (theL3.pt >  cut_pt_hlt[i_ptcut] && 
                      fabs(theL3.dxy/theL3.dxy_error) >  cut_dxy_s[i_dxycut] &&
                      passcomb[k]==false) {
                        pt_dxy_hlt[k] -> Fill(theL3.pt);
                        passcomb[k] = true;
                  }
                  k++;
                }
              }
            }
          }
        }
      } //end loop on L1 muons       
    }
  
  }


  outfile            -> cd();
  lumis              -> Write();  
  L3muon_size        -> Write();
  L1muon_size        -> Write();

  for (int i=0; i < nhist_pt; i++) {
    pt_hlt[i]            -> Write();
  }
  for (int i=0; i < nhist_dxy; i++) {
    dxy_err_hlt[i]       -> Write();
  }
  {
    int k=0;
    for (int i=0; i < nhist_pt; i++) {
      for (int j=0; j < nhist_dxy; j++) {
        pt_dxy_hlt[k]       -> Write();
        k++;
      }
    }
  }
  return;
}



bool matchToL1(HLTMuonCand mu, L1MuonCand L1){

  bool match = false;
  float minDR = 0.3;
  float theDR = 100;
  theDR = deltaR(L1.etaAtVtx, L1.phiAtVtx, mu.eta, mu.phi);
  if (theDR < minDR){
    minDR = theDR;
    match = true;
  }
  return match;
}


bool matchToL1(HLTJetCand jet, L1MuonCand L1){

  bool match  = false;
  float minDR = 0.3;
  float theDR = 100;
  theDR = deltaR(L1.eta, L1.phi, jet.eta, jet.phi);
  if (theDR < minDR){
    minDR = theDR;
    match = true;
  }
  return match;
}

