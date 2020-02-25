#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLT/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"


void dorate          (TChain*, TFile*, int, TFile* );
bool matchToL1       (HLTMuonCand, L1MuonCand );
bool matchGenToL1Muon(HLTMuonCand , std::vector<GenParticleCand> , std::vector<GenParticleCand> );
   
float cut_pt_hlt[60];
const int nhist_pt  = 60;

float cut_dxy_s[60];
const int nhist_dxy  = 60;

float cut_pt_l1[]   = {22,    18.,   14.,  12.,  10.,   8.,   7.  , 6.};//10.5, 14, 16.5, 18.5, 19.5, 20, 20.5, 21., 22}; 
float cut_eta_l1[]  = {2.5,   1.5,   1.5,  1.5,  1.5,  1.5,   1.5 , 1.5};//10.5, 14, 16.5, 18.5, 19.5, 20, 20.5, 21., 22}; 
int   min_pu[]      = {46 ,   42 ,   37 ,  31 ,  25 ,   25,   25  , 25 };
int   max_pu[]      = {60 ,   52,    47,   41,   35,    35,   35  , 35 };


int qualcut        = 11;
int n_lumi_columns = 1;
int mylumicol      = 7;

void readNtuples_ROCPt_QCD(){

  for (int i=0; i < nhist_dxy; i++){
    cut_dxy_s[i] = 0.2*i - 0.2;
  }
  for (int i=0; i < nhist_pt; i++){
    cut_pt_hlt[i] = 0.5*i;
  }
  std::string thestr = "_ROC_rate_qscale_17GeV_5e33";

  TChain* tree1   = new TChain("muonNtuples/muonTree");
  TFile* tfile1   = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2018/forKstee/QCD_Pt-15to3000_TuneCP5_Flat_13TeV_pythia8/crab_QCD_15_3000_muNtuples/skim_hasL1Mu_qscale17_30.root", "r");
  TFile* outfile1 = TFile::Open(Form("QCD15to3000%s.root", thestr.c_str()),"RECREATE");
  tree1 -> Add(tfile1->GetName());
  dorate(tree1, outfile1, 0, tfile1);
  tree1    -> Delete();   tfile1   -> Close();   outfile1 -> Close();  

//   TFile* tfile1   = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2017/forKstee/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/crab_QCD_15_30_muNtuples_addPtHatInfo/skim_hasL1Mu_filterqScale_cut17GeV_fixNden.root", "r");
//   TFile* outfile1 = TFile::Open(Form("ptplot_QCD/QCD15to30%s.root", thestr.c_str()),"RECREATE");
//   tree1 -> Add(tfile1->GetName());
//   dorate(tree1, outfile1, 0, tfile1);
//   tree1    -> Delete();   tfile1   -> Close();   outfile1 -> Close();  
//   
//   TChain* tree0   = new TChain("muonNtuples/muonTree");
//   TFile* tfile0   = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2017/forKstee/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/crab_QCD_30_50_muNtuples_addPtHatInfo/skim_hasL1Mu_filterqScale_cut20GeV_fixNden.root", "r");
//   TFile* outfile0 = TFile::Open(Form("ptplot_QCD/QCD30to50%s.root", thestr.c_str()),"RECREATE");
//   tree0 -> Add(tfile0->GetName());
//   dorate(tree0, outfile0, 0, tfile0);
//   tree0    -> Delete();  tfile0   -> Close();  outfile0 -> Close();  
//   
//   TChain* tree2   = new TChain("muonNtuples/muonTree");
//   TFile* tfile2   = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2017/forKstee/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/crab_QCD_50_80_muNtuples_addPtHatInfo/skim_hasL1Mu_filterqScale_cut20GeV_fixNden.root", "r");
//   TFile* outfile2 = TFile::Open(Form("ptplot_QCD/QCD50to80%s.root", thestr.c_str()),"RECREATE");
//   tree2 -> Add(tfile2->GetName());
//   dorate(tree2, outfile2, 0, tfile2);
//   tree2 -> Delete();   tfile2 -> Close();  outfile2 -> Close();  
// 
//   TChain* tree3   = new TChain("muonNtuples/muonTree");
//   TFile* tfile3   = new TFile("/eos/cms/store/group/phys_bphys/fiorendi/13TeV/data2017/forKstee/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/crab_QCD_80_120_muNtuples_addPtHatInfo/skim_hasL1Mu_filterqScale_cut20GeV_fixNden.root", "r");
//   TFile* outfile3 = TFile::Open(Form("ptplot_QCD/QCD80to120%s.root", thestr.c_str()),"RECREATE");
//   tree3 -> Add(tfile3->GetName());
//   dorate(tree3, outfile3, 0, tfile3);
//   tree3 -> Delete();  tfile3 -> Close();  outfile3 -> Close();  

  return;

}

void dorate( TChain* tree,  TFile* outfile, int reducedEntries, TFile* infile){

  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TH1F* muon_size = new TH1F("muon_size" ,"muon_size" ,   10,  0, 10 );

  TH1F *pt_hlt            [nhist_pt];
  TH1F *pt_hlt_fromb      [nhist_pt];
  TH1F *dxy_err_hlt       [nhist_dxy];
  TH1F *dxy_err_hlt_fromb [nhist_dxy];
  TH1F *pt_dxy_hlt        [nhist_dxy*nhist_pt];
  TH1F *pt_dxy_hlt_fromb  [nhist_dxy*nhist_pt];

  for (int i=0; i < nhist_pt; i++) {
    pt_hlt[i]            = new TH1F(Form("pt_hlt_%d"  , i),     "", 6000,  0,  600);
    pt_hlt_fromb[i]      = new TH1F(Form("pt_hlt_fromb_%d", i), "", 6000,  0,  600);
  }  
  for (int i=0; i < nhist_dxy; i++) {
    dxy_err_hlt[i]       = new TH1F(Form("dxy_err_hlt_%d"  , i),     "", 1000,  0,  500);
    dxy_err_hlt_fromb[i] = new TH1F(Form("dxy_err_hlt_fromb_%d", i), "", 1000,  0,  500);
  }

  int k = 0;
  for (int i=0; i < nhist_pt; i++) {
      for (int j=0; j < nhist_dxy; j++) { 
          pt_dxy_hlt[k]       = new TH1F(Form("pt_dxy_hlt_%d"       , k),     "", 1000,  0,  500);
          pt_dxy_hlt_fromb[k] = new TH1F(Form("pt_dxy_hlt_fromb_%d" , k),     "", 1000,  0,  500);
          k++;
      }
  }
  
  bool passpt[nhist_pt];
  bool passdxy[nhist_dxy];
  bool passcomb[nhist_dxy*nhist_pt];

  bool passpt_nob[nhist_pt];
  bool passdxy_nob[nhist_dxy];
  bool passcomb_nob[nhist_dxy*nhist_pt];

  MuonEvent* ev = new MuonEvent();
  tree -> SetBranchAddress("event", &ev);

  int count = 0;
  int nentries = tree->GetEntries();
  if (reducedEntries ==1) nentries = nentries/100;
  std::cout << "Number of entries = " << nentries << std::endl;


  for (int ilumi = 0; ilumi < n_lumi_columns; ilumi++){
  
    ilumi = mylumicol;
    
    for (Int_t eventNo=0; eventNo < nentries; eventNo++)
    {
  
      // reset all booleans
      for (int ibool = 0; ibool < nhist_pt ; ibool++){
        passpt[ibool]= false;
        passpt_nob[ibool]= false;
      }
      for (int ibool = 0; ibool < nhist_dxy ; ibool++){
        passdxy[ibool] = false;
        passdxy_nob[ibool] = false;
      }
      int k = 0;
      for (int ibool = 0; ibool < nhist_pt ; ibool++){
        for (int jbool = 0; jbool < nhist_dxy ; jbool++){
          passcomb_nob[k] = false;
          passcomb[k] = false;
          k++;
        }
      }

      Int_t IgetEvent   = tree   -> GetEvent(eventNo);
      // define pileup range being considered
      if ( ev -> trueNI < min_pu[ilumi] || ev -> trueNI > max_pu[ilumi] ) continue;
      count++;
  
      unsigned int nL1muons  = ev->L1muons.size();
      unsigned int nhltmuons = ev->iterL3muons.size();
      muon_size -> Fill( nhltmuons );
      if (nL1muons == 0 || nhltmuons == 0) continue;
      
      // loop on L1 muons and find those passing the L1 selection
      for (int il1=0; il1 < nL1muons; il1++){
        L1MuonCand theL1 = ev -> L1muons.at(il1);
        if ( (theL1.pt >= cut_pt_l1[ilumi] && fabs(theL1.eta) < cut_eta_l1[ilumi] && theL1.quality > qualcut) ||
             (theL1.pt >= cut_pt_l1[0]     && fabs(theL1.eta) < cut_eta_l1[0]     && theL1.quality > qualcut))
        {
  
          // loop on L3 muons and find those matched to the L1 and passing the possible L3 selections
          for (int il3=0; il3 < nhltmuons; il3++){
            HLTMuonCand theL3 = ev -> iterL3muons.at(il3);

            if (matchToL1(theL3,theL1)){   
              // apply pt cut only
              for (int i_ptcut = 0; i_ptcut < nhist_pt; i_ptcut++){  
                if (theL3.pt >  cut_pt_hlt[i_ptcut] && passpt_nob[i_ptcut]==false ){ 
                  pt_hlt[i_ptcut]      -> Fill(theL3.pt);
                  passpt_nob[i_ptcut]=true;
                }
              }
              // apply dxy cut only
              for (int i_dxycut = 0; i_dxycut < nhist_dxy; i_dxycut++){  
                if (fabs(theL3.dxy/theL3.dxy_error) >  cut_dxy_s[i_dxycut] && passdxy_nob[i_dxycut]==false){ 
                   dxy_err_hlt[i_dxycut] -> Fill(fabs(theL3.dxy/theL3.dxy_error));
                   passdxy_nob[i_dxycut]=true;
                }
              }
                    
              // apply pt and dxy cuts in combination
              int k = 0;
              for (int i_ptcut = 0; i_ptcut < nhist_pt; i_ptcut++){  
                for (int i_dxycut = 0; i_dxycut < nhist_dxy; i_dxycut++){  
                  if (theL3.pt >  cut_pt_hlt[i_ptcut] && 
                      fabs(theL3.dxy/theL3.dxy_error) >  cut_dxy_s[i_dxycut] &&
                      passcomb_nob[k]==false) {
                    pt_dxy_hlt[k] -> Fill(theL3.pt);
                    passcomb_nob[k]=true;
                  }
                  k++;
                }
              }
            }
                  
            // now require that the muon is coming from a generic B      
            if ( (matchGenToL1Muon(theL3, ev -> genMuons, ev -> genBs))){
              if (matchToL1(theL3,theL1)){   
                for (int i_ptcut = 0; i_ptcut < nhist_pt; i_ptcut++){  
                  if (theL3.pt >  cut_pt_hlt[i_ptcut] && passpt[i_ptcut]==false ){ 
                    pt_hlt_fromb[i_ptcut]      -> Fill(theL3.pt);
                    passpt[i_ptcut]=true;
                  }
                }
                for (int i_dxycut = 0; i_dxycut < nhist_dxy; i_dxycut++){  
                  if (fabs(theL3.dxy/theL3.dxy_error) >  cut_dxy_s[i_dxycut] && passdxy[i_dxycut]==false){ 
                    dxy_err_hlt_fromb[i_dxycut] -> Fill(fabs(theL3.dxy/theL3.dxy_error));
                    passdxy[i_dxycut]=true;
                  }
                }
                int k = 0;
                for (int i_ptcut = 0; i_ptcut < nhist_pt; i_ptcut++){  
                  for (int i_dxycut = 0; i_dxycut < nhist_dxy; i_dxycut++){  
                    if (theL3.pt >  cut_pt_hlt[i_ptcut] && 
                        fabs(theL3.dxy/theL3.dxy_error) >  cut_dxy_s[i_dxycut] &&
                        passcomb[k]==false){
                      pt_dxy_hlt_fromb[k] -> Fill(theL3.pt);
                      passcomb[k]=true;
                    }
                    k++;
                  }
                }
              }
            }
            
          } // end loop on L3 muons       
        } 
      } // end loop on L1 muons
    }  
  }

  std::cout << "events in PU range: " << count << std::endl;
  outfile            -> cd();
  
  muon_size                -> Write();
  for (int i=0; i < nhist_pt; i++) {
    pt_hlt[i]            -> Write();
    pt_hlt_fromb[i]      -> Write();
  }
  for (int i=0; i < nhist_dxy; i++) {
    dxy_err_hlt[i]       -> Write();
    dxy_err_hlt_fromb[i] -> Write();
  }
  {
    int k=0;
    for (int i=0; i < nhist_pt; i++) {
      for (int j=0; j < nhist_dxy; j++) {
        pt_dxy_hlt[k]       -> Write();
        pt_dxy_hlt_fromb[k] -> Write();
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


bool matchGenToL1Muon(HLTMuonCand mu, std::vector<GenParticleCand> gen, std::vector<GenParticleCand> genB){

  bool match = false;
  int ngen = gen.size();

  bool isfromB     = false; 
  bool hasAnotherB = false;
  float genB_muonMother_eta, genB_muonMother_phi;
  float minDR = 0.1;
  float theDR = 100;
  float mu_b_dr = 100;
  for ( std::vector<GenParticleCand>::const_iterator it = gen.begin(); it != gen.end(); ++it ) {
    for (int imother = 0; imother < it ->pdgMotherId.size(); imother++){
      if (fabs(it -> pdgMotherId.at(imother)) > 500 && fabs(it -> pdgMotherId.at(imother)) < 600) {
        isfromB = true;
        genB_muonMother_eta = it -> pdgMotherEta.at(imother);
        genB_muonMother_phi = it -> pdgMotherPhi.at(imother);
        break;
      }
    }
    if (!isfromB) continue;
    for ( std::vector<GenParticleCand>::const_iterator b_it = genB.begin(); b_it != genB.end(); ++b_it ) {
      if (fabs(b_it->eta) > 2.5) continue;
      if (b_it -> pt      < 3  ) continue;
      mu_b_dr = deltaR(genB_muonMother_eta, genB_muonMother_phi, b_it -> eta, b_it -> phi);
      if (mu_b_dr > 0.05) hasAnotherB = true;
    }
    if (!hasAnotherB) continue;
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
      break;
    }
  }
  return match;
}
