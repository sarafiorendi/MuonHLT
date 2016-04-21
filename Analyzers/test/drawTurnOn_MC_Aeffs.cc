/* 
g++ -o drawTurnOn_MC_Aeffs drawTurnOn_MC_Aeffs.cc `root-config --cflags --libs` -lFoam -lMinuit
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>   
#include <sstream>
#include "Riostream.h"
#include <string>   
#include <vector>
#include <math.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h> 
#include <cassert>

void drawCurve() ;

void    drawEffIso ( TFile*, float*, float*, std::stringstream &, float, float);
TGraph* doGraph    ( int   , std::stringstream &     );
void    setIso     ( float*, float*);

std::stringstream ss;
std::stringstream isoss, savess, legss;
std::string isotype;
std::string strsave;
bool savePdf  = false;
bool debug    = false;
bool isBarrel = false;
bool isEndcap = false;
 
const int num  =  50 ;
double minPt   =    0.;
double maxPt   = 1000.;
bool   limitPt = false;

//define wp to draw
int ebwp    = 11;
int eewp    =  8;
int hbwp    = 21;
int hewp    = 22;
int tbwp    =  9;
int tewp    =  9;

int ebnwp   = 11;
int eenwp   =  8;
int hbnwp   = 21;
int henwp   = 22;
int tbnwp   =  9;
int tenwp   =  9;

const int nfile = 12;
TFile *fEff[12];
float iso_x[num]; 
float iso_x_err[num]; 

int main(int argc, char**argv)
{
  if (argc != 3) { 
    std::cout << "type ecal_iso, hcal_iso, trk_iso, iso..., 1 for saving pdf and 0 for not" << std::endl; 
    return 0;
  }
  isotype = argv[1];
  strsave = argv[2];
  if (strsave.find("1")!=std::string::npos) savePdf = true;

  TApplication app("App",&argc, argv);
  drawCurve() ;
  std::cout << "Program finished" << std::endl;
  app.Run (); 
  return 0;
}

//----------------------------------------
void drawCurve()  
{
  gStyle->SetTextFont(42);
  gStyle->SetLegendFont(42);
// gROOT->SetBatch(true);

  for (int i = 0; i < nfile; i++){
    fEff[i]    = new TFile( Form("eff_aeff/efficiency_Aeff%d.root", i)                                  , "READ");
  }

  isoss.str(""); isoss << isotype;
  if (isotype.find("barrel")!=std::string::npos) {isBarrel = true; std::cout << "is barrel" << std::endl;}
  if (isotype.find("endcap")!=std::string::npos) {isEndcap = true; std::cout << "is endcap" << std::endl;}

  int thewp  = 0;
  int thenwp = 0;
  if (isotype.find("hcal")!=std::string::npos && isBarrel) { thewp = hbwp; thenwp = hbnwp; } 
  if (isotype.find("hcal")!=std::string::npos && isEndcap) { thewp = hewp; thenwp = henwp; }
  if (isotype.find("ecal")!=std::string::npos && isBarrel) { thewp = ebwp; thenwp = ebnwp; }
  if (isotype.find("ecal")!=std::string::npos && isEndcap) { thewp = eewp; thenwp = eenwp; }
  if (isotype.find("trk") !=std::string::npos && isBarrel) { thewp = tbwp; thenwp = tbnwp; }
  if (isotype.find("trk") !=std::string::npos && isEndcap) { thewp = tewp; thenwp = tenwp; }
  if (debug) std::cout << "iso wp is " << thewp << std::endl;

  setIso(iso_x, iso_x_err);

  TGraph* theRocs[12];
  TCanvas *c = new TCanvas("c","c",500,500);
  TMultiGraph *mg = new TMultiGraph();

  for (int i = 0; i < nfile; i ++){
    theRocs[i] = doGraph(i, isoss);
    mg->Add(theRocs[i] );
  }

  mg->Draw("APC");
  mg->GetYaxis()->SetTitleOffset(1.3);
  mg->GetXaxis()->SetTitle("Isolation cut");
  mg->GetYaxis()->SetTitle("Isolation efficiency on muons from DY");
  mg->GetXaxis()->SetRangeUser(0.0  ,0.15);
  mg->GetYaxis()->SetRangeUser(0.6 ,1.0);
  c->Modified();

   //build legend
  legss.str(""); 
  if (isotype.find("ecal")!=std::string::npos) {legss << "ECAL";}
  if (isotype.find("hcal")!=std::string::npos) {legss << "HCAL";}
  if (isotype.find("trk") !=std::string::npos) {legss << "TRK" ;}

  if      (isBarrel) legss << " Barrel";
  else if (isEndcap) legss << " Endcap";
  else               legss << " ";
  if      (limitPt)  legss << ", " << minPt << " < p_{T} < " << maxPt << " GeV" ;
  
  TLegend *leg = new TLegend(0.45, 0.2, .86, .62, legss.str().c_str());
  leg -> SetLineColor(0);
  leg -> SetFillColor(0);
  leg -> SetTextSize(0.027);
  float theq;
  for (int i = 0; i < nfile; i ++){
    theq = 0.95 - i*0.05;
    leg->AddEntry(theRocs[i], Form("%.2f quantile", theq), "pl");
  }
  leg->Draw("same");

  gPad -> SetGrid();
  c    -> Update();   
  savess.str("");
  savess << "effVsIso_" << isoss.str().c_str() << "_AEffs.pdf"; 
  if (savePdf) c->SaveAs(savess.str().c_str());
}


// -----------------------------------------
TGraph* doGraph(int i, std::stringstream & isoss ){

  std::cout << "--------------------- prob: " << 0.95 - i*0.05 << "---------------------" << std::endl; 
  float isoEff0[num], errisoEff0[num];
  drawEffIso (fEff[i], isoEff0, errisoEff0, isoss, minPt, maxPt) ;

  TGraph *g0 = new TGraph(num, iso_x, isoEff0);
  g0 -> SetTitle(" ");
  g0 -> SetMarkerStyle(23);
  g0 -> SetMarkerSize(0.1);
  g0 -> SetLineWidth(2);
  g0 -> SetLineColor  (kAzure + i);
  g0 -> SetMarkerColor(kAzure + i);
  if (i > 4) {
    g0 -> SetLineColor  (kGreen + i - 5);
    g0 -> SetMarkerColor(kGreen + i - 5);
  }
  if (i > 9) {
    g0 -> SetLineColor  (kOrange + i - 9);
    g0 -> SetMarkerColor(kOrange + i - 9);
  }
  g0 -> GetXaxis() -> SetRangeUser(0.5,1);

  return g0;

}


//--------------------------------
void drawEffIso(TFile* fEff, float* eff, float* errEff, std::stringstream & isoss, float minPt, float maxPt){

  try{
    fEff -> GetName();
    if (debug) std::cout << "found file" << std::endl;
  }
  catch (...) {
    std::cout << "efficiency file does not exist" << std::endl;
    assert(0);   
  }    

  ss.str("");   
  ss << "HMuonPt_" << isoss.str().c_str()  << "0";
  TEfficiency* tempEff;
  TH1F*        temph;
  try{
    tempEff = (TEfficiency*) fEff    -> Get(ss.str().c_str());
    temph   = (TH1F*)        tempEff -> GetTotalHistogram();
    if (debug) std::cout << "found hist" << std::endl;
  }
  catch (...) {
    std::cout << "TEfficiency denominator" << ss.str().c_str() << "does not exist" << std::endl;
    assert(0);   
  }    


  int bin_max, bin_min;
  try{
    bin_min = temph -> FindBin(minPt);
    bin_max = temph -> FindBin(maxPt);
  }
  catch (...) {
    std::cout << "bins not found" << std::endl;
    assert(0);   
  }    
  

  for(int i = 0; i < num; i++)
  {
	ss.str("");   
	ss << "HMuonPt_" << isoss.str().c_str()  << i;
	TEfficiency * hEff = (TEfficiency*) fEff->Get(ss.str().c_str());
	
	float isoEff_den  = ((TH1F*)hEff -> GetTotalHistogram() )  -> Integral(bin_min, bin_max);
	float isoEff_num  = ((TH1F*)hEff -> GetPassedHistogram())  -> Integral(bin_min, bin_max);
	eff[i]    = isoEff_num/isoEff_den;
	errEff[i] = sqrt(isoEff_num)/isoEff_den;
	
	if (eff[i] > 0.97 && eff[i] < 0.99) std::cout << i << ": " << eff[i] << std::endl;
//     hEff -> Draw("AP");
	if (debug) std::cout << fEff->GetName() << " - isoEff " << i << ": " << eff[i] << " +- " << errEff[i] << std::endl;
  }
}

//--------------------------------
void setIso(float* iso, float* iso_err){
  
  for(int i =0; i< num; i++)
  {
    iso[i]     = i * 0.01; 
    iso_err[i] = 0.005; 
  }
}


