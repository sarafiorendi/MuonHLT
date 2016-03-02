#include <iostream>
#include <sstream>


float lumi_scale   =  1.  ;
int   ls_len       =  23.3; 
int   n_ls         =  10  ; 

struct rateError{
    float rate;
    float error;
    float purity;
};

rateError evalOneRate(TFile* );

void computeRateOnData(){
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile *data_file        = TFile::Open("rates/rate_DYLL_pu30_pt28_newtight.root", "r");
    rateError rate_data     = evalOneRate(data_file);
    
    float sumRate = rate_data.rate;
    float errRate = rate_data.error; 
    
    std::cout << "total rate : " << sumRate << " +/- " << errRate << " Hz" <<  std::endl;

}

// ********************************************************************


rateError evalOneRate(TFile* file){

  std::stringstream sname;
  sname << file->GetName();
  
  float xsec;
  int nevents;
  
  std::cout <<                                      std::endl;
  std::cout << "input file: " << file->GetName() << std::endl;

  TH1F* h = (TH1F*)file->  Get("wp_passing");
  int npass = h -> GetEntries();


  float totalRate = lumi_scale * npass / (ls_len * n_ls) ;
  float errorRate = totalRate/sqrt(npass);
  std::cout << "n passing:   " << npass      << std::endl;
  std::cout << "n ls:        " << n_ls       << std::endl;
  std::cout << "rate:        " << totalRate  << " +/- " << errorRate << " Hz" << std::endl;
    
  rateError theRateAndError;
  theRateAndError.rate    = totalRate;
  theRateAndError.error   = errorRate;
//     theRateAndError.purity = purity;
//     cout << "rate=" << sumRate << " error=" << error << endl;
    
    return theRateAndError;
}


