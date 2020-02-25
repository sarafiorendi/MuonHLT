#ifndef  MuonTree_h
#define  MuonTree_h

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>



class GenParticleCand {
public:
  Int_t   pdgId; 
  Int_t   status; 
  Float_t energy; 
  Float_t pt; 
  Float_t eta; 
  Float_t phi; 
  bool isPromptDecayed; 
  std::vector<Int_t>    pdgMotherId; 
  std::vector<Float_t>  pdgMotherPt; 
  std::vector<Float_t>  pdgMotherEta; 
  std::vector<Float_t>  pdgMotherPhi; 
  std::vector<Int_t>    pdgRealMother; 

  GenParticleCand(){};
  virtual ~GenParticleCand(){};
  
  ClassDef(GenParticleCand,1)
};




class MuonCand {
public:

  Float_t pt;  
  Float_t eta; 
  Float_t phi; 
  Int_t   charge;    

  Int_t   isGlobal;
  Int_t   isTracker;

  Int_t   isLoose;
  Int_t   isMedium;
  Int_t   isTight;
  Int_t   isSoft;
  
  MuonCand(){};
  virtual ~MuonCand(){};

  ClassDef(MuonCand,1)
};


class HLTMuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t trkpt;         // pt of the track of the hlt muon [GeV]
  Int_t   charge;         // pt of the track of the hlt muon [GeV]
  
  Float_t dz;
  Float_t dxy;
  Float_t dxy_error;

  HLTMuonCand(){};
  virtual ~HLTMuonCand(){};

  ClassDef(HLTMuonCand,1)

};


class L2MuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Int_t   charge;      

  Int_t   minStations;      
  Int_t   minHits    ;      
  
  L2MuonCand(){};
  virtual ~L2MuonCand(){};

  ClassDef(L2MuonCand,1)

};


class L1MuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t etaAtVtx;          
  Float_t phiAtVtx;          
  Int_t   charge;      
  Int_t   quality;      
  
  L1MuonCand(){};
  virtual ~L1MuonCand(){};

  ClassDef(L1MuonCand,1)

};


class HLTJetCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Float_t trkpt;         // pt of the track of the hlt muon [GeV]
  Int_t   charge;         // pt of the track of the hlt muon [GeV]
  
  HLTJetCand(){};
  virtual ~HLTJetCand(){};

  ClassDef(HLTJetCand,1)

};


class HLTObjCand {
public:

  std::string filterTag; // name of filter passed by the object
  Float_t pt;            // pt of the object passing the filter [GeV]
  Float_t eta;           // eta of the object passing the filter
  Float_t phi;           // phi of the object passing the filter
  
  HLTObjCand(){};
  virtual ~HLTObjCand(){};

  ClassDef(HLTObjCand,1)

};





class HLTInfo {
public:
  std::vector<std::string>  triggers;  
  std::vector<HLTObjCand>   objects;   


  HLTInfo(){};
  virtual ~HLTInfo(){};
  bool match( const std::string & path ) {
	if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )  return true;
//     if (! iname.compare("HLT_Mu20_v1") == 0) continue;
	return false;
  }

  bool find( const std::string & path ) {
	for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
//       std::cout << *it << std::endl;
      if ( it-> compare(path) == 0) return true;
//       if ( it->find ( path ) != std::string::npos ) return true;
	}
	return false;
  }

  ClassDef(HLTInfo,1)

};


class MuonEvent {
public:

  Int_t   runNumber;             
  Int_t   luminosityBlockNumber; 
  Int_t   eventNumber;           

  Int_t   nVtx;                    
  Float_t primaryVertex[3];        
  Float_t cov_primaryVertex[3][3]; 

  Float_t trueNI;   
  Float_t rho; 
  Float_t qScale;
  Float_t max_pt_hats;
  
  Float_t bxId;
  Float_t instLumi; 

  std::vector <GenParticleCand> genMuons; 
  std::vector <GenParticleCand> genBs; 
//   std::vector <GenParticleCand> genParticles; 
  std::vector <MuonCand>        muons;         
  std::vector <HLTMuonCand>     tkmuons;      
  std::vector <HLTMuonCand>     hltmuons;      
  std::vector <HLTMuonCand>     iterL3muons;      

  std::vector <HLTMuonCand>     L3OImuons;      
  std::vector <HLTMuonCand>     L3IOmuons;      
  std::vector <HLTMuonCand>     L3L1muons;      


  std::vector <L2MuonCand>      L2muons;      
  std::vector <L1MuonCand>      L1muons;      
  std::vector <HLTJetCand>      Jets;
  HLTInfo                       hlt;           
  HLTInfo                       hltTag;            

  MuonEvent(){};
  virtual ~MuonEvent(){};

  ClassDef(MuonEvent,1)
};


#endif

