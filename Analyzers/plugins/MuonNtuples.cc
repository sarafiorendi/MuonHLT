/** \class MuonNtuples
 */
      
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

#include "MuonHLT/Analyzers/src/MuonTree.h"


class MuonNtuples : public edm::EDAnalyzer {

 public:
  MuonNtuples(const edm::ParameterSet& cfg);
  virtual ~MuonNtuples() {};

  virtual void analyze (const edm::Event& event, const edm::EventSetup & eventSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginEvent();
  virtual void beginRun(const edm::Run & run,    const edm::EventSetup & eventSetup);
  virtual void endRun  (const edm::Run & run,    const edm::EventSetup & eventSetup);

 private:

  void fillHlt(const edm::Handle<edm::TriggerResults> &, 
	           const edm::Handle<trigger::TriggerEvent> &,
	           const edm::TriggerNames &,
	           const edm::Event &);

  void fillMuons(const edm::Handle<reco::MuonCollection> &,
		         const reco::Vertex &, 
		         const edm::Event   & );

  void fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> &,
		            const reco::Vertex &, 
		            const edm::Event   & );

  void MonteCarloStudies(const edm::Event&);
  

//   virtual void endEvent();

  edm::InputTag offlinePVTag_;
  edm::InputTag offlineMuonTag_;
  /// file service
  edm::Service<TFileService> outfile_;

  // Trigger process
  std::string triggerProcess_;

  // Input tags
  edm::InputTag l3candTag_;
  edm::InputTag chargedDepTag_;
  edm::InputTag neutralDepTag_;
  edm::InputTag photonsDepTag_;
  edm::InputTag neutralDepTag05_;
  edm::InputTag photonsDepTag05_;
  edm::InputTag neutralDepTag1_;
  edm::InputTag photonsDepTag1_;
  edm::InputTag rhoCorrectionTag_;
  edm::InputTag rhoCorrectionOfflineTag_;
  edm::InputTag beamSpotTag_;
  edm::InputTag offlineECalPFTag03_;
  edm::InputTag offlineHCalPFTag03_;
  edm::InputTag offlineECalPFTag04_;
  edm::InputTag offlineHCalPFTag04_;
  
  MuonEvent event_;
  std::map<std::string,TTree*> tree_;
  
  unsigned int nGoodVtx; 


};

/// default constructor
MuonNtuples::MuonNtuples(const edm::ParameterSet& cfg): 
  offlinePVTag_           (cfg.getParameter<edm::InputTag>("offlineVtx")), 
  offlineMuonTag_         (cfg.getParameter<edm::InputTag>("offlineMuons")), 
  triggerProcess_         (cfg.getParameter<std::string>("triggerProcess")), 
  l3candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L3Candidates")),
  chargedDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("ChargedDeposit")), 
  neutralDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("NeutralDeposit")), 
  photonsDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDeposit")), 
  neutralDepTag05_        (cfg.getUntrackedParameter<edm::InputTag>("NeutralDeposit05")), 
  photonsDepTag05_        (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDeposit05")), 
  neutralDepTag1_         (cfg.getUntrackedParameter<edm::InputTag>("NeutralDeposit1" )), 
  photonsDepTag1_         (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDeposit1" )), 
  rhoCorrectionTag_       (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionOnline")), 
  rhoCorrectionOfflineTag_(cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionOffline")), 
  beamSpotTag_            (cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag")),
  offlineECalPFTag03_     (cfg.getUntrackedParameter<edm::InputTag>("offlineECalPFIso03")),
  offlineHCalPFTag03_     (cfg.getUntrackedParameter<edm::InputTag>("offlineHCalPFIso03")),
  offlineECalPFTag04_     (cfg.getUntrackedParameter<edm::InputTag>("offlineECalPFIso04")),
  offlineHCalPFTag04_     (cfg.getUntrackedParameter<edm::InputTag>("offlineHCalPFIso04"))
{}

void MuonNtuples::beginJob() {

  TH1::SetDefaultSumw2() ;
  tree_["muonTree"] = outfile_-> make<TTree>("muonTree","muonTree");
  tree_["muonTree"] -> Branch("event" ,&event_, 64000,2);

}    

void MuonNtuples::endJob() {}

void MuonNtuples::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}

void MuonNtuples::endRun  (const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonNtuples::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) {

  beginEvent();

  // Fill general info
  event_.runNumber             = event.id().run();
  event_.luminosityBlockNumber = event.id().luminosityBlock();
  event_.eventNumber           = event.id().event();


  // Fill vertex info
  edm::Handle<reco::VertexCollection> vertices; 
  event.getByLabel(offlinePVTag_, vertices);
  for(reco::VertexCollection::const_iterator it = vertices->begin(); it != vertices->end(); ++it) {
    if( !it->isValid())  continue;
    nGoodVtx++;
  }

  event_.nVtx = nGoodVtx;
  const reco::Vertex           & pv      = vertices->at(0);


  // Fill offline rho info
  edm::Handle <double>  rhoCollectionOffline;
  event.getByLabel(rhoCorrectionOfflineTag_, rhoCollectionOffline);
  if (rhoCollectionOffline.isValid()) event_.rho     = *(rhoCollectionOffline.product());


  // Fill PU info
  if (!event.isRealData()) {
	edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
	if ( event.getByLabel("addPileupInfo",puInfo)){
	  std::vector<PileupSummaryInfo>::const_iterator PVI;
	  for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) 
	  {
		if(PVI->getBunchCrossing()==0){
		  event_.trueNI   = PVI->getTrueNumInteractions();
		  continue;
		}
	  }
	} 
	else  
	  edm::LogError("") << "PU collection not found !!!";
  }
  
  // Fill MC GEN info
  if (!event.isRealData()) 
    MonteCarloStudies(event);

  
  // Fill trigger information
  edm::Handle<edm::TriggerResults>   triggerResults;
  edm::Handle<trigger::TriggerEvent> triggerEvent;
      
  if (event.getByLabel(edm::InputTag("TriggerResults"      , "", triggerProcess_), triggerResults) &&
	  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent)) {
	  
    edm::TriggerNames triggerNames_ = event.triggerNames(*triggerResults);
	fillHlt(triggerResults, triggerEvent, triggerNames_, event);
  }
  else 
	edm::LogError("") << "Trigger collection not found !!!";


  // Handle the offline muon collection and fill offline muons
  edm::Handle<std::vector<reco::Muon> > muons;
  event.getByLabel(offlineMuonTag_, muons);
  fillMuons(muons, pv, event);


  // Handle the online muon collection and fill online muons
  edm::Handle<reco::RecoChargedCandidateCollection> l3cands;
  if (event.getByLabel(l3candTag_, l3cands))
    fillHltMuons(l3cands, pv, event);
  else
    edm::LogWarning("") << "Online muon collection not found !!!";
  
//   endEvent();
  tree_["muonTree"] -> Fill();
}

//------------------------------------------------------------------------
void MuonNtuples::MonteCarloStudies(const edm::Event& event)
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByLabel("genParticles", genParticles);
  int ZId   =    23;
  int muId  =    13;

  for ( size_t i=0; i< genParticles->size(); ++i) 
  { 
    const reco::GenParticle &p = (*genParticles)[i];
  	if(fabs(p.pdgId()) != ZId ) 	continue; 

    // return all daughters of Z 
    if (p.numberOfDaughters()!= 2 ) continue;  
    for ( size_t ides=0; ides < p.numberOfDaughters(); ++ides ) 
    {
  	  const reco::Candidate *des = p.daughter(ides);
	  int dauId = des->pdgId();
      if( fabs(dauId) != muId ) continue;
      
      GenParticleCand theGen;
      theGen.pdgId  = des -> pdgId();
      theGen.pt     = des -> pt() ;
      theGen.eta    = des -> eta();
      theGen.phi    = des -> phi();
      theGen.energy = des -> energy();
      theGen.status = des -> status();

      event_.genParticles.push_back(theGen);

	} // end for des
  }  // end for genParticles
}

// --------------------------------------------------------------------
void MuonNtuples::fillHlt(const edm::Handle<edm::TriggerResults>   & triggerResults, 
				          const edm::Handle<trigger::TriggerEvent> & triggerEvent  ,
				          const edm::TriggerNames                  & triggerNames  ,
				          const edm::Event                         & event         )
{    
   
  for (unsigned int itrig=0; itrig < triggerNames.size(); ++itrig) 
  {
    LogDebug ("triggers") << triggerNames.triggerName(itrig) ;
    if (triggerResults->accept(itrig)) 
	{
	  std::string pathName = triggerNames.triggerName(itrig);
	  event_.hlt.triggers.push_back(pathName);
	}
  }
     
     
  const trigger::size_type nFilters(triggerEvent->sizeFilters());
  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
  {
    std::string filterTag = triggerEvent->filterTag(iFilter).encode();

    trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
    const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
	
    for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
	{  
	  trigger::size_type objKey = objectKeys.at(iKey);
	  const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
	  
	  HLTObjCand hltObj;
	  
	  hltObj.filterTag = filterTag;

	  hltObj.pt  = triggerObj.pt();
	  hltObj.eta = triggerObj.eta();
	  hltObj.phi = triggerObj.phi();
	  
	  event_.hlt.objects.push_back(hltObj);
	  
	}       
  }
  
  // fill hlt rho information
  edm::Handle <double>  hltRhoCollection;
  if (event.getByLabel(rhoCorrectionTag_, hltRhoCollection) && hltRhoCollection.isValid())
    event_.hlt.rho = *(hltRhoCollection.product());

}

// ---------------------------------------------------------------------
void MuonNtuples::fillMuons(const edm::Handle<reco::MuonCollection>       & muons ,
				            const reco::Vertex                            & pv    ,
				            const edm::Event                              & event )
{

  int n_mu = 0;
  for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) 
  { 

    n_mu++;
    MuonCand theMu;

    theMu.pt      = mu1 -> pt();
    theMu.eta     = mu1 -> eta();
    theMu.phi     = mu1 -> phi();
    theMu.charge  = mu1 -> charge();

    mu1 -> isGlobalMuon  () ? theMu.isGlobal  = 1 : theMu.isGlobal   = 0;
    mu1 -> isTrackerMuon () ? theMu.isTracker = 1 : theMu.isTracker  = 0;

    muon::isTightMuon ( (*mu1), pv ) ? theMu.isTight  = 1 : theMu.isTight  = 0;
    muon::isLooseMuon ( (*mu1)     ) ? theMu.isLoose  = 1 : theMu.isLoose  = 0;
    muon::isMediumMuon( (*mu1)     ) ? theMu.isMedium = 1 : theMu.isMedium = 0;
    
    theMu.chargedDep_dR03 = mu1->pfIsolationR03().sumChargedHadronPt ;
    theMu.neutralDep_dR03 = mu1->pfIsolationR03().sumNeutralHadronEt ;
    theMu.photonDep_dR03  = mu1->pfIsolationR03().sumPhotonEt        ;
    theMu.puPt_dR03       = mu1->pfIsolationR03().sumPUPt            ;

    theMu.chargedDep_dR04 = mu1->pfIsolationR04().sumChargedHadronPt ;
    theMu.neutralDep_dR04 = mu1->pfIsolationR04().sumNeutralHadronEt ;
    theMu.photonDep_dR04  = mu1->pfIsolationR04().sumPhotonEt        ;
    theMu.puPt_dR04       = mu1->pfIsolationR04().sumPUPt            ;


    // now fill PF clusters
	reco::MuonRef nmuonRef = reco::MuonRef(muons, n_mu-1);
    
    edm::Handle< edm::ValueMap<float> > muECalIsoMap03;
    edm::Handle< edm::ValueMap<float> > muHCalIsoMap03;
    if (event.getByLabel(offlineECalPFTag03_, muECalIsoMap03) &&
        event.getByLabel(offlineHCalPFTag03_, muHCalIsoMap03))
    { 
      const edm::ValueMap<float> muECalIso03 = *(muECalIsoMap03);
      const edm::ValueMap<float> muHCalIso03 = *(muHCalIsoMap03);
      theMu.ecalPFCluster_dR03 = muECalIso03[nmuonRef];
      theMu.hcalPFCluster_dR03 = muHCalIso03[nmuonRef];
    }
    else {
	  edm::LogWarning("") << "Offline PF cluster in dR 03 collection not found !!!";
      theMu.ecalPFCluster_dR03 = -9999 ;
      theMu.hcalPFCluster_dR03 = -9999 ;
    }

    edm::Handle< edm::ValueMap<float> > muECalIsoMap04;
    edm::Handle< edm::ValueMap<float> > muHCalIsoMap04;
    if (event.getByLabel(offlineECalPFTag04_, muECalIsoMap04) &&
        event.getByLabel(offlineHCalPFTag04_, muHCalIsoMap04) 
    ){ 
      const edm::ValueMap<float> muECalIso04 = *(muECalIsoMap04);
      const edm::ValueMap<float> muHCalIso04 = *(muHCalIsoMap04);
      theMu.ecalPFCluster_dR04 = muECalIso04[nmuonRef];
      theMu.hcalPFCluster_dR04 = muHCalIso04[nmuonRef];
    }
    else {
	  edm::LogWarning("") << "Offline PF cluster in dR 04 collection not found !!!";
      theMu.ecalPFCluster_dR04 = -9999 ;
      theMu.hcalPFCluster_dR04 = -9999 ;
    }


    event_.muons.push_back(theMu);
  }
}

// ---------------------------------------------------------------------
void MuonNtuples::fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
				               const reco::Vertex                                         & pv      ,
				               const edm::Event                                           & event )
{

  edm::Handle<reco::IsoDepositMap> trkDepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap;

  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap05;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap05;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap1 ;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap1 ;


  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
    HLTMuonCand theL3Mu;

	reco::RecoChargedCandidateRef candref(l3cands, il3);
    theL3Mu.pt      = candref -> pt();
    theL3Mu.eta     = candref -> eta();
    theL3Mu.phi     = candref -> phi();
    theL3Mu.charge  = candref -> charge();

    reco::TrackRef trkmu = candref->track();
    theL3Mu.trkpt   = trkmu -> pt();


    if (event.getByLabel(chargedDepTag_, trkDepMap)     &&
        event.getByLabel(neutralDepTag_, neutralDepMap) &&
        event.getByLabel(photonsDepTag_, photonsDepMap) ){
        
      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi = (*neutralDepMap).find( candref );
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi = (*photonsDepMap).find( candref );
	  reco::IsoDeposit theTkIsolation = (*trkDepMap)[candref];

      theL3Mu.hcalDep = hcal_mapi->val;
      theL3Mu.ecalDep = ecal_mapi->val; 
      theL3Mu.trkDep  = theTkIsolation.depositWithin(0.3);
    }
    else {
	  edm::LogWarning("") << "Online PF cluster collection not found !!!";
      theL3Mu.hcalDep =  -9999 ;
      theL3Mu.ecalDep =  -9999 ; 
      theL3Mu.trkDep  =  -9999 ;
    }


    // fill deposits with veto cones
    if (event.getByLabel(neutralDepTag05_, neutralDepMap05) &&
        event.getByLabel(photonsDepTag05_, photonsDepMap05) &&
        event.getByLabel(neutralDepTag1_,  neutralDepMap1 ) &&
        event.getByLabel(photonsDepTag1_,  photonsDepMap1 ) 
        ){
        
      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi05 = (*neutralDepMap05).find( candref );
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi05 = (*photonsDepMap05).find( candref );
      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi1  = (*neutralDepMap1 ).find( candref );
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi1  = (*photonsDepMap1 ).find( candref );

      theL3Mu.hcalDep05 = hcal_mapi05->val;
      theL3Mu.ecalDep05 = ecal_mapi05->val; 
      theL3Mu.hcalDep1  = hcal_mapi1 ->val;
      theL3Mu.ecalDep1  = ecal_mapi1 ->val; 
    }
    else {
	  edm::LogWarning("") << "Online PF cluster collection not found !!!";
      theL3Mu.hcalDep05 =  -9999 ;
      theL3Mu.ecalDep05 =  -9999 ; 
      theL3Mu.hcalDep1  =  -9999 ;
      theL3Mu.ecalDep1  =  -9999 ;
    }

    event_.hltmuons.push_back(theL3Mu);
  }
}


//---------------------------------------------
void MuonNtuples::beginEvent()
{

  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();
  event_.hlt.rho = -1;

  event_.genParticles.clear();
  event_.muons.clear();
  event_.hltmuons.clear();
  
  for (unsigned int ix=0; ix<3; ++ix) {
    event_.primaryVertex[ix] = 0.;
    for (unsigned int iy=0; iy<3; ++iy) {
      event_.cov_primaryVertex[ix][iy] = 0.;
    }
  }
  event_.nVtx   = -1;
  event_.trueNI = -1;
  event_.rho    = -1;
  
  nGoodVtx = 0; 
}



// define this as a plug-in
DEFINE_FWK_MODULE(MuonNtuples);
