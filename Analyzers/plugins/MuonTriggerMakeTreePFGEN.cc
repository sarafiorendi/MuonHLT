/** \class MuonTriggerMakeTreePFGEN
 *  Class to measure muon trigger efficiencies (very rough)
 *
 *  $Date: 2012/11/16 19:13:44 $
 *  $Revision: 1.3 $
 *  \authors D. Trocino - Northeastern University <daniele.trocino@cern.ch>
 */
      
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "TLorentzVector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TTree.h"

class MuonTriggerMakeTreePFGEN : public edm::EDAnalyzer {

 public:
  MuonTriggerMakeTreePFGEN(const edm::ParameterSet& cfg);
  virtual ~MuonTriggerMakeTreePFGEN() {};

 private:
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void beginEvent();
  virtual void endEvent();

  edm::InputTag vertexes_;
  edm::Service<TFileService> outfile_;

  // Trigger process
  std::string tagTriggerProcess_;
  std::string triggerProcess_;
  // Trigger names
  std::string tagTriggerName_;
  std::string triggerName_;
  std::string probeFilterDen_;
  std::string probeFilterNum_;
  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  // HLTConfig
  HLTConfigProvider hltConfig_;
  HLTConfigProvider hltConfigTag_;

  TLorentzVector p_mup,p_mum, p_z;
  bool boolGEN = false;
  void MonteCarloStudies(const edm::Event&);
  
  TTree* outTree_;
  int T_Run,  T_Lumi, T_Event,T_Nprim, T_TrueNI; 
  double T_hlt_rho;
  std::vector<int>    *T_GenPdgId;
  std::vector<double> *T_GenPt, *T_GenEta, *T_GenPhi;
  std::vector<double> *T_hltTagMatch, *T_hltMatch;
  std::vector<double> *T_hlt_MuPt, *T_hlt_MuTrkPt, *T_hlt_MuEta, *T_hlt_MuPhi, *T_hlt_trkDep;
  std::vector<float>  *T_hlt_HcalClusterDep,   *T_hlt_EcalClusterDep,   *T_hlt_HcalClusterDepV1,   *T_hlt_HcalClusterDepV5,   *T_hlt_EcalClusterDepV1,   *T_hlt_EcalClusterDepV5;
  std::vector<float>  *T_hlt_HcalClusterDepEt, *T_hlt_EcalClusterDepEt, *T_hlt_HcalClusterDepV1Et, *T_hlt_HcalClusterDepV5Et, *T_hlt_EcalClusterDepV1Et, *T_hlt_EcalClusterDepV5Et;
  std::vector<float>  *T_GenHLTMatch;
  

  //Deposits tag & l3 muons
  edm::InputTag l3candTag_;
  edm::InputTag chargedDepTag_;
  edm::InputTag neutralDepTag_;
  edm::InputTag neutralDepTagV1_;
  edm::InputTag neutralDepTagV5_;
  edm::InputTag photonsDepTag_;
  edm::InputTag photonsDepTagV1_;
  edm::InputTag photonsDepTagV5_;
  edm::InputTag neutralDepTagEt_;
  edm::InputTag neutralDepTagV1Et_;
  edm::InputTag neutralDepTagV5Et_;
  edm::InputTag photonsDepTagEt_;
  edm::InputTag photonsDepTagV1Et_;
  edm::InputTag photonsDepTagV5Et_;
  edm::InputTag rhoCorrectionTag_;
  edm::InputTag VtxTag_;
  edm::InputTag AllVtxTag_;
  
  //Debug 
  bool debug_;

};

/// default constructor
MuonTriggerMakeTreePFGEN::MuonTriggerMakeTreePFGEN(const edm::ParameterSet& cfg): 
  vertexes_               (cfg.getParameter<edm::InputTag>("vertexes")), 
  tagTriggerProcess_      (cfg.getParameter<std::string>("tagTriggerProcess")), 
  triggerProcess_         (cfg.getParameter<std::string>("triggerProcess")), 
  tagTriggerName_         (cfg.getParameter<std::string>("tagTriggerName")), 
  triggerName_            (cfg.getParameter<std::string>("triggerName")), 
  probeFilterDen_         (cfg.getParameter<std::string>("probeFilterDen")),
  probeFilterNum_         (cfg.getParameter<std::string>("probeFilterNum")),
  l3candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L3CandidatesLabel"      )),
  chargedDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("ChargedDepositLabel")), 
  neutralDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabel")), 
  neutralDepTagV1_        (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelVeto1")), 
  neutralDepTagV5_        (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelVeto5")), 
  photonsDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabel")), 
  photonsDepTagV1_        (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelVeto1")), 
  photonsDepTagV5_        (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelVeto5")), 
  neutralDepTagEt_        (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelEt")), 
  neutralDepTagV1Et_      (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelVeto1Et")), 
  neutralDepTagV5Et_      (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelVeto5Et")), 
  photonsDepTagEt_        (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelEt")), 
  photonsDepTagV1Et_      (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelVeto1Et")), 
  photonsDepTagV5Et_      (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelVeto5Et")), 
  rhoCorrectionTag_       (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionLabel")), 
  debug_                  (cfg.getUntrackedParameter<bool>("boolDebug")) 
{}

void MuonTriggerMakeTreePFGEN::beginJob() {

  outTree_ = outfile_->make<TTree>("outTree_", "outTree_");

  outTree_->Branch("Run"                           , &T_Run                         );
  outTree_->Branch("Lumi"                          , &T_Lumi                        );
  outTree_->Branch("Event"                         , &T_Event                       );
  outTree_->Branch("Nprim"                         , &T_Nprim                       );
  outTree_->Branch("TrueNI"                        , &T_TrueNI                      );
  
  outTree_->Branch("T_GenPdgId"                    , &T_GenPdgId                    );
  outTree_->Branch("T_GenPt"                       , &T_GenPt                       );
  outTree_->Branch("T_GenEta"                      , &T_GenEta                      );
  outTree_->Branch("T_GenPhi"                      , &T_GenPhi                      );
  outTree_->Branch("T_GenHLTMatch"                 , &T_GenHLTMatch                 );
 
  outTree_->Branch("T_hlt_rho"                     , &T_hlt_rho                     );
  outTree_->Branch("T_hlt_MuPt"                    , &T_hlt_MuPt                    );
  outTree_->Branch("T_hlt_MuTrkPt"                 , &T_hlt_MuTrkPt                 );
  outTree_->Branch("T_hlt_MuEta"                   , &T_hlt_MuEta                   );
  outTree_->Branch("T_hlt_MuPhi"                   , &T_hlt_MuPhi                   );
  outTree_->Branch("T_hlt_trkDep"                  , &T_hlt_trkDep                  );
  outTree_->Branch("T_hlt_EcalClusterDep"          , &T_hlt_EcalClusterDep          );
  outTree_->Branch("T_hlt_EcalClusterDepVeto0p1"   , &T_hlt_EcalClusterDepV1        );
  outTree_->Branch("T_hlt_EcalClusterDepVeto0p05"  , &T_hlt_EcalClusterDepV5        );
  outTree_->Branch("T_hlt_HcalClusterDep"          , &T_hlt_HcalClusterDep          );
  outTree_->Branch("T_hlt_HcalClusterDepVeto0p1"   , &T_hlt_HcalClusterDepV1        );
  outTree_->Branch("T_hlt_HcalClusterDepVeto0p05"  , &T_hlt_HcalClusterDepV5        );
  outTree_->Branch("T_hlt_EcalClusterDepEt"        , &T_hlt_EcalClusterDepEt        );
  outTree_->Branch("T_hlt_EcalClusterDepVeto0p1Et" , &T_hlt_EcalClusterDepV1Et      );
  outTree_->Branch("T_hlt_EcalClusterDepVeto0p05Et", &T_hlt_EcalClusterDepV5Et      );
  outTree_->Branch("T_hlt_HcalClusterDepEt"        , &T_hlt_HcalClusterDepEt        );
  outTree_->Branch("T_hlt_HcalClusterDepVeto0p1Et" , &T_hlt_HcalClusterDepV1Et      );
  outTree_->Branch("T_hlt_HcalClusterDepVeto0p05Et", &T_hlt_HcalClusterDepV5Et      );
  outTree_->Branch("T_hltMatch"                    , &T_hltMatch                    );
  outTree_->Branch("T_hltTagMatch"                 , &T_hltTagMatch                 );
  
}    

void MuonTriggerMakeTreePFGEN::endJob() {
}


void MuonTriggerMakeTreePFGEN::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonTriggerMakeTreePFGEN::analyze(const edm::Event &event, const edm::EventSetup &eventSetup) {

  using reco::Muon;
  
  beginEvent();

//   edm::Handle<reco::VertexCollection> pvHandle; 
//   event.getByLabel(vertexes_, pvHandle);
//   const reco::VertexCollection & vertices = *pvHandle.product();

  // Handle to the rho collection
  edm::Handle <double>  rhoCollection;
  event.getByLabel(rhoCorrectionTag_, rhoCollection);

  T_Run    = event.id().run();
  T_Lumi   = event.id().luminosityBlock();
  T_Event  = event.id().event();
//   T_Nprim  = vertices.size();
  
  if (rhoCollection.isValid()) T_hlt_rho = *(rhoCollection.product());

  edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
  try {
    event.getByLabel("addPileupInfo",puInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) 
    {
      if(PVI->getBunchCrossing()==0){
        T_TrueNI    = PVI->getTrueNumInteractions();
        continue;
      }
    }
  } catch (...) {}


  // GEN particles	
  MonteCarloStudies(event);
  if (!boolGEN) return;


  // Handle to the hlt muon collection
  edm::Handle<reco::RecoChargedCandidateCollection> l3cands;
  if (!event.getByLabel(l3candTag_, l3cands)){ 
	edm::LogError("") << " Muon collection does not exist !!!";
  }
  if (!l3cands.isValid()){
//      std::cout << "now filling" << std::endl;
	 T_hlt_MuTrkPt            ->push_back(-100); 
	 T_hlt_MuPt               ->push_back(-100); 
	 T_hlt_MuEta              ->push_back(-100); 
	 T_hlt_MuPhi              ->push_back(-100); 
	 T_hlt_EcalClusterDep     ->push_back(-100);
	 T_hlt_EcalClusterDepV1   ->push_back(-100);
	 T_hlt_EcalClusterDepV5   ->push_back(-100);
	 T_hlt_HcalClusterDep     ->push_back(-100);
	 T_hlt_HcalClusterDepV1   ->push_back(-100);
	 T_hlt_HcalClusterDepV5   ->push_back(-100);
	 T_hlt_EcalClusterDepEt   ->push_back(-100);
	 T_hlt_EcalClusterDepV1Et ->push_back(-100);
	 T_hlt_EcalClusterDepV5Et ->push_back(-100);
	 T_hlt_HcalClusterDepEt   ->push_back(-100);
	 T_hlt_HcalClusterDepV1Et ->push_back(-100);
	 T_hlt_HcalClusterDepV5Et ->push_back(-100);
	 T_hlt_trkDep             ->push_back(-100);
	 T_hltMatch               ->push_back(-100);
	 T_GenHLTMatch            ->push_back(-100);
     outTree_ -> Fill();
     endEvent();
     return;
  }
  
  
  edm::Handle<reco::IsoDepositMap> TrkDepMap;
  event.getByLabel(chargedDepTag_, TrkDepMap);

  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap;
  event.getByLabel(neutralDepTag_, neutralDepMap);
  event.getByLabel(photonsDepTag_, photonsDepMap);

  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV1;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV5;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV1;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV5;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapEt;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapEt;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV1Et;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV5Et;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV1Et;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV5Et;
  event.getByLabel(neutralDepTagV1_  , neutralDepMapV1  );
  event.getByLabel(neutralDepTagV5_  , neutralDepMapV5  );
  event.getByLabel(photonsDepTagV1_  , photonsDepMapV1  );
  event.getByLabel(photonsDepTagV5_  , photonsDepMapV5  );
  event.getByLabel(neutralDepTagEt_  , neutralDepMapEt  );
  event.getByLabel(photonsDepTagEt_  , photonsDepMapEt  );
  event.getByLabel(neutralDepTagV1Et_, neutralDepMapV1Et);
  event.getByLabel(neutralDepTagV5Et_, neutralDepMapV5Et);
  event.getByLabel(photonsDepTagV1Et_, photonsDepMapV1Et);
  event.getByLabel(photonsDepTagV5Et_, photonsDepMapV5Et);


  // Get trigger results
  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(edm::InputTag("TriggerResults", "", triggerProcess_), triggerResults);
  if(!triggerResults.isValid()) {
    std::cout << "Trigger results not valid" << std::endl;
    return;
  } 
  // Get trigger results
  edm::Handle<edm::TriggerResults> tagTriggerResults;
  event.getByLabel(edm::InputTag("TriggerResults", "", tagTriggerProcess_), tagTriggerResults);
  if(!tagTriggerResults.isValid()) {
    std::cout << "Trigger tag results not valid" << std::endl;
    return;
  } 

  //Print trigger in the event
  if (debug_){
	edm::TriggerNames triggerNames_ = event.triggerNames(*triggerResults);
    int ntrigs=triggerResults->size();
	std::vector<std::string> triggernames = triggerNames_.triggerNames();
	for(int itrig = 0; itrig != ntrigs; ++itrig)
	{
	  std::cout << triggerNames_.triggerName(itrig) << std::endl;
	}

	edm::TriggerNames tagTriggerNames_ = event.triggerNames(*tagTriggerResults);
    int ntrigsTag=tagTriggerResults->size();
	std::vector<std::string> triggernamestag = tagTriggerNames_.triggerNames();
	for(int itrig = 0; itrig != ntrigsTag; ++itrig)
	{
	  std::cout << "tag:" << tagTriggerNames_.triggerName(itrig) << std::endl;
	}
  }

//   if( !tagTriggerResults->accept(tagTriggerIndex_) ) return; // there are no tags


  // Get trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent);
  if(!triggerEvent.isValid()) { 
    std::cout << "TriggerEvent not valid" << std::endl;
    return;
  }

  assert(triggerResults->size()==hltConfig_.size());
  const trigger::TriggerObjectCollection & toc = triggerEvent->getObjects();

  // Loop muon collection and fill histograms
  if (l3cands->size() !=0){
    for(unsigned int il3=0; il3<l3cands->size(); ++il3) {
  
      reco::RecoChargedCandidateRef candref(l3cands, il3);
      // Check if it is matched to the last-1 hlt filter
      bool isGood = false; 
      const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_));
      const unsigned int m(hltConfig_.size(triggerIndex_));
      assert( moduleLabels.size()==m );
      const unsigned int lastModuleIndex(triggerResults->index(triggerIndex_));
      assert(lastModuleIndex<m);
    
      //If no filter is specified for the denominator, take all offline muons
      if( probeFilterDen_.length()==0 ) isGood = true; 
      else {   
        unsigned int filterIndex = triggerEvent->sizeFilters();
        for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
          if( probeFilterDen_.compare(moduleLabels[j])!=0 ) continue; 
          filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterDen_, "", triggerProcess_));
          break; 
        }
        if( filterIndex<triggerEvent->sizeFilters() ) {
          const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
          const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
          assert( vids.size()==keys.size() );
          const unsigned int nProbeTrig(vids.size());
          for(unsigned int i=0; i<nProbeTrig; ++i) {
            const trigger::TriggerObject & to = toc[keys[i]];
            if (fabs(to.pt() - candref->pt()) > 0.1 ) continue; 
            isGood = true;
            break;
          }
        } // end if( filterIndex<triggerEvent->sizeFilters() )
      } // end else

      if(!isGood) {
        T_hlt_MuTrkPt            ->push_back(-100); 
        T_hlt_MuPt               ->push_back(-100); 
        T_hlt_MuEta              ->push_back(-100); 
        T_hlt_MuPhi              ->push_back(-100); 
	    T_hlt_EcalClusterDep     ->push_back(-100);
	    T_hlt_EcalClusterDepV1   ->push_back(-100);
	    T_hlt_EcalClusterDepV5   ->push_back(-100);
	    T_hlt_HcalClusterDep     ->push_back(-100);
	    T_hlt_HcalClusterDepV1   ->push_back(-100);
	    T_hlt_HcalClusterDepV5   ->push_back(-100);
	    T_hlt_EcalClusterDepEt   ->push_back(-100);
	    T_hlt_EcalClusterDepV1Et ->push_back(-100);
	    T_hlt_EcalClusterDepV5Et ->push_back(-100);
	    T_hlt_HcalClusterDepEt   ->push_back(-100);
	    T_hlt_HcalClusterDepV1Et ->push_back(-100);
	    T_hlt_HcalClusterDepV5Et ->push_back(-100);
        T_hlt_trkDep             ->push_back(-100);
        T_hltMatch               ->push_back(-100);
        T_GenHLTMatch            ->push_back(-100);
        continue; 
      }
        
      float matchMuP = deltaR(candref->eta(),candref->phi(),p_mup.Eta(),p_mup.Phi());
      float matchMuM = deltaR(candref->eta(),candref->phi(),p_mum.Eta(),p_mum.Phi());
      if (matchMuP < matchMuM) T_GenHLTMatch->push_back(matchMuP);
      else                     T_GenHLTMatch->push_back(matchMuM);
  
      unsigned int filterIndex = triggerEvent->sizeFilters(); 
      for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
        if( probeFilterNum_.compare(moduleLabels[j])!=0 ) continue; 
        filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterNum_, "", triggerProcess_));
        break; 
      }
      if( filterIndex < triggerEvent->sizeFilters() ) {
        const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
        const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
        assert( vids.size()==keys.size() );
        const unsigned int nProbeTrig(vids.size());
  
        T_hlt_MuPt ->push_back(candref->pt()); 
        T_hlt_MuEta->push_back(candref->eta()); 
        T_hlt_MuPhi->push_back(candref->phi()); 
        reco::TrackRef trkmu = candref->track();
        T_hlt_MuTrkPt ->push_back(trkmu->pt());
        
        reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi = (*neutralDepMap).find( candref );
        T_hlt_HcalClusterDep -> push_back(hcal_mapi->val);
        reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi = (*photonsDepMap).find( candref );
        T_hlt_EcalClusterDep -> push_back(ecal_mapi->val); 

	    reco::IsoDeposit theTkIsolation = (*TrkDepMap)[candref];
        T_hlt_trkDep->push_back(theTkIsolation.depositWithin(0.3));

		reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV1 = (*neutralDepMapV1).find( candref );
		T_hlt_HcalClusterDepV1 -> push_back(hcal_mapiV1->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV5 = (*neutralDepMapV5).find( candref );
		T_hlt_HcalClusterDepV5 -> push_back(hcal_mapiV5->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV1 = (*photonsDepMapV1).find( candref );
		T_hlt_EcalClusterDepV1 -> push_back(ecal_mapiV1->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV5 = (*photonsDepMapV5).find( candref );
		T_hlt_EcalClusterDepV5 -> push_back(ecal_mapiV5->val);

		reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiet = (*neutralDepMapEt).find( candref );
		T_hlt_HcalClusterDepEt -> push_back(hcal_mapiet->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiet = (*photonsDepMapEt).find( candref );
		T_hlt_EcalClusterDepEt -> push_back(ecal_mapiet->val); 
		reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV1et = (*neutralDepMapV1Et).find( candref );
		T_hlt_HcalClusterDepV1Et -> push_back(hcal_mapiV1et->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV5et = (*neutralDepMapV5Et).find( candref );
		T_hlt_HcalClusterDepV5Et -> push_back(hcal_mapiV5et->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV1et = (*photonsDepMapV1Et).find( candref );
		T_hlt_EcalClusterDepV1Et -> push_back(ecal_mapiV1et->val);
		reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV5et = (*photonsDepMapV5Et).find( candref );
		T_hlt_EcalClusterDepV5Et -> push_back(ecal_mapiV5et->val);

        bool found = false;
        for(unsigned int i=0; i<nProbeTrig; ++i) {
          found = false;
          const trigger::TriggerObject & to = toc[keys[i]];
          if (fabs(to.pt() - candref->pt()) < 0.1) { 
             T_hltMatch->push_back(1);
             found = true;
             break;
          }
        }
        if (found==false) T_hltMatch->push_back(0);
      }
    }
  }//if l3 size != 0
  else {
    T_hlt_MuTrkPt            ->push_back(-100); 
    T_hlt_MuPt               ->push_back(-100); 
    T_hlt_MuEta              ->push_back(-100); 
    T_hlt_MuPhi              ->push_back(-100); 
	T_hlt_EcalClusterDep     ->push_back(-100);
	T_hlt_EcalClusterDepV1   ->push_back(-100);
	T_hlt_EcalClusterDepV5   ->push_back(-100);
	T_hlt_HcalClusterDep     ->push_back(-100);
	T_hlt_HcalClusterDepV1   ->push_back(-100);
	T_hlt_HcalClusterDepV5   ->push_back(-100);
	T_hlt_EcalClusterDepEt   ->push_back(-100);
	T_hlt_EcalClusterDepV1Et ->push_back(-100);
	T_hlt_EcalClusterDepV5Et ->push_back(-100);
	T_hlt_HcalClusterDepEt   ->push_back(-100);
	T_hlt_HcalClusterDepV1Et ->push_back(-100);
	T_hlt_HcalClusterDepV5Et ->push_back(-100);
    T_hlt_trkDep             ->push_back(-100);
    T_hltMatch               ->push_back(-100);
    T_GenHLTMatch            ->push_back(-100);
  }

  outTree_ -> Fill();
  endEvent();
}

//------------------------------------------------------------------------
void MuonTriggerMakeTreePFGEN::MonteCarloStudies(const edm::Event& iEvent)
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  int ZId   =    23;
  int muId  =    13;

  for ( size_t i=0; i< genParticles->size(); ++i) 
  { 
    const reco::GenParticle &p = (*genParticles)[i];
    int id = p.pdgId();
  	if(fabs(id) != ZId ) 	continue; 

    // return all daughters of Bc 
    if (p.numberOfDaughters()!=2 ) continue;  
	bool boolMuP = false;
    bool boolMuM = false;
    for ( size_t ides=0; ides < p.numberOfDaughters(); ++ides ) 
    {
  	  const reco::Candidate *des = p.daughter(ides);
	  int dauId = des->pdgId();
	  if( dauId == muId && des->status() == 1) 
	  {
	    p_mum.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
	    boolMuP = true;
      }
      else if (dauId == -1*muId && des->status() == 1)
	  {
	    p_mup.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
	    boolMuM = true;
      }
      T_GenPdgId ->push_back(des->pdgId());
      T_GenPt    ->push_back(des->pt() );
      T_GenEta   ->push_back(des->eta());
      T_GenPhi   ->push_back(des->phi());
	} // end for des
	if (!(boolMuP && boolMuM)) continue;
	boolGEN = true;
	p_z.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
  }
}



//----------------------------------------------------------------------
void MuonTriggerMakeTreePFGEN::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {

  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    assert(false);
  }
  bool changedT = true;
  if( hltConfigTag_.init(run, eventSetup, tagTriggerProcess_, changedT) ) {
  }
  else {
    std::cout << "Warning, didn't find tag process " << tagTriggerProcess_.c_str() << std::endl;
    assert(false);
  }

  triggerIndex_ = -1; 
  tagTriggerIndex_ = -1; 
  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.compare(triggerName_) == 0) {
      triggerIndex_ = int(iHltPath);
    }
    if( triggerIndex_>-1  ) break; 
  }
  for(unsigned iHltPath=0; iHltPath<hltConfigTag_.size(); ++iHltPath) {
    std::string tempName = hltConfigTag_.triggerName(iHltPath);
    if(tempName.compare(tagTriggerName_) == 0) {
      tagTriggerIndex_ = int(iHltPath);
    }
    if(tagTriggerIndex_>-1 ) break; 
  }

  if( triggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find trigger " <<  triggerName_.c_str() << std::endl;
  }
  if( tagTriggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find tag trigger " <<  tagTriggerName_.c_str() << std::endl;
  }
  
  
}


//---------------------------------------------
void MuonTriggerMakeTreePFGEN::beginEvent()
{

T_GenPdgId                = new std::vector<int>   ;
T_GenPt                   = new std::vector<double>;
T_GenEta                  = new std::vector<double>;
T_GenPhi                  = new std::vector<double>;
T_hltMatch                = new std::vector<double>;
T_hltTagMatch             = new std::vector<double>;
T_hlt_MuPt                = new std::vector<double>;
T_hlt_MuTrkPt             = new std::vector<double>;
T_hlt_MuEta               = new std::vector<double>;
T_hlt_MuPhi               = new std::vector<double>;
T_hlt_HcalClusterDep      = new std::vector<float> ;
T_hlt_HcalClusterDepV1    = new std::vector<float> ;
T_hlt_HcalClusterDepV5    = new std::vector<float> ;
T_hlt_EcalClusterDep      = new std::vector<float> ;
T_hlt_EcalClusterDepV1    = new std::vector<float> ;
T_hlt_EcalClusterDepV5    = new std::vector<float> ;
T_hlt_HcalClusterDepEt    = new std::vector<float> ;
T_hlt_HcalClusterDepV1Et  = new std::vector<float> ;
T_hlt_HcalClusterDepV5Et  = new std::vector<float> ;
T_hlt_EcalClusterDepEt    = new std::vector<float> ;
T_hlt_EcalClusterDepV1Et  = new std::vector<float> ;
T_hlt_EcalClusterDepV5Et  = new std::vector<float> ;
T_hlt_trkDep              = new std::vector<double>;
T_GenHLTMatch             = new std::vector<float> ;

}


void MuonTriggerMakeTreePFGEN::endEvent()
{

delete T_GenPdgId                ;
delete T_GenPt                   ;
delete T_GenEta                  ;
delete T_GenPhi                  ;
delete T_hltMatch                ;
delete T_hltTagMatch             ;
delete T_hlt_MuPt                ;
delete T_hlt_MuTrkPt             ;
delete T_hlt_MuEta               ;
delete T_hlt_MuPhi               ;           
delete T_hlt_HcalClusterDep      ;
delete T_hlt_HcalClusterDepV1    ;
delete T_hlt_HcalClusterDepV5    ;
delete T_hlt_EcalClusterDep      ;
delete T_hlt_EcalClusterDepV1    ;
delete T_hlt_EcalClusterDepV5    ;
delete T_hlt_HcalClusterDepEt    ;
delete T_hlt_HcalClusterDepV1Et  ;
delete T_hlt_HcalClusterDepV5Et  ;
delete T_hlt_EcalClusterDepEt    ;
delete T_hlt_EcalClusterDepV1Et  ;
delete T_hlt_EcalClusterDepV5Et  ;
delete T_hlt_trkDep              ;
delete T_GenHLTMatch             ;

}



// define this as a plug-in
DEFINE_FWK_MODULE(MuonTriggerMakeTreePFGEN);
