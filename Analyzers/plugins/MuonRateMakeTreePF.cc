/** \class MuonRateMakeTreePF
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
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"	

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"


#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

typedef std::vector< edm::Handle<trigger::TriggerFilterObjectWithRefs> > TrigFiltVect;

class MuonRateMakeTreePF : public edm::EDAnalyzer {

 public:
  MuonRateMakeTreePF(const edm::ParameterSet& cfg);
  virtual ~MuonRateMakeTreePF() {};

 private:
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  virtual void beginEvent();
  virtual void endEvent();

  /// input tag for mouns
  edm::InputTag vertexes_;
  edm::InputTag muons_;
  /// file service
  edm::Service<TFileService> outfile_;
  /// histograms
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TH2*> hist2D_;

  // Trigger process
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

  // Max number of offline muons allowed in the event
  unsigned int nMaxMuons_;

  TTree* outTree_;
  int T_Run,  T_Lumi, T_Event, T_NTrks, T_Nprim, T_TrueNI; 
  double T_hlt_rho;
  std::vector<double> *T_hltTagMatch, *T_hltMatch;
  std::vector<double> *T_hlt_MuPt, *T_hlt_MuTrkPt, *T_hlt_MuEta, *T_hlt_MuPhi, *T_hlt_trkDep;
  std::vector<float>  *T_hlt_HcalClusterDep,   *T_hlt_EcalClusterDep;
  std::vector<float>  *T_hlt_HcalClusterDepV1, *T_hlt_EcalClusterDepV1;
  std::vector<float>  *T_hlt_HcalClusterDepV5, *T_hlt_EcalClusterDepV5;
  std::vector<float>  *T_hlt_HcalClusterDepEtMin,   *T_hlt_EcalClusterDepEtMin;
  std::vector<float>  *T_hlt_HcalClusterDepV1EtMin, *T_hlt_EcalClusterDepV1EtMin;
  std::vector<float>  *T_hlt_HcalClusterDepV5EtMin, *T_hlt_EcalClusterDepV5EtMin;
  std::vector<int>    *T_hlt_missingL3;


  //Deposits tag & l3 muons
  edm::InputTag l3candTag_;
  edm::InputTag chargedDepTag_;
  edm::InputTag neutralDepTag_;
  edm::InputTag photonsDepTag_;
  edm::InputTag neutralDepTagV1_;
  edm::InputTag photonsDepTagV1_;
  edm::InputTag neutralDepTagV5_;
  edm::InputTag photonsDepTagV5_;
  edm::InputTag neutralDepEtMinTag_;
  edm::InputTag photonsDepEtMinTag_;
  edm::InputTag neutralDepV1EtMinTag_;
  edm::InputTag photonsDepV1EtMinTag_;
  edm::InputTag neutralDepV5EtMinTag_;
  edm::InputTag photonsDepV5EtMinTag_;
  edm::InputTag rhoCorrectionTag_;
  edm::InputTag VtxTag_;
  edm::InputTag AllVtxTag_;

  // Services
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<MuonDetLayerGeometry> detLayerGeometry_;
};

/// default constructor
MuonRateMakeTreePF::MuonRateMakeTreePF(const edm::ParameterSet& cfg): 
  vertexes_              (cfg.getParameter<edm::InputTag>("vertexes")), 
  muons_                 (cfg.getParameter<edm::InputTag>("muons")), 
  triggerProcess_        (cfg.getParameter<std::string>("triggerProcess")), 
  tagTriggerName_        (cfg.getParameter<std::string>("tagTriggerName")), 
  triggerName_           (cfg.getParameter<std::string>("triggerName")), 
  probeFilterDen_        (cfg.getParameter<std::string>("probeFilterDen")),
  probeFilterNum_        (cfg.getParameter<std::string>("probeFilterNum")),
  nMaxMuons_             (cfg.getUntrackedParameter<unsigned int>("maxNumberMuons", 999999)),
  l3candTag_             (cfg.getUntrackedParameter<edm::InputTag>("L3CandidatesLabel")),
  chargedDepTag_         (cfg.getUntrackedParameter<edm::InputTag>("ChargedDepositLabel")), 
  neutralDepTag_         (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabel")), 
  photonsDepTag_         (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabel")), 
  neutralDepTagV1_       (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelV1")), 
  photonsDepTagV1_       (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelV1")), 
  neutralDepTagV5_       (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelV5")), 
  photonsDepTagV5_       (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelV5")), 
  neutralDepEtMinTag_    (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelEtMin")), 
  photonsDepEtMinTag_    (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelEtMin")), 
  neutralDepV1EtMinTag_  (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelV1EtMin")), 
  photonsDepV1EtMinTag_  (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelV1EtMin")), 
  neutralDepV5EtMinTag_  (cfg.getUntrackedParameter<edm::InputTag>("NeutralDepositLabelV5EtMin")), 
  photonsDepV5EtMinTag_  (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDepositLabelV5EtMin")), 
  rhoCorrectionTag_      (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionLabel")), 
  VtxTag_                (cfg.getUntrackedParameter<edm::InputTag>("VtxLabel")), 
  AllVtxTag_             (cfg.getUntrackedParameter<edm::InputTag>("AllVtxLabel")) 
{}

void MuonRateMakeTreePF::beginJob() {

  TH1::SetDefaultSumw2() ;

  outTree_ = outfile_->make<TTree>("outTree_", "outTree_");

  outTree_->Branch("Run"                           , &T_Run                         );
  outTree_->Branch("Lumi"                          , &T_Lumi                        );
  outTree_->Branch("Event"                         , &T_Event                       );
  outTree_->Branch("Nprim"                         , &T_Nprim                       );
  outTree_->Branch("TrueNI"                        , &T_TrueNI                      );
            
  outTree_->Branch("T_hlt_MuPt"                    , &T_hlt_MuPt                    );
  outTree_->Branch("T_hlt_MuTrkPt"                 , &T_hlt_MuTrkPt                 );
  outTree_->Branch("T_hlt_MuEta"                   , &T_hlt_MuEta                   );
  outTree_->Branch("T_hlt_MuPhi"                   , &T_hlt_MuPhi                   );
  outTree_->Branch("T_hlt_trkDep"                  , &T_hlt_trkDep                  );
  outTree_->Branch("T_hlt_EcalClusterDep"          , &T_hlt_EcalClusterDep          );
  outTree_->Branch("T_hlt_HcalClusterDep"          , &T_hlt_HcalClusterDep          );
  outTree_->Branch("T_hlt_EcalClusterDepV1"        , &T_hlt_EcalClusterDepV1        );
  outTree_->Branch("T_hlt_HcalClusterDepV1"        , &T_hlt_HcalClusterDepV1        );
  outTree_->Branch("T_hlt_EcalClusterDepV5"        , &T_hlt_EcalClusterDepV5        );
  outTree_->Branch("T_hlt_HcalClusterDepV5"        , &T_hlt_HcalClusterDepV5        );

  outTree_->Branch("T_hlt_EcalClusterDepEtMin"     , &T_hlt_EcalClusterDepEtMin     );
  outTree_->Branch("T_hlt_HcalClusterDepEtMin"     , &T_hlt_HcalClusterDepEtMin     );
  outTree_->Branch("T_hlt_EcalClusterDepV1EtMin"   , &T_hlt_EcalClusterDepV1EtMin   );
  outTree_->Branch("T_hlt_HcalClusterDepV1EtMin"   , &T_hlt_HcalClusterDepV1EtMin   );
  outTree_->Branch("T_hlt_EcalClusterDepV5EtMin"   , &T_hlt_EcalClusterDepV5EtMin   );
  outTree_->Branch("T_hlt_HcalClusterDepV5EtMin"   , &T_hlt_HcalClusterDepV5EtMin   );

  outTree_->Branch("T_hlt_rho"                     , &T_hlt_rho                     );

  hists_["Mu24_sizeToc"    ]   = outfile_->make<TH1F>("Mu24_sizeToc"    , "Mu24_sizeToc"     ,  30,  0., 30.);
  hists_["IsoMu24_sizeToc" ]   = outfile_->make<TH1F>("IsoMu24_sizeToc" , "IsoMu24_sizeToc"  ,  30,  0., 30.);

  hists_["Mu24_pt"   ]      = outfile_->make<TH1F>("Mu24_pt"     , "pt"  ,  100,  0., 300.);
  hists_["IsoMu24_pt"]      = outfile_->make<TH1F>("IsoMu24_pt"  , "pt"  ,  100,  0., 300.);
}

void MuonRateMakeTreePF::endJob() {
}


void MuonRateMakeTreePF::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonRateMakeTreePF::analyze(const edm::Event &event, const edm::EventSetup &eventSetup) {

  using reco::Muon;

  beginEvent();

  // Get trigger results
  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(edm::InputTag("TriggerResults", "", triggerProcess_), triggerResults);

  if(!triggerResults.isValid()) {
    std::cout << "Trigger results not valid" << std::endl;
    return;
  } 

  if( !triggerResults->accept(tagTriggerIndex_) ) return; // there are no tags

  // Get trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent);

  if(!triggerEvent.isValid()) { 
    std::cout << "TriggerEvent not valid" << std::endl;
    return;
  }

  // Sanity check
  assert(triggerResults->size()==hltConfig_.size());

  TrigFiltVect hltFilters( 1 );
  event.getByLabel(probeFilterDen_, hltFilters[0]);


  // Handle to the rho collection
  edm::Handle <double>  rhoCollection;
  event.getByLabel(rhoCorrectionTag_, rhoCollection);
  T_Run    = event.id().run();
  T_Lumi   = event.id().luminosityBlock();
  T_Event  = event.id().event();
  if (rhoCollection.isValid()) T_hlt_rho = *(rhoCollection.product());

  // Handle to the online deposits 
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV1;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV1;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV5;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV5;

  event.getByLabel(neutralDepTag_, neutralDepMap);
  event.getByLabel(photonsDepTag_, photonsDepMap);
  event.getByLabel(neutralDepTagV1_, neutralDepMapV1);
  event.getByLabel(photonsDepTagV1_, photonsDepMapV1);
  event.getByLabel(neutralDepTagV5_, neutralDepMapV5);
  event.getByLabel(photonsDepTagV5_, photonsDepMapV5);

  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapEtMin;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapEtMin;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV1EtMin;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV1EtMin;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMapV5EtMin;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMapV5EtMin;

  event.getByLabel(neutralDepEtMinTag_, neutralDepMapEtMin);
  event.getByLabel(photonsDepEtMinTag_, photonsDepMapEtMin);
  event.getByLabel(neutralDepV1EtMinTag_, neutralDepMapV1EtMin);
  event.getByLabel(photonsDepV1EtMinTag_, photonsDepMapV1EtMin);
  event.getByLabel(neutralDepV5EtMinTag_, neutralDepMapV5EtMin);
  event.getByLabel(photonsDepV5EtMinTag_, photonsDepMapV5EtMin);

  edm::Handle<reco::IsoDepositMap> TrkDepMap;
  event.getByLabel(chargedDepTag_, TrkDepMap);

  // Get trigger objects from trigger summary
  const trigger::TriggerObjectCollection & toc = triggerEvent->getObjects();

  // Modules in tag trigger path
  const std::vector<std::string>& tagModuleLabels(hltConfig_.moduleLabels(tagTriggerIndex_));
  assert( tagModuleLabels.size()==hltConfig_.size(tagTriggerIndex_) );
  const unsigned int tagModuleIndex( hltConfig_.size(tagTriggerIndex_)-2 ); // index of last filter (excluding HLTEndBool)
  const unsigned int tagFilterIndex( triggerEvent->filterIndex( edm::InputTag( tagModuleLabels[tagModuleIndex], "", triggerProcess_) ) );
  assert( tagFilterIndex < triggerEvent->sizeFilters() );
  const trigger::Vids & tagVids( triggerEvent->filterIds(tagFilterIndex) );
  const trigger::Keys & tagKeys( triggerEvent->filterKeys(tagFilterIndex) );
  assert( tagVids.size()==tagKeys.size() );
  const unsigned int nTagTrig(tagVids.size());

  hists_["Mu24_sizeToc"]->Fill(nTagTrig); 
//   const trigger::TriggerObject & tagTo = toc[tagKeys[0]];
//   hists_["Mu24_pt"]->Fill(tagTo.pt()); 

  for(unsigned int i=0; i!=nTagTrig; ++i) 
  {
	const trigger::TriggerObject & tagTo = toc[tagKeys[i]];
	hists_["Mu24_pt"]->Fill(tagTo.pt()); 
  }


  // Modules in probe trigger path
  const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_));
  const unsigned int m(hltConfig_.size(triggerIndex_));
  assert( moduleLabels.size()==m );
  const unsigned int lastModuleIndex(triggerResults->index(triggerIndex_));
  assert(lastModuleIndex<m);

  unsigned int filterIndex = triggerEvent->sizeFilters();

  for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
	if( probeFilterDen_.compare(moduleLabels[j])!=0 ) continue; 
	filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterNum_, "", triggerProcess_));
	break; 
  }


  std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
  reco::RecoChargedCandidateRef ref;
  hltFilters[0]->getObjects(trigger::TriggerMuon, prevMuonRefs);
  

//   ref = prevMuonRefs[j];
//   reco::TrackRef mu = ref->track();
  if( filterIndex<triggerEvent->sizeFilters() ) 
  {
	const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
	const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
	assert( vids.size()==keys.size() );
	const unsigned int nProbeTrig(vids.size());
	hists_["IsoMu24_sizeToc"]->Fill(nProbeTrig); 

    if (nProbeTrig != prevMuonRefs.size()) std::cout << "different sizes" << std::endl;


	for(unsigned int i=0; i<nProbeTrig; ++i) 
    {
	  const trigger::TriggerObject & to = toc[keys[i]];
	  hists_["IsoMu24_pt"]->Fill(to.pt()); 

      T_hlt_MuPt ->push_back(to.pt()); 
      T_hlt_MuEta->push_back(to.eta()); 
      T_hlt_MuPhi->push_back(to.phi()); 

      ref = prevMuonRefs[i];
      reco::TrackRef trkmu = ref->track();
      T_hlt_MuTrkPt ->push_back(trkmu->pt());
//       reco::TrackRef mu = ref->track();
      if (fabs(to.pt() - ref->pt()) > 0.0001) std::cout << "different pt: " << to.pt() << " != " <<  ref->pt() << std::endl;

      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi = (*neutralDepMap).find( ref );
      T_hlt_HcalClusterDep -> push_back(hcal_mapi->val);
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi = (*photonsDepMap).find( ref );
      T_hlt_EcalClusterDep -> push_back(ecal_mapi->val); 

      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV1 = (*neutralDepMapV1).find( ref );
      T_hlt_HcalClusterDepV1 -> push_back(hcal_mapiV1->val);
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV1 = (*photonsDepMapV1).find( ref );
      T_hlt_EcalClusterDepV1 -> push_back(ecal_mapiV1->val); 

      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV5 = (*neutralDepMapV5).find( ref );
      T_hlt_HcalClusterDepV5 -> push_back(hcal_mapiV5->val);
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV5 = (*photonsDepMapV5).find( ref );
      T_hlt_EcalClusterDepV5 -> push_back(ecal_mapiV5->val); 


      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi_et = (*neutralDepMapEtMin).find( ref );
      T_hlt_HcalClusterDepEtMin -> push_back(hcal_mapi_et->val);
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi_et = (*photonsDepMapEtMin).find( ref );
      T_hlt_EcalClusterDepEtMin -> push_back(ecal_mapi_et->val); 

      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV1_et = (*neutralDepMapV1EtMin).find( ref );
      T_hlt_HcalClusterDepV1EtMin -> push_back(hcal_mapiV1_et->val);
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV1_et = (*photonsDepMapV1EtMin).find( ref );
      T_hlt_EcalClusterDepV1EtMin -> push_back(ecal_mapiV1_et->val); 

      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapiV5_et = (*neutralDepMapV5EtMin).find( ref );
      T_hlt_HcalClusterDepV5EtMin -> push_back(hcal_mapiV5_et->val);
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapiV5_et = (*photonsDepMapV5EtMin).find( ref );
      T_hlt_EcalClusterDepV5EtMin -> push_back(ecal_mapiV5_et->val); 



	  reco::IsoDeposit theTkIsolation = (*TrkDepMap)[ref];
	  T_hlt_trkDep->push_back(theTkIsolation.depositWithin(0.3));


	}
  } // end if( filterIndex<triggerEvent->sizeFilters() )

  outTree_ -> Fill();
  endEvent();

}


//---------------------------------------------
void MuonRateMakeTreePF::beginEvent()
{

T_hlt_MuPt                    = new std::vector<double>;
T_hlt_MuTrkPt                 = new std::vector<double>;
T_hlt_MuEta                   = new std::vector<double>;
T_hlt_MuPhi                   = new std::vector<double>;
T_hlt_HcalClusterDep          = new std::vector<float>;
T_hlt_EcalClusterDep          = new std::vector<float>;
T_hlt_HcalClusterDepV1        = new std::vector<float>;
T_hlt_EcalClusterDepV1        = new std::vector<float>;
T_hlt_HcalClusterDepV5        = new std::vector<float>;
T_hlt_EcalClusterDepV5        = new std::vector<float>;
T_hlt_HcalClusterDepEtMin     = new std::vector<float>;
T_hlt_EcalClusterDepEtMin     = new std::vector<float>;
T_hlt_HcalClusterDepV1EtMin   = new std::vector<float>;
T_hlt_EcalClusterDepV1EtMin   = new std::vector<float>;
T_hlt_HcalClusterDepV5EtMin   = new std::vector<float>;
T_hlt_EcalClusterDepV5EtMin   = new std::vector<float>;
T_hlt_trkDep                  = new std::vector<double>;

}


void MuonRateMakeTreePF::endEvent()
{

delete T_hlt_MuPt                  ;
delete T_hlt_MuTrkPt               ;
delete T_hlt_MuEta                 ;
delete T_hlt_MuPhi                 ;
delete T_hlt_HcalClusterDep        ;
delete T_hlt_EcalClusterDep        ;
delete T_hlt_HcalClusterDepV1      ;
delete T_hlt_EcalClusterDepV1      ;
delete T_hlt_HcalClusterDepV5      ;
delete T_hlt_EcalClusterDepV5      ;
delete T_hlt_HcalClusterDepEtMin   ;
delete T_hlt_EcalClusterDepEtMin   ;
delete T_hlt_HcalClusterDepV1EtMin ;
delete T_hlt_EcalClusterDepV1EtMin ;
delete T_hlt_HcalClusterDepV5EtMin ;
delete T_hlt_EcalClusterDepV5EtMin ;
delete T_hlt_trkDep                ;

}


//---------------------------------------------
void MuonRateMakeTreePF::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {
  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    // Now crash
    assert(false);
  }

  triggerIndex_ = -1; 
  tagTriggerIndex_ = -1; 

  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.compare(triggerName_) == 0) {
      triggerIndex_ = int(iHltPath);
    }
    if(tempName.compare(tagTriggerName_) == 0) {
      tagTriggerIndex_ = int(iHltPath);
    }
    if( triggerIndex_>-1 && tagTriggerIndex_>-1 ) break; 
  } // end for each path

  if( triggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find trigger " <<  triggerName_.c_str() << std::endl;
//     assert(false);    
  }
  if( tagTriggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find tag trigger " <<  tagTriggerName_.c_str() << std::endl;
//     assert(false);    
  }
}


// define this as a plug-in
DEFINE_FWK_MODULE(MuonRateMakeTreePF);
