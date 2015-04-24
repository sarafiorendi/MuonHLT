/** \class MuonTriggerMakeTree
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
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "TLorentzVector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class MuonTriggerMakeTree : public edm::EDAnalyzer {

 public:
  /// default constructor
  MuonTriggerMakeTree(const edm::ParameterSet& cfg);
  /// default destructor
  virtual ~MuonTriggerMakeTree() {};
  /// everything that needs to be done before the event loop

 private:
  virtual void beginJob();
  /// everything that needs to be done after the event loop
  virtual void endJob();
  /// everything that needs to be done before each run
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done after each run
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  virtual void beginEvent();
  virtual void endEvent();

  /// input tag for mouns
  edm::InputTag vertexes_;
  /// input tag for mouns
  edm::InputTag muons_;
  /// file service
  edm::Service<TFileService> outfile_;
  /// histograms
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TH2*> hist2D_;

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
  // Max number of offline muons allowed in the event
  unsigned int nMaxMuons_;
  
  TTree* outTree_;
  int T_Run,  T_Lumi, T_Event, T_NTrks, T_Nprim, T_TrueNI; 
  double T_hlt_rho, T_Rho;
  //muon reco
  std::vector<double> *T_MuPt,  *T_MuEta,  *T_MuPhi, *T_MuEcalDep, *T_MuHcalDep, *T_MuPUDep, *T_MuTrkDep;
  std::vector<double> *T_GenPt, *T_GenEta, *T_GenPhi;
  std::vector<double> *T_hltTagMatch, *T_hltMatch;
  std::vector<double> *T_hlt_MuPt, *T_hlt_MuEta, *T_hlt_MuPhi, *T_hlt_MuTrkPt;
  std::vector<double> *T_hlt_ECalDep, *T_hlt_HCalDep, *T_hlt_TotIso, *T_hlt_CaloDep, *T_hlt_TrkDep;
  std::vector<float>  *T_GenRecoMatch;
  std::vector<int>    *T_hlt_missingL3, *T_IsTight;
  

  //Deposits tag & l3 muons
  edm::InputTag l3candTag_;
  edm::InputTag caloTag_;
  edm::InputTag ecalTag_;
  edm::InputTag hcalTag_;
  edm::InputTag trkisoTag_;
  edm::InputTag isoTag_;
  edm::InputTag rhoCorrectionTag_;
  edm::InputTag rhoCorrectionOfflineTag_;
  edm::InputTag VtxTag_;
  edm::InputTag AllVtxTag_;
  
  //Cuts
  double offlineIsoCut_;
  //Debug 
  bool debug_;

  // Services
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<MuonDetLayerGeometry> detLayerGeometry_;
};

/// default constructor
MuonTriggerMakeTree::MuonTriggerMakeTree(const edm::ParameterSet& cfg): 
  vertexes_         (cfg.getParameter<edm::InputTag>("vertexes")), 
  muons_            (cfg.getParameter<edm::InputTag>("muons")), 
  tagTriggerProcess_(cfg.getParameter<std::string>("tagTriggerProcess")), 
  triggerProcess_   (cfg.getParameter<std::string>("triggerProcess")), 
  tagTriggerName_   (cfg.getParameter<std::string>("tagTriggerName")), 
  triggerName_      (cfg.getParameter<std::string>("triggerName")), 
  probeFilterDen_   (cfg.getParameter<std::string>("probeFilterDen")),
  probeFilterNum_   (cfg.getParameter<std::string>("probeFilterNum")),
  nMaxMuons_        (cfg.getUntrackedParameter<unsigned int>("maxNumberMuons", 999999)),
  l3candTag_        (cfg.getUntrackedParameter<edm::InputTag>("L3CandidatesLabel"      )),
  caloTag_          (cfg.getUntrackedParameter<edm::InputTag>("CaloDepositLabel")), 
  ecalTag_          (cfg.getUntrackedParameter<edm::InputTag>("ECalDepositLabel")), 
  hcalTag_          (cfg.getUntrackedParameter<edm::InputTag>("HCalDepositLabel")), 
  trkisoTag_        (cfg.getUntrackedParameter<edm::InputTag>("TrkIsoLabel")), 
  isoTag_           (cfg.getUntrackedParameter<edm::InputTag>("IsoDepositLabel")), 
  rhoCorrectionTag_ (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionLabel")), 
  rhoCorrectionOfflineTag_(cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionOfflineLabel")), 
  VtxTag_           (cfg.getUntrackedParameter<edm::InputTag>("VtxLabel")), 
  offlineIsoCut_    (cfg.getUntrackedParameter<double>("OfflineIsoCut")), 
  debug_            (cfg.getUntrackedParameter<bool>("boolDebug")) 
{}

void MuonTriggerMakeTree::beginJob() {

  TH1::SetDefaultSumw2() ;

  double eta_bins[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4}; 
  int eta_bin_n = sizeof(eta_bins)/sizeof(double); 
  
  outTree_ = outfile_->make<TTree>("outTree_", "outTree_");

  outTree_->Branch("Run"                           , &T_Run                         );
  outTree_->Branch("Lumi"                          , &T_Lumi                        );
  outTree_->Branch("Event"                         , &T_Event                       );
  outTree_->Branch("Nprim"                         , &T_Nprim                       );
  outTree_->Branch("TrueNI"                        , &T_TrueNI                      );
  outTree_->Branch("T_Rho"                         , &T_Rho                         );
            
  outTree_->Branch("T_IsTight"                     , &T_IsTight                     );
  outTree_->Branch("T_MuPt"                        , &T_MuPt                        );
  outTree_->Branch("T_MuEta"                       , &T_MuEta                       );
  outTree_->Branch("T_MuPhi"                       , &T_MuPhi                       );
  outTree_->Branch("T_MuEcalDep"                   , &T_MuEcalDep                   );
  outTree_->Branch("T_MuHcalDep"                   , &T_MuHcalDep                   );
  outTree_->Branch("T_MuTrkDep"                    , &T_MuTrkDep                    );
  outTree_->Branch("T_MuPUDep"                     , &T_MuPUDep                     );
  
  outTree_->Branch("T_GenPt"                       , &T_GenPt                       );
  outTree_->Branch("T_GenEta"                      , &T_GenEta                      );
  outTree_->Branch("T_GenPhi"                      , &T_GenPhi                      );
  outTree_->Branch("T_GenRecoMatch"                , &T_GenRecoMatch                );

  outTree_->Branch("T_hlt_MuPt"                    , &T_hlt_MuPt                    );
  outTree_->Branch("T_hlt_MuTrkPt"                 , &T_hlt_MuTrkPt                 );
  outTree_->Branch("T_hlt_MuEta"                   , &T_hlt_MuEta                   );
  outTree_->Branch("T_hlt_MuPhi"                   , &T_hlt_MuPhi                   );
  outTree_->Branch("T_hlt_ECalDep"                 , &T_hlt_ECalDep                 );
  outTree_->Branch("T_hlt_HCalDep"                 , &T_hlt_HCalDep                 );
  outTree_->Branch("T_hlt_CaloDep"                 , &T_hlt_CaloDep                 );
  outTree_->Branch("T_hlt_TrkDep"                  , &T_hlt_TrkDep                  );
  outTree_->Branch("T_hlt_TotIso"                  , &T_hlt_TotIso                  );
  outTree_->Branch("T_hlt_rho"                     , &T_hlt_rho                     );
  outTree_->Branch("T_hltMatch"                    , &T_hltMatch                    );
  outTree_->Branch("T_hltTagMatch"                 , &T_hltTagMatch                 );
  outTree_->Branch("T_hlt_missingL3"               , &T_hlt_missingL3               );
  
  
  hists_["muonPt_offline"   ] = outfile_->make<TH1F>("muonPt_offline"   , "pt"     , 100        ,  0., 300.);
  hists_["muonEta_offline"  ] = outfile_->make<TH1F>("muonEta_offline"  , "eta"    , eta_bin_n-1, eta_bins );
  hists_["muonPhi_offline"  ] = outfile_->make<TH1F>("muonPhi_offline"  , "phi"    , 100        , -5.,   5.); 
  hists_["muonNvtx_offline" ] = outfile_->make<TH1F>("muonNvtx_offline" , "nvtx"   , 100        , 0.5, 100 );

  hists_["deltaR_trobj_probe"	 ] = outfile_->make<TH1F>("deltaR_trobj_probe" , "#DeltaR(trig,#mu)" , 600, 0., 6.0); 

}    

void MuonTriggerMakeTree::endJob() {
}


void MuonTriggerMakeTree::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonTriggerMakeTree::analyze(const edm::Event &event, const edm::EventSetup &eventSetup) {

  using reco::Muon;
  
  beginEvent();

  edm::Handle<reco::VertexCollection> pvHandle; 
  event.getByLabel(vertexes_, pvHandle);
  const reco::VertexCollection & vertices = *pvHandle.product();
  unsigned int nGoodVtx = 0; 
  for(reco::VertexCollection::const_iterator it=vertices.begin(); it!=vertices.end(); ++it) {
    if( it->ndof()>4                     && 
	(fabs(it->z())<=24.)             && 
	(fabs(it->position().rho())<=2.)   ) 
      nGoodVtx++;
  }
  if( nGoodVtx==0 ) return;
  const reco::Vertex & pv = vertices[0];

  // Handle to the rho collection
  edm::Handle <double>  rhoCollection;
  event.getByLabel(rhoCorrectionTag_, rhoCollection);

  // Handle to the offline rho collection
  edm::Handle <double>  rhoCollectionOffline;
  event.getByLabel(rhoCorrectionOfflineTag_, rhoCollectionOffline);

  T_Run    = event.id().run();
  T_Lumi   = event.id().luminosityBlock();
  T_Event  = event.id().event();
  T_Nprim  = vertices.size();
  
  if (rhoCollection.isValid()) T_hlt_rho = *(rhoCollection.product());
  if (rhoCollectionOffline.isValid()) T_Rho     = *(rhoCollectionOffline.product());


//   float truePu=0.;
  edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
  try {
    event.getByLabel("addPileupInfo",puInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    //The in-time crossing is getBunchCrossing = 0; negative ones are early, positive ones are late.
    for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) 
    {
      if(PVI->getBunchCrossing()==0){
//         T_Event_nPU = PVI->getPU_NumInteractions();
        T_TrueNI    = PVI->getTrueNumInteractions();
      }
    }
  } catch (...) {}


  // Handle to the muon collection
  edm::Handle<std::vector<Muon> > muons;
  event.getByLabel(muons_, muons);
  if( nMaxMuons_>0 && muons->size()>nMaxMuons_ ) return; 


  // Handle to the online deposits 
  edm::Handle <reco::IsoDepositMap > caloDepMap;       
  event.getByLabel(caloTag_, caloDepMap);

  edm::Handle<reco::IsoDepositMap> IsoCaloMapEcal;
  event.getByLabel(ecalTag_,IsoCaloMapEcal);

  edm::Handle<reco::IsoDepositMap> IsoCaloMapHcal;
  event.getByLabel(hcalTag_,IsoCaloMapHcal);

  edm::Handle<reco::IsoDepositMap> IsoMapTrk;
  event.getByLabel(trkisoTag_, IsoMapTrk);

  edm::Handle < edm::ValueMap<double>> isoDepMap;       
  event.getByLabel(isoTag_, isoDepMap);

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

  if( !tagTriggerResults->accept(tagTriggerIndex_) ) return; // there are no tags


  // GEN particles	
  MonteCarloStudies(event);
  if (!boolGEN) return;

  // Get trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent);
  if(!triggerEvent.isValid()) { 
    std::cout << "TriggerEvent not valid" << std::endl;
    return;
  }

  // Sanity check
  assert(triggerResults->size()==hltConfig_.size());

  // Get trigger objects from trigger summary
  const trigger::TriggerObjectCollection & toc = triggerEvent->getObjects();

  // Loop muon collection and fill histograms
  for(std::vector<Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) 
  { 
    if( muon::isTightMuon( (*mu1), pv ) ) T_IsTight -> push_back(1);
    else                                  T_IsTight -> push_back(0);
//     if(!(muon::isLooseMuon(*mu1))       ) continue; 
    if(! ((*mu1).pt()>25)               ) continue;
    if( fabs((*mu1).eta())  > 2.1       ) continue;

    T_MuPt ->push_back((*mu1).pt());
    T_MuEta->push_back((*mu1).eta());
    T_MuPhi->push_back((*mu1).phi());
    T_MuTrkDep  -> push_back(mu1->pfIsolationR03().sumChargedHadronPt);
    T_MuHcalDep -> push_back(mu1->pfIsolationR03().sumNeutralHadronEt);
    T_MuEcalDep -> push_back(mu1->pfIsolationR03().sumPhotonEt       );
    T_MuPUDep   -> push_back(mu1->pfIsolationR03().sumPUPt           );

   //Define offline isolation for tag muon
// 	float neutralDeposits = mu1->pfIsolationR04().sumNeutralHadronEt + mu1->pfIsolationR04().sumPhotonEt - 0.5*mu1->pfIsolationR04().sumPUPt; 
// 	if (neutralDeposits < 0)  neutralDeposits = 0 ;
// 	float m1PFiso = (mu1->pfIsolationR04().sumChargedHadronPt + neutralDeposits)/mu1->pt();

    //Match to GEN
    float matchMu = 10000;
    if ( mu1->charge() > 0 )    matchMu = deltaR(mu1->eta(),mu1->phi(),p_mup.Eta(),p_mup.Phi());
	else				        matchMu = deltaR(mu1->eta(),mu1->phi(),p_mum.Eta(),p_mum.Phi());
	T_GenRecoMatch->push_back(matchMu);

	hists_["muonPt_offline"  ]->Fill( mu1->pt () );
	hists_["muonEta_offline" ]->Fill( mu1->eta() );
	hists_["muonPhi_offline" ]->Fill( mu1->phi() );
	hists_["muonNvtx_offline"]->Fill( nGoodVtx   );

	// Check if it is matched to the last-1 hlt filter
	bool isGood = false; 
	const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_));
	const unsigned int m(hltConfig_.size(triggerIndex_));
	assert( moduleLabels.size()==m );
	const unsigned int lastModuleIndex(triggerResults->index(triggerIndex_));
	assert(lastModuleIndex<m);

    double muEta = (*mu1).eta();
    double muPhi = (*mu1).phi();
    double tmpProbeDeltaR;
	//If no filter is specified for the denominator, take all offline muons
	if( probeFilterDen_.length()==0 ) isGood = true; 
	else {   
	  unsigned int filterIndex = triggerEvent->sizeFilters();
	  double maxProbeDeltaR = 100;
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
		  tmpProbeDeltaR = deltaR( muEta, muPhi, to.eta(), to.phi() ); 
		  if( tmpProbeDeltaR < maxProbeDeltaR ) {
			isGood = true;
			break;
		  }
		}
	  } // end if( filterIndex<triggerEvent->sizeFilters() )
	} // end else

    T_hltTagMatch->push_back(tmpProbeDeltaR);
	if(!isGood) continue; 

	// Check if probe passes
	unsigned int filterIndex = triggerEvent->sizeFilters(); 
	double finProbeDeltaR = 999999.; 

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

	  edm::Handle<reco::RecoChargedCandidateCollection> l3cands;
	  event.getByLabel(l3candTag_, l3cands);
	  int thisL3 = -1; 
	  double minL3DeltaR = 0.1; 
	  for(unsigned int il3=0; il3<l3cands->size(); ++il3) {
	    double tmpL3DeltaR = deltaR( muEta, muPhi, l3cands->at(il3).eta(), l3cands->at(il3).phi() ); 
	    if(tmpL3DeltaR<minL3DeltaR) {
		  minL3DeltaR = tmpL3DeltaR; 
		  thisL3 = il3; 
	    }
	  }
	  if(thisL3>-1) {
	    reco::RecoChargedCandidateRef candref(l3cands, (unsigned int)thisL3);
        T_hlt_MuPt ->push_back(candref->pt()); 
        T_hlt_MuEta->push_back(candref->eta()); 
        T_hlt_MuPhi->push_back(candref->phi()); 
        reco::TrackRef trkmu = candref->track();
	    T_hlt_MuTrkPt ->push_back(trkmu->pt());
	    
        if (caloDepMap.isValid()){
          reco::IsoDeposit theIsolation = (*caloDepMap)[candref];
          T_hlt_CaloDep ->push_back( theIsolation.depositWithin(0.3));
        }
        else T_hlt_CaloDep ->push_back( -10);

        if (IsoCaloMapEcal.isValid()){
          reco::IsoDeposit theEcalIsolation = (*IsoCaloMapEcal)[candref];
          T_hlt_ECalDep ->push_back( theEcalIsolation.depositWithin(0.3));
        }
        else T_hlt_ECalDep ->push_back( -10);
       
        if (IsoCaloMapHcal.isValid()){
          reco::IsoDeposit theHcalIsolation = (*IsoCaloMapHcal)[candref];
          T_hlt_HCalDep ->push_back( theHcalIsolation.depositWithin(0.3));
        }
        else T_hlt_HCalDep ->push_back( -10);

        if (IsoMapTrk.isValid()){
          reco::IsoDeposit theTrkIsolation = (*IsoMapTrk)[candref];
          T_hlt_TrkDep ->push_back( theTrkIsolation.depositWithin(0.3));
        }
        else T_hlt_TrkDep ->push_back( -10);

	    const edm::ValueMap<double>::value_type & isoDep = (*isoDepMap)[candref];
	    T_hlt_TotIso->push_back(isoDep);

        T_hlt_missingL3->push_back(0);

// 	    const edm::ValueMap<double>::value_type & isohlt = (*isoDepMap)[candref];
// 	    
// 	    reco::IsoDeposit theTkIsolation = (*TrkDepMap)[candref];
//         T_hlt_trkDep->push_back(theTkIsolation.depositWithin(0.3));
                
//         reco::TrackRef muref = ref->track();
//         trkDeps[j] = trkExtractor->deposit(event, eventSetup, *muref);
// 

// 		  const edm::ValueMap<double>::value_type & chargedDep = (*chargedDepMap)[candref];
// 		  const edm::ValueMap<double>::value_type & neutralDep = (*neutralDepMap)[candref];
// 		  const edm::ValueMap<double>::value_type & photonsDep = (*photonsDepMap)[candref];

	  }
	  else {
	    T_hlt_missingL3->push_back(1);
        T_hlt_MuPt    ->push_back(-100); 
        T_hlt_MuTrkPt ->push_back(-100); 
        T_hlt_MuEta   ->push_back(-100); 
        T_hlt_MuPhi   ->push_back(-100); 
	    T_hlt_CaloDep ->push_back(-100);
	    T_hlt_TotIso  ->push_back(-100);
        T_hlt_TrkDep  ->push_back(-100);
        T_hlt_ECalDep ->push_back(-100);
        T_hlt_HCalDep ->push_back(-100);
	  }
	  for(unsigned int i=0; i<nProbeTrig; ++i) {
		finProbeDeltaR = 999999.;
		const trigger::TriggerObject & to = toc[keys[i]];
		double tmpProbeDeltaR = deltaR( muEta, muPhi, to.eta(), to.phi() ); 
		if( tmpProbeDeltaR<finProbeDeltaR ) {
		  finProbeDeltaR = tmpProbeDeltaR;
		}
	  }
	  T_hltMatch->push_back(tmpProbeDeltaR);
	} // end if( filterIndex<triggerEvent->sizeFilters() )

	hists_["deltaR_trobj_probe"]->Fill(finProbeDeltaR); 


  }
  outTree_ -> Fill();
  endEvent();
}

//------------------------------------------------------------------------
void MuonTriggerMakeTree::MonteCarloStudies(const edm::Event& iEvent)
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
	  if( dauId == muId ) 
	  {
	    p_mum.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
	    boolMuP = true;
      }
      else if (dauId == -1*muId)
	  {
	    p_mup.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
	    boolMuM = true;
      }
      T_GenPt ->push_back(des->pt() );
      T_GenEta->push_back(des->eta());
      T_GenPhi->push_back(des->phi());
	} // end for des
	if (!(boolMuP && boolMuM)) continue;
	boolGEN = true;
	p_z.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
  }
}



//----------------------------------------------------------------------
void MuonTriggerMakeTree::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {

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
void MuonTriggerMakeTree::beginEvent()
{

T_MuPt         = new std::vector<double>;
T_MuEta        = new std::vector<double>;
T_MuPhi        = new std::vector<double>;
T_MuTrkDep     = new std::vector<double>;
T_MuEcalDep    = new std::vector<double>;
T_MuHcalDep    = new std::vector<double>;
T_MuPUDep      = new std::vector<double>;
T_GenPt        = new std::vector<double>;
T_GenEta       = new std::vector<double>;
T_GenPhi       = new std::vector<double>;
T_hltTagMatch  = new std::vector<double>;
T_hltMatch     = new std::vector<double>;
T_hlt_MuPt     = new std::vector<double>;
T_hlt_MuTrkPt  = new std::vector<double>;
T_hlt_MuEta    = new std::vector<double>;
T_hlt_MuPhi    = new std::vector<double>;
T_hlt_TotIso   = new std::vector<double>;
T_hlt_CaloDep  = new std::vector<double>;
T_hlt_TrkDep   = new std::vector<double>;
T_hlt_ECalDep  = new std::vector<double>;
T_hlt_HCalDep  = new std::vector<double>;
T_GenRecoMatch = new std::vector<float> ;
T_hlt_missingL3= new std::vector<int>   ;
T_IsTight      = new std::vector<int>   ;

}


void MuonTriggerMakeTree::endEvent()
{

delete T_MuPt         ;
delete T_MuEta        ;
delete T_MuPhi        ;
delete T_IsTight      ;
delete T_MuEcalDep    ;
delete T_MuHcalDep    ;
delete T_MuTrkDep     ;
delete T_MuPUDep      ;
delete T_GenPt        ;
delete T_GenEta       ;
delete T_GenPhi       ;
delete T_hltTagMatch  ;
delete T_hltMatch     ;
delete T_hlt_MuPt     ;
delete T_hlt_MuTrkPt  ;
delete T_hlt_MuEta    ;
delete T_hlt_MuPhi    ;
delete T_hlt_TotIso   ;
delete T_hlt_CaloDep  ;
delete T_hlt_ECalDep  ;
delete T_hlt_HCalDep  ;
delete T_hlt_TrkDep   ;
delete T_GenRecoMatch ;
delete T_hlt_missingL3;

}


// define this as a plug-in
DEFINE_FWK_MODULE(MuonTriggerMakeTree);
