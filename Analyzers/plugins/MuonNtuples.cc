/** \class MuonNtuples
 */
      
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

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
               const edm::Event &,
               bool 
              );
  
  void fillMuons(const edm::Handle<reco::MuonCollection> &,
                 const reco::Vertex &, 
                 const edm::Event   & 
                );

  void fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> &,
                    const edm::Event   &,
                    bool isL3           ,
                    bool isIterL3       ,
                    bool isTk           
                   );

  void fillHltMuonsSlim(const edm::Handle<reco::RecoChargedCandidateCollection> &,
                        const edm::Event   &,
                        int step
                       );

  void fillL2Muons(const edm::Handle<reco::RecoChargedCandidateCollection> &,
                    const edm::Event   &
                   );

  void fillL1Muons(const edm::Handle<l1t::MuonBxCollection> &,
                    const edm::Event   &
                   );

  void fillHltJets(const edm::Handle<reco::PFJetCollection> &,
                        const edm::Event   &
                       );
  void MonteCarloStudies(const edm::Event&);
  

//   virtual void endEvent();

  edm::InputTag offlinePVTag_;
  edm::EDGetTokenT<reco::VertexCollection> offlinePVToken_;
  edm::InputTag offlineMuonTag_;
  edm::EDGetTokenT<std::vector<reco::Muon>> offlineMuonToken_;
  /// file service
  edm::Service<TFileService> outfile_;

  // Trigger process
  edm::InputTag triggerResultTag_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;
  edm::InputTag triggerSummTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummToken_;
  edm::InputTag tagTriggerResultTag_;
  edm::EDGetTokenT<edm::TriggerResults>   tagTriggerResultToken_;
  edm::InputTag tagTriggerSummTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> tagTriggerSummToken_;

  // Input tags
  edm::InputTag l3candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l3candToken_; 
  edm::InputTag iterl3candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> iterl3candToken_; 
  edm::InputTag l2candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l2candToken_; 
  edm::InputTag l1candTag_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; 
  edm::InputTag tkMucandTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> tkMucandToken_; 
  edm::InputTag iterl3OIcandTag_; 
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> iterl3OIcandToken_; 
  edm::InputTag iterl3IOcandTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> iterl3IOcandToken_; 
  edm::InputTag iterl3l1candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> iterl3l1candToken_; 

  edm::InputTag jetTag_;
  edm::EDGetTokenT<reco::PFJetCollection> jetToken_; 

  edm::InputTag chargedDepTag_;
  edm::EDGetTokenT<reco::IsoDepositMap> chargedDepToken_;
  edm::InputTag chargedDep2016Tag_;
  edm::EDGetTokenT<reco::IsoDepositMap> chargedDep2016Token_;

  edm::InputTag neutralDepTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> neutralDepToken_;
  edm::InputTag photonsDepTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> photonsDepToken_;
  edm::InputTag neutralDepTag05_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> neutralDepToken05_;
  edm::InputTag photonsDepTag05_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> photonsDepToken05_;
  edm::InputTag neutralDepTag1_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> neutralDepToken1_;
  edm::InputTag photonsDepTag1_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> photonsDepToken1_;

  edm::InputTag neutralM2DepTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> neutralM2DepToken_;
  edm::InputTag photonsMFDepTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateIsolationMap> photonsMFDepToken_;

  
  edm::InputTag rhoCorrectionTag_;
  edm::EDGetTokenT<double> rhoCorrectionToken_;
  edm::InputTag rhoCorrectionOfflineTag_;
  edm::EDGetTokenT<double> rhoCorrectionOfflineToken_;
 
  edm::InputTag rhoCorrectionOnline05Tag_;
  edm::EDGetTokenT<double> rhoCorrectionOnline05Token_;
  edm::InputTag rhoCorrectionM2Tag_;
  edm::EDGetTokenT<double> rhoCorrectionM2Token_;
  edm::InputTag rhoCorrectionMFTag_;
  edm::EDGetTokenT<double> rhoCorrectionMFToken_;
  edm::InputTag rhoCorrectionM2MFTag_;
  edm::EDGetTokenT<double> rhoCorrectionM2MFToken_;

  edm::InputTag rhoCorrectionECALTag_;
  edm::EDGetTokenT<double> rhoCorrectionECALToken_;
  edm::InputTag rhoCorrectionECALMFTag_;
  edm::EDGetTokenT<double> rhoCorrectionECALMFToken_;
  edm::InputTag rhoCorrectionHCALTag_;
  edm::EDGetTokenT<double> rhoCorrectionHCALToken_;
  edm::InputTag rhoCorrectionHCALM2Tag_;
  edm::EDGetTokenT<double> rhoCorrectionHCALM2Token_;
  edm::InputTag rhoCorrectionECAL05Tag_;
  edm::EDGetTokenT<double> rhoCorrectionECAL05Token_;
  edm::InputTag rhoCorrectionHCAL05Tag_;
  edm::EDGetTokenT<double> rhoCorrectionHCAL05Token_;

  edm::InputTag offlineECalPFTag03_;
  edm::EDGetTokenT<edm::ValueMap<float>> offlineECalPFToken03_;
  edm::InputTag offlineHCalPFTag03_;
  edm::EDGetTokenT<edm::ValueMap<float>> offlineHCalPFToken03_;
  edm::InputTag offlineECalPFTag04_;
  edm::EDGetTokenT<edm::ValueMap<float>> offlineECalPFToken04_;
  edm::InputTag offlineHCalPFTag04_;
  edm::EDGetTokenT<edm::ValueMap<float>> offlineHCalPFToken04_;

  edm::InputTag lumiScalerTag_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalerToken_;

  edm::InputTag puTag_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;

  edm::InputTag genTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

  edm::InputTag beamspotTag_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  bool doOffline_;

  MuonEvent event_;
  std::map<std::string,TTree*> tree_;
  
  unsigned int nGoodVtx; 


};

/// default constructor
MuonNtuples::MuonNtuples(const edm::ParameterSet& cfg): 
  offlinePVTag_           (cfg.getParameter<edm::InputTag>("offlineVtx")), 
    offlinePVToken_         (consumes<reco::VertexCollection>(offlinePVTag_)), 
  offlineMuonTag_         (cfg.getParameter<edm::InputTag>("offlineMuons")), 
    offlineMuonToken_       (consumes<std::vector<reco::Muon>>(offlineMuonTag_)), 

  triggerResultTag_       (cfg.getUntrackedParameter<edm::InputTag>("triggerResult")), 
    triggerResultToken_     (consumes<edm::TriggerResults>(triggerResultTag_)),
  triggerSummTag_         (cfg.getUntrackedParameter<edm::InputTag>("triggerSummary")), 
    triggerSummToken_       (consumes<trigger::TriggerEvent>(triggerSummTag_)),

  tagTriggerResultTag_    (cfg.getUntrackedParameter<edm::InputTag>("tagTriggerResult")), 
    tagTriggerResultToken_  (consumes<edm::TriggerResults>(tagTriggerResultTag_)),
  tagTriggerSummTag_      (cfg.getUntrackedParameter<edm::InputTag>("tagTriggerSummary")), 
    tagTriggerSummToken_    (consumes<trigger::TriggerEvent>(tagTriggerSummTag_)),

  l3candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L3Candidates")),
    l3candToken_            (consumes<reco::RecoChargedCandidateCollection>(l3candTag_)),
  iterl3candTag_          (cfg.getUntrackedParameter<edm::InputTag>("iterL3Candidates")),
    iterl3candToken_        (consumes<reco::RecoChargedCandidateCollection>(iterl3candTag_)),
  l2candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L2Candidates")),
    l2candToken_            (consumes<reco::RecoChargedCandidateCollection>(l2candTag_)),
  l1candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L1Candidates")),
    l1candToken_            (consumes<l1t::MuonBxCollection>(l1candTag_)),
  tkMucandTag_            (cfg.getUntrackedParameter<edm::InputTag>("TkMuCandidates")),
    tkMucandToken_          (consumes<reco::RecoChargedCandidateCollection>(tkMucandTag_)),
  iterl3OIcandTag_        (cfg.getUntrackedParameter<edm::InputTag>("iterL3OICandidates")),
    iterl3OIcandToken_      (consumes<reco::RecoChargedCandidateCollection>(iterl3OIcandTag_)),
  iterl3IOcandTag_        (cfg.getUntrackedParameter<edm::InputTag>("iterL3IOCandidates")),
    iterl3IOcandToken_      (consumes<reco::RecoChargedCandidateCollection>(iterl3IOcandTag_)),
  iterl3l1candTag_        (cfg.getUntrackedParameter<edm::InputTag>("iterL3L1Candidates")),
    iterl3l1candToken_      (consumes<reco::RecoChargedCandidateCollection>(iterl3l1candTag_)),

  jetTag_        (cfg.getUntrackedParameter<edm::InputTag>("JetCandidates")),
    jetToken_      (consumes<reco::PFJetCollection>(jetTag_)),

  chargedDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("ChargedDeposit")), 
    chargedDepToken_        (consumes<reco::IsoDepositMap>(chargedDepTag_)), 
  chargedDep2016Tag_      (cfg.getUntrackedParameter<edm::InputTag>("ChargedDeposit2016")), 
    chargedDep2016Token_    (consumes<reco::IsoDepositMap>(chargedDep2016Tag_)), 

  neutralDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("NeutralDeposit")), 
    neutralDepToken_        (consumes<reco::RecoChargedCandidateIsolationMap>(neutralDepTag_)), 
  photonsDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDeposit")), 
    photonsDepToken_        (consumes<reco::RecoChargedCandidateIsolationMap>(photonsDepTag_)), 
  neutralDepTag05_        (cfg.getUntrackedParameter<edm::InputTag>("NeutralDeposit05")), 
    neutralDepToken05_      (consumes<reco::RecoChargedCandidateIsolationMap>(neutralDepTag05_)), 
  photonsDepTag05_        (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDeposit05")), 
    photonsDepToken05_      (consumes<reco::RecoChargedCandidateIsolationMap>(photonsDepTag05_)), 
  neutralDepTag1_         (cfg.getUntrackedParameter<edm::InputTag>("NeutralDeposit1" )), 
    neutralDepToken1_       (consumes<reco::RecoChargedCandidateIsolationMap>(neutralDepTag1_ )), 
  photonsDepTag1_         (cfg.getUntrackedParameter<edm::InputTag>("PhotonsDeposit1" )), 
    photonsDepToken1_       (consumes<reco::RecoChargedCandidateIsolationMap>(photonsDepTag1_)),

  neutralM2DepTag_          (cfg.getUntrackedParameter<edm::InputTag>("NeutralM2Deposit")), 
    neutralM2DepToken_        (consumes<reco::RecoChargedCandidateIsolationMap>(neutralM2DepTag_)), 
  photonsMFDepTag_          (cfg.getUntrackedParameter<edm::InputTag>("PhotonsMFDeposit")), 
    photonsMFDepToken_        (consumes<reco::RecoChargedCandidateIsolationMap>(photonsMFDepTag_)), 

  rhoCorrectionTag_       (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionOnline")), 
    rhoCorrectionToken_     (consumes<double>(rhoCorrectionTag_)), 
  rhoCorrectionOfflineTag_(cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionOffline")), 
    rhoCorrectionOfflineToken_(consumes<double>(rhoCorrectionOfflineTag_)), 

  rhoCorrectionOnline05Tag_   (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionOnline05")),
    rhoCorrectionOnline05Token_ (consumes<double>(rhoCorrectionOnline05Tag_)),
  rhoCorrectionM2Tag_         (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionM2")),
    rhoCorrectionM2Token_       (consumes<double>(rhoCorrectionM2Tag_)),
  rhoCorrectionMFTag_         (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionMF")),
    rhoCorrectionMFToken_       (consumes<double>(rhoCorrectionMFTag_)),
  rhoCorrectionM2MFTag_       (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionM2MF")),
    rhoCorrectionM2MFToken_     (consumes<double>(rhoCorrectionM2MFTag_)),

  rhoCorrectionECALTag_       (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionECAL")),
    rhoCorrectionECALToken_     (consumes<double>(rhoCorrectionECALTag_)),
  rhoCorrectionECALMFTag_     (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionECALMF")),
    rhoCorrectionECALMFToken_   (consumes<double>(rhoCorrectionECALMFTag_)),
  rhoCorrectionHCALTag_       (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionHCAL")),
    rhoCorrectionHCALToken_     (consumes<double>(rhoCorrectionHCALTag_)),
  rhoCorrectionHCALM2Tag_     (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionHCALM2")),
    rhoCorrectionHCALM2Token_   (consumes<double>(rhoCorrectionHCALM2Tag_)),
  rhoCorrectionECAL05Tag_     (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionECAL05")),
    rhoCorrectionECAL05Token_   (consumes<double>(rhoCorrectionECAL05Tag_)),
  rhoCorrectionHCAL05Tag_     (cfg.getUntrackedParameter<edm::InputTag>("RhoCorrectionHCAL05")),
    rhoCorrectionHCAL05Token_   (consumes<double>(rhoCorrectionHCAL05Tag_)),
  
  offlineECalPFTag03_     (cfg.getUntrackedParameter<edm::InputTag>("offlineECalPFIso03")),
    offlineECalPFToken03_   (consumes<edm::ValueMap<float>>(offlineECalPFTag03_)), 
  offlineHCalPFTag03_     (cfg.getUntrackedParameter<edm::InputTag>("offlineHCalPFIso03")),
    offlineHCalPFToken03_   (consumes<edm::ValueMap<float>>(offlineHCalPFTag03_)), 
  offlineECalPFTag04_     (cfg.getUntrackedParameter<edm::InputTag>("offlineECalPFIso04")),
    offlineECalPFToken04_   (consumes<edm::ValueMap<float>>(offlineECalPFTag04_)), 
  offlineHCalPFTag04_     (cfg.getUntrackedParameter<edm::InputTag>("offlineHCalPFIso04")),
    offlineHCalPFToken04_   (consumes<edm::ValueMap<float>>(offlineHCalPFTag04_)), 

  lumiScalerTag_          (cfg.getUntrackedParameter<edm::InputTag>("lumiScalerTag")),
    lumiScalerToken_        (consumes<LumiScalersCollection>(lumiScalerTag_)), 
  puTag_                  (cfg.getUntrackedParameter<edm::InputTag>("puInfoTag")),
    puToken_                (consumes<std::vector< PileupSummaryInfo>>(puTag_)), 

  genTag_                 (cfg.getUntrackedParameter<edm::InputTag>("genParticlesTag")),
    genToken_               (consumes<reco::GenParticleCollection>(genTag_)), 

  beamspotTag_            (cfg.getParameter<edm::InputTag>("onlineBeamspot")), 
    beamspotToken_         (consumes<reco::BeamSpot>(beamspotTag_)), 

  doOffline_                 (cfg.getUntrackedParameter<bool>("doOffline"))
{
}

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
  if (doOffline_){
    edm::Handle<reco::VertexCollection> vertices; 
    event.getByToken(offlinePVToken_, vertices);
    for(reco::VertexCollection::const_iterator it = vertices->begin(); it != vertices->end(); ++it) {
      if( !it->isValid())  continue;
      nGoodVtx++;
    }
    event_.nVtx = nGoodVtx;
    const reco::Vertex           & pv      = vertices->at(0);


    // Fill offline rho info
    edm::Handle <double>  rhoCollectionOffline;
    event.getByToken(rhoCorrectionOfflineToken_, rhoCollectionOffline);
    if (rhoCollectionOffline.isValid()) event_.rho     = *(rhoCollectionOffline.product());

  // Handle the offline muon collection and fill offline muons
    edm::Handle<std::vector<reco::Muon> > muons;
    event.getByToken(offlineMuonToken_, muons);
    fillMuons(muons, pv, event);

  // Fill bx and inst lumi info
    if (event.isRealData()) {
      event_.bxId  = event.bunchCrossing();

      if (lumiScalerTag_.label() != "none")
      {
        edm::Handle<LumiScalersCollection> lumiScaler;
        event.getByToken(lumiScalerToken_, lumiScaler);

        if (lumiScaler->begin() != lumiScaler->end())
          event_.instLumi = lumiScaler->begin()->instantLumi();
      } 
    }
  }


  // Fill PU info
  if (!event.isRealData()) {
    edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
    if ( event.getByToken(puToken_,puInfo)){
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

  
  // Fill trigger information for probe muon
  edm::Handle<edm::TriggerResults>   triggerResults;
  edm::Handle<trigger::TriggerEvent> triggerEvent;

  if (event.getByToken(triggerResultToken_, triggerResults) &&
      event.getByToken(triggerSummToken_  , triggerEvent)) {
      
    edm::TriggerNames triggerNames_ = event.triggerNames(*triggerResults);
    fillHlt(triggerResults, triggerEvent, triggerNames_, event, false);
  }
  else 
    edm::LogError("") << "Trigger collection for probe muon not found !!!";

  // Fill trigger information for tag muon
  edm::Handle<edm::TriggerResults>   tagTriggerResults;
  edm::Handle<trigger::TriggerEvent> tagTriggerEvent;
      
  if (event.getByToken(tagTriggerResultToken_, tagTriggerResults) &&
      event.getByToken(tagTriggerSummToken_  , tagTriggerEvent)) {
      
    edm::TriggerNames tagTriggerNames_ = event.triggerNames(*tagTriggerResults);
    fillHlt(tagTriggerResults, tagTriggerEvent, tagTriggerNames_, event, true);
  }
  else 
    edm::LogError("") << "Trigger collection for tag muon not found !!!";



 // Handle the online muon collection and fill online muons
  edm::Handle<reco::RecoChargedCandidateCollection> iterl3cands;
  if (event.getByToken(iterl3candToken_, iterl3cands))
    fillHltMuons(iterl3cands, event, false, true, false);//isL3, isIterL3,  isTk,  
  else
    edm::LogWarning("") << "Online iter L3 collection not found !!!";

 // Handle the online muon collection and fill online muons
//   edm::Handle<reco::RecoChargedCandidateCollection> iterl3OIcands;
//   if (event.getByToken(iterl3OIcandToken_, iterl3OIcands))
//     fillHltMuonsSlim(iterl3OIcands, event, 0);  
//   else
//     edm::LogWarning("") << "Iter L3 Step 0 collection not found !!!";
// 
//  // Handle the online muon collection and fill online muons
//   edm::Handle<reco::RecoChargedCandidateCollection> iterl3IOcands;
//   if (event.getByToken(iterl3IOcandToken_, iterl3IOcands))
//     fillHltMuonsSlim(iterl3IOcands, event, 1);
//   else
//     edm::LogWarning("") << "Iter L3 Step 1 collection not found !!!";
//  
//  // Handle the online muon collection and fill online muons
//   edm::Handle<reco::RecoChargedCandidateCollection> iterl3l1cands;
//   if (event.getByToken(iterl3l1candToken_, iterl3l1cands))
//     fillHltMuonsSlim(iterl3l1cands, event, 2);//isL3, isIterL3,  isTk,  
//   else
//     edm::LogWarning("") << "Iter L3 Step 2 collection not found !!!";



 // Handle the online jet collection and fill online muons
  edm::Handle<reco::PFJetCollection> jetcands;
  if (event.getByToken(jetToken_, jetcands))
    fillHltJets(jetcands, event);//isL3, isIterL3,  isTk,  
  else
    edm::LogWarning("") << "IJet collection not found !!!";


 // Handle the online muon collection and fill online muons
  edm::Handle<reco::RecoChargedCandidateCollection> l3cands;
  if (event.getByToken(l3candToken_, l3cands))
    fillHltMuons(l3cands, event,  true, false, false);//isL3, isIterL3,  isTk,  
  else
    edm::LogWarning("") << "Online muon collection not found !!!";

//   // Handle the online tk muon collection and fill online muons
//   edm::Handle<reco::RecoChargedCandidateCollection> tkMucands;
//   if (event.getByToken(tkMucandToken_, tkMucands))
//     fillHltMuons(tkMucands, event, false, false, true);//isL3, isIterL3,  isTk, 
//   else
//     edm::LogWarning("") << "Online tracker muon collection not found !!!";

  // Handle the online muon collection and fill L2 muons
  edm::Handle<reco::RecoChargedCandidateCollection> l2cands;
  if (event.getByToken(l2candToken_, l2cands))
    fillL2Muons(l2cands, event);
  else
    edm::LogWarning("") << "Online L2 muon collection not found !!!";

  // Handle the online muon collection and fill L1 muons
  edm::Handle<l1t::MuonBxCollection> l1cands;
  if (event.getByToken(l1candToken_, l1cands))
    fillL1Muons(l1cands, event);
  else
    edm::LogWarning("") << "Online L1 muon collection not found !!!";
  

  // endEvent();
  tree_["muonTree"] -> Fill();
}



//------------------------------------------------------------------------
void MuonNtuples::MonteCarloStudies(const edm::Event& event)
{
  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(genToken_, genParticles);
  int muId  =    13;

  for ( size_t i=0; i< genParticles->size(); ++i) 
  { 
    const reco::GenParticle &p = (*genParticles)[i];
//     if(fabs(p.pdgId()) < 500 ||  abs(p.pdgId()) > 600  )     continue; 
    
    // only save muons
    if(fabs(p.pdgId()) == muId ){ 
      
      GenParticleCand theGen;
      theGen.pdgId  = p.pdgId();
      theGen.pt     = p.pt() ;
      theGen.eta    = p.eta();
      theGen.phi    = p.phi();
      theGen.energy = p.energy();
      theGen.status = p.status();
    
      unsigned int n_moms = p.numberOfMothers();
      if (n_moms == 0 ){
        theGen.pdgMotherId.push_back(0);
        theGen.pdgRealMother.push_back(0);
      }
      else {
        for (unsigned int im=0; im < n_moms; ++im){
          
          theGen.pdgMotherId .push_back(p.motherRef(im)->pdgId());
          theGen.pdgMotherPt .push_back(p.motherRef(im)->pt());
          theGen.pdgMotherEta.push_back(p.motherRef(im)->eta());
          theGen.pdgMotherPhi.push_back(p.motherRef(im)->phi());
          // if coming from a muon, go back one step ** to be improved **
          if(n_moms == 1 && fabs(p.motherRef(0)->pdgId()) == muId){
            for (unsigned int igm = 0; igm < p.motherRef(0)->numberOfMothers(); igm++){
              theGen.pdgRealMother.push_back(p.motherRef(0)->motherRef(igm)->pdgId());
            }
          }
          else{
            theGen.pdgRealMother.push_back(0);
//             theGen.pdgMotherPt .push_back(-999.);
//             theGen.pdgMotherEta.push_back(-999.);
//             theGen.pdgMotherPhi.push_back(-999.);
          }
        }
      }

    event_.genMuons.push_back(theGen);
    }

    else if(fabs(p.pdgId()) > 500 &&  abs(p.pdgId()) < 600  )     { 

      GenParticleCand theGen;
      theGen.pdgId  = p.pdgId();
      theGen.pt     = p.pt() ;
      theGen.eta    = p.eta();
      theGen.phi    = p.phi();
      theGen.energy = p.energy();
      theGen.status = p.status();

      unsigned int n_moms = p.numberOfMothers();
      if (n_moms == 0 ){
        theGen.pdgMotherId.push_back(0);
        theGen.pdgRealMother.push_back(0);
      }
      else {
        for (unsigned int im=0; im < n_moms; ++im){
          
          theGen.pdgMotherId .push_back(p.motherRef(im)->pdgId());
          theGen.pdgMotherPt .push_back(p.motherRef(im)->pt());
          theGen.pdgMotherEta.push_back(p.motherRef(im)->eta());
          theGen.pdgMotherPhi.push_back(p.motherRef(im)->phi());
        }
      }  

      event_.genBs.push_back(theGen);

    }

  }  // end for genParticles
}



// --------------------------------------------------------------------
void MuonNtuples::fillHlt(const edm::Handle<edm::TriggerResults>   & triggerResults, 
                          const edm::Handle<trigger::TriggerEvent> & triggerEvent  ,
                          const edm::TriggerNames                  & triggerNames  ,
                          const edm::Event                         & event         ,
                          bool                                       isTag         )
{    
   
  for (unsigned int itrig=0; itrig < triggerNames.size(); ++itrig) 
  {
    LogDebug ("triggers") << triggerNames.triggerName(itrig) ;
    if (triggerResults->accept(itrig)) 
    {
      std::string pathName = triggerNames.triggerName(itrig);
      if ( pathName.find ("HLT_IsoMu"  ) !=std::string::npos ||
           pathName.find ("HLT_Mu45"   ) !=std::string::npos ||
           pathName.find ("HLT_Mu"     ) !=std::string::npos ||
           pathName.find ("HLT_TkMu"   ) !=std::string::npos ||
           pathName.find ("HLT_IsoTkMu") !=std::string::npos ||
           pathName.find ("HLT_Mu17"   ) !=std::string::npos ||
           pathName.find ("HLT_Zero"   ) !=std::string::npos ||
           pathName.find ("HLT_L2Mu"   ) !=std::string::npos ||
           pathName.find ("HLT_DoubleL2" ) !=std::string::npos ||
           pathName.find ("HLT_Mu8_"  ) !=std::string::npos
      ){
        if (isTag) event_.hltTag.triggers.push_back(pathName);
        else       event_.hlt   .triggers.push_back(pathName);
      }
    }
  }
     
     
  const trigger::size_type nFilters(triggerEvent->sizeFilters());
  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter)
  {
    std::string filterTag = triggerEvent->filterTag(iFilter).encode();

    if ( ( filterTag.find ("sMu"     ) !=std::string::npos ||
           filterTag.find ("SingleMu") !=std::string::npos ||
           filterTag.find ("SingleMu") !=std::string::npos ||
           filterTag.find ("L2Mu"    ) !=std::string::npos ||
           filterTag.find ("DiMuon"  ) !=std::string::npos 
//            filterTag.find ("DiMuonGlb" ) !=std::string::npos
           ) &&
           filterTag.find ("Tau"       ) ==std::string::npos   &&
           filterTag.find ("EG"        ) ==std::string::npos   &&
           filterTag.find ("MultiFit"  ) ==std::string::npos
       )
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
        
        if (isTag)       event_.hltTag.objects.push_back(hltObj);
        else             event_.hlt   .objects.push_back(hltObj);
        
      }
    }         
  }
  
  // fill hlt rho information
  if (!isTag){
    edm::Handle <double>  hltRhoCollection;
    if (event.getByToken(rhoCorrectionToken_, hltRhoCollection) && hltRhoCollection.isValid()){
       event_.hlt   .rho = *(hltRhoCollection.product());
    } 

    edm::Handle <double>  hltRhoCollection05;
    if (event.getByToken(rhoCorrectionOnline05Token_, hltRhoCollection05) && hltRhoCollection05.isValid()){
      event_.hlt   .rho05 = *(hltRhoCollection05.product());
    }

    edm::Handle <double>  hltRhoCollectionM2;
    if (event.getByToken(rhoCorrectionM2Token_, hltRhoCollectionM2) && hltRhoCollectionM2.isValid()){
       event_.hlt   .rhoM2 = *(hltRhoCollectionM2.product());
    } 
    edm::Handle <double>  hltRhoCollectionMF;
    if (event.getByToken(rhoCorrectionMFToken_, hltRhoCollectionMF) && hltRhoCollectionMF.isValid()){
       event_.hlt   .rhoMF = *(hltRhoCollectionMF.product());
    } 
    edm::Handle <double>  hltRhoCollectionM2MF;
    if (event.getByToken(rhoCorrectionM2MFToken_, hltRhoCollectionM2MF) && hltRhoCollectionM2MF.isValid()){
       event_.hlt   .rhoM2MF = *(hltRhoCollectionM2MF.product());
    } 
  
    // rho ecal
    edm::Handle <double>  hltRhoECALCollection;
    if (event.getByToken(rhoCorrectionECALToken_, hltRhoECALCollection) && hltRhoECALCollection.isValid()){
      event_.hlt   .rho_ecal = *(hltRhoECALCollection.product());
    }
    edm::Handle <double>  hltRhoECALMFCollection;
    if (event.getByToken(rhoCorrectionECALMFToken_, hltRhoECALMFCollection) && hltRhoECALMFCollection.isValid()){
      event_.hlt   .rho_ecalMF = *(hltRhoECALMFCollection.product());
    }
    // rho hcal
    edm::Handle <double>  hltRhoHCALCollection;
    if (event.getByToken(rhoCorrectionHCALToken_, hltRhoHCALCollection) && hltRhoHCALCollection.isValid()){
      event_.hlt   .rho_hcal = *(hltRhoHCALCollection.product());
    }
    edm::Handle <double>  hltRhoHCALM2Collection;
    if (event.getByToken(rhoCorrectionHCALM2Token_, hltRhoHCALM2Collection) && hltRhoHCALM2Collection.isValid()){
      event_.hlt   .rho_hcalM2 = *(hltRhoHCALM2Collection.product());
    }
  
    // rho ecal ( eta up to 5.0 )
    edm::Handle <double>  hltRhoECAL05Collection;
    if (event.getByToken(rhoCorrectionECAL05Token_, hltRhoECAL05Collection) && hltRhoECAL05Collection.isValid()){
      event_.hlt   .rho_ecal05 = *(hltRhoECAL05Collection.product());
    }
    // rho hcal ( eta up to 5.0 )
    edm::Handle <double>  hltRhoHCAL05Collection;
    if (event.getByToken(rhoCorrectionHCAL05Token_, hltRhoHCAL05Collection) && hltRhoHCAL05Collection.isValid()){
      event_.hlt   .rho_hcal05 = *(hltRhoHCAL05Collection.product());
    }
  }
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
    muon::isSoftMuon  ( (*mu1), pv ) ? theMu.isSoft   = 1 : theMu.isSoft   = 0;
    
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
    if (event.getByToken(offlineECalPFToken03_, muECalIsoMap03) &&
        event.getByToken(offlineHCalPFToken03_, muHCalIsoMap03))
    { 
      const edm::ValueMap<float> muECalIso03 = *(muECalIsoMap03);
      const edm::ValueMap<float> muHCalIso03 = *(muHCalIsoMap03);
      theMu.ecalPFCluster_dR03 = muECalIso03[nmuonRef];
      theMu.hcalPFCluster_dR03 = muHCalIso03[nmuonRef];
    }
    else {
//       edm::LogWarning("") << "Offline PF cluster in dR 03 collection not found !!!";
      theMu.ecalPFCluster_dR03 = -9999 ;
      theMu.hcalPFCluster_dR03 = -9999 ;
    }

    edm::Handle< edm::ValueMap<float> > muECalIsoMap04;
    edm::Handle< edm::ValueMap<float> > muHCalIsoMap04;
    if (event.getByToken(offlineECalPFToken04_, muECalIsoMap04) &&
        event.getByToken(offlineHCalPFToken04_, muHCalIsoMap04) 
    ){ 
      const edm::ValueMap<float> muECalIso04 = *(muECalIsoMap04);
      const edm::ValueMap<float> muHCalIso04 = *(muHCalIsoMap04);
      theMu.ecalPFCluster_dR04 = muECalIso04[nmuonRef];
      theMu.hcalPFCluster_dR04 = muHCalIso04[nmuonRef];
    }
    else {
//       edm::LogWarning("") << "Offline PF cluster in dR 04 collection not found !!!";
      theMu.ecalPFCluster_dR04 = -9999 ;
      theMu.hcalPFCluster_dR04 = -9999 ;
    }


    event_.muons.push_back(theMu);
  }
}




// ---------------------------------------------------------------------
void MuonNtuples::fillHltMuons(const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
                               const edm::Event                                        & event   , 
                               bool isL3           ,
                               bool isIterL3       ,
                               bool isTk           
                               )
{

  edm::Handle<reco::IsoDepositMap> trkDepMap;
  edm::Handle<reco::IsoDepositMap> trkDep2016Map;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap;

//   edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap05;
//   edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap05;
//   edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralDepMap1 ;
//   edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsDepMap1 ;

  edm::Handle<reco::RecoChargedCandidateIsolationMap> neutralM2DepMap;
  edm::Handle<reco::RecoChargedCandidateIsolationMap> photonsMFDepMap;

  edm::Handle<reco::BeamSpot>                   theBeamSpot    ;
  event.getByToken(beamspotToken_,              theBeamSpot)   ;
  reco::BeamSpot bs = *theBeamSpot;

  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
    HLTMuonCand theL3Mu;

    reco::RecoChargedCandidateRef candref(l3cands, il3);
    theL3Mu.pt      = candref -> pt();
    theL3Mu.eta     = candref -> eta();
    theL3Mu.phi     = candref -> phi();
    theL3Mu.charge  = candref -> charge();
    
    theL3Mu.dz = candref -> vz();

    reco::TrackRef trkmu = candref->track();
    theL3Mu.trkpt     = trkmu -> pt();
    theL3Mu.dxy       = trkmu -> dxy(bs.position());
    theL3Mu.dxy_error = trkmu -> dxyError();

    if (isIterL3 && event.getByToken(neutralDepToken_, neutralDepMap) ){
      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi = (*neutralDepMap).find( candref );
      theL3Mu.hcalDep    = hcal_mapi->val;
    }
    else    theL3Mu.hcalDep    = -9999;
    
    if (isIterL3 && event.getByToken(photonsDepToken_, photonsDepMap) ){
      reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi = (*photonsDepMap).find( candref );
      theL3Mu.ecalDep    = ecal_mapi->val;
    }
    else    theL3Mu.ecalDep    = -9999;

    if (isIterL3 && event.getByToken(chargedDepToken_, trkDepMap) ){
      reco::IsoDeposit theTkIsolation     = (*trkDepMap)[candref];
      theL3Mu.trkDep    = theTkIsolation.depositWithin(0.3);
    }
    else    theL3Mu.trkDep    = -9999;

//     if (event.getByToken(chargedDep2016Token_, trkDep2016Map) ){
//       reco::IsoDeposit theTkIsolation2016     = (*trkDep2016Map)[candref];
//       theL3Mu.trkDep2016    = theTkIsolation2016.depositWithin(0.3);
//     }
//     else    theL3Mu.trkDep2016    = -9999;
    

    if (isIterL3 && event.getByToken(neutralM2DepToken_, neutralM2DepMap) ){
      reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi = (*neutralM2DepMap).find( candref );
      theL3Mu.hcalM2Dep    = hcal_mapi->val;
    }
    else    theL3Mu.hcalM2Dep    = -9999;

    if (isIterL3 &&  event.getByToken(photonsMFDepToken_, photonsMFDepMap) ){
	  reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi = (*photonsMFDepMap).find( candref );
	  theL3Mu.ecalMFDep    = ecal_mapi->val;
	}
	else    theL3Mu.ecalMFDep    = -9999;


    // fill deposits with veto cones
//     if (event.getByToken(neutralDepToken05_, neutralDepMap05) &&
//         event.getByToken(photonsDepToken05_, photonsDepMap05) &&
//         event.getByToken(neutralDepToken1_,  neutralDepMap1 ) &&
//         event.getByToken(photonsDepToken1_,  photonsDepMap1 ) &&
//         isL3  
//         ){
//         
//       reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi05 = (*neutralDepMap05).find( candref );
//       reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi05 = (*photonsDepMap05).find( candref );
//       reco::RecoChargedCandidateIsolationMap::const_iterator hcal_mapi1  = (*neutralDepMap1 ).find( candref );
//       reco::RecoChargedCandidateIsolationMap::const_iterator ecal_mapi1  = (*photonsDepMap1 ).find( candref );
// 
//       theL3Mu.hcalDep05 = hcal_mapi05->val;
//       theL3Mu.ecalDep05 = ecal_mapi05->val; 
//       theL3Mu.hcalDep1  = hcal_mapi1 ->val;
//       theL3Mu.ecalDep1  = ecal_mapi1 ->val; 
//     }
//     else {
// //       edm::LogWarning("") << "Online PF cluster collection not found !!!";
	theL3Mu.hcalDep05 =  -9999 ;
	theL3Mu.ecalDep05 =  -9999 ; 
	theL3Mu.hcalDep1  =  -9999 ;
	theL3Mu.ecalDep1  =  -9999 ;
//     }


    if       (isL3)      event_.hltmuons   .push_back(theL3Mu);
    else if  (isIterL3)  event_.iterL3muons.push_back(theL3Mu); 
    else if  (isTk)      event_.tkmuons    .push_back(theL3Mu);
  }
}


// ---------------------------------------------------------------------
void MuonNtuples::fillHltMuonsSlim(const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
                                   const edm::Event                                        & event   , 
                                   int step       
                               )
{

  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
    HLTMuonCand theL3Mu;

    reco::RecoChargedCandidateRef candref(l3cands, il3);
    theL3Mu.pt      = candref -> pt();
    theL3Mu.eta     = candref -> eta();
    theL3Mu.phi     = candref -> phi();
    theL3Mu.charge  = candref -> charge();
    
    theL3Mu.dz = candref -> vz();

    reco::TrackRef trkmu = candref->track();
    theL3Mu.trkpt   = trkmu -> pt();

	theL3Mu.ecalDep   =  -9999 ;
	theL3Mu.hcalDep   =  -9999 ;
	theL3Mu.ecalMFDep =  -9999 ;
	theL3Mu.hcalM2Dep =  -9999 ;
	theL3Mu.hcalDep05 =  -9999 ;
	theL3Mu.ecalDep05 =  -9999 ; 
	theL3Mu.hcalDep1  =  -9999 ;
	theL3Mu.ecalDep1  =  -9999 ;
	theL3Mu.trkDep    =  -9999 ;


    if       (step==0)      event_.L3OImuons   .push_back(theL3Mu);
    else if  (step==1)      event_.L3IOmuons   .push_back(theL3Mu); 
    else if  (step==2)      event_.L3L1muons   .push_back(theL3Mu);
  }
}





// ---------------------------------------------------------------------
void MuonNtuples::fillHltJets(const edm::Handle<reco::PFJetCollection> & jets ,
                              const edm::Event          & event)
{

  for (std::vector<reco::PFJet>::const_iterator jets_iter = jets->begin(); jets_iter != jets->end(); ++jets_iter) 
  {
    if (jets_iter -> pt() < 10 ) continue;
    HLTJetCand theJet;

    theJet.pt      = jets_iter -> pt();
    theJet.eta     = jets_iter -> eta();
    theJet.phi     = jets_iter -> phi();
//     theJet.charge  = jets_iter -> charge();
    
//     reco::TrackRef trkmu = candref->track();
//     theL3Mu.trkpt   = trkmu -> pt();
    event_.Jets   .push_back(theJet);
  }
}





// ---------------------------------------------------------------------
void MuonNtuples::fillL1Muons(const edm::Handle<l1t::MuonBxCollection> & l1cands ,
                              const edm::Event                         & event    
                              )
{

  for (int ibx = l1cands->getFirstBX(); ibx <= l1cands->getLastBX(); ++ibx) {
    if (ibx != 0) continue;
    for (auto it = l1cands->begin(ibx); it != l1cands->end(ibx); it++){

      l1t::MuonRef muon(l1cands, distance(l1cands->begin(l1cands->getFirstBX()),it) );

      L1MuonCand theL1Mu;

      theL1Mu.pt       = muon -> pt();
      theL1Mu.eta      = muon -> eta();
      theL1Mu.phi      = muon -> phi();
//       theL1Mu.etaAtVtx = muon -> etaAtVtx();
//       theL1Mu.phiAtVtx = muon -> phiAtVtx();
      theL1Mu.charge   = muon -> charge();
      theL1Mu.quality  = muon -> hwQual();

      event_.L1muons.push_back(theL1Mu);
    }
  }
}

//---------------------------------------------------
void MuonNtuples::fillL2Muons (const edm::Handle<reco::RecoChargedCandidateCollection> & l2cands ,
                               const edm::Event                                        & event   
                               )
{
  for( unsigned int il2 = 0; il2 < l2cands->size(); ++il2) 
  {
    L2MuonCand theL2Mu;

    reco::RecoChargedCandidateRef candref(l2cands, il2);
    theL2Mu.pt      = candref -> pt();
    theL2Mu.eta     = candref -> eta();
    theL2Mu.phi     = candref -> phi();
    theL2Mu.charge  = candref -> charge();
    
    reco::TrackRef mu = candref->get<reco::TrackRef>();
    theL2Mu.minStations = mu->hitPattern().muonStationsWithAnyHits();
    theL2Mu.minHits     = mu->numberOfValidHits();


    event_.L2muons.push_back(theL2Mu);
  }
}

//---------------------------------------------
void MuonNtuples::beginEvent()
{

  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();
  event_.hlt.rho        = -1;
  event_.hlt.rho05      = -1;
  event_.hlt.rhoM2      = -1;
  event_.hlt.rhoMF      = -1;
  event_.hlt.rhoM2MF    = -1;
  event_.hlt.rho_ecal   = -1;
  event_.hlt.rho_hcal   = -1;
  event_.hlt.rho_ecalMF = -1;
  event_.hlt.rho_hcalM2 = -1;
  event_.hlt.rho_ecal05 = -1;
  event_.hlt.rho_hcal05 = -1;  
  
  event_.hltTag.triggers.clear();
  event_.hltTag.objects.clear();
  event_.hltTag.rho = -1;

//   event_.genParticles.clear();
  event_.genMuons.clear();
  event_.genBs.clear();
  event_.muons.clear();
  event_.hltmuons.clear();
  event_.iterL3muons.clear();
  event_.L2muons.clear();
  event_.L1muons.clear();
  event_.tkmuons.clear();

  event_.L3OImuons.clear();
  event_.L3IOmuons.clear();
  event_.L3L1muons.clear();

  event_.Jets.clear();

  for (unsigned int ix=0; ix<3; ++ix) {
    event_.primaryVertex[ix] = 0.;
    for (unsigned int iy=0; iy<3; ++iy) {
      event_.cov_primaryVertex[ix][iy] = 0.;
    }
  }
  event_.nVtx       = -1;
  event_.trueNI     = -1;
  event_.rho        = -1;
  event_.bxId       = -1;
  event_.instLumi   = -1;
  
  nGoodVtx = 0; 
}


// define this as a plug-in
DEFINE_FWK_MODULE(MuonNtuples);
