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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


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

  edm::InputTag lumiScalerTag_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalerToken_;

  edm::InputTag puTag_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;

  edm::InputTag genTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::InputTag genInfoTag_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

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

  lumiScalerTag_          (cfg.getUntrackedParameter<edm::InputTag>("lumiScalerTag")),
    lumiScalerToken_        (consumes<LumiScalersCollection>(lumiScalerTag_)), 
  puTag_                  (cfg.getUntrackedParameter<edm::InputTag>("puInfoTag")),
    puToken_                (consumes<std::vector< PileupSummaryInfo>>(puTag_)), 

  genTag_                 (cfg.getUntrackedParameter<edm::InputTag>("genParticlesTag")),
    genToken_               (consumes<reco::GenParticleCollection>(genTag_)), 
  genInfoTag_             (cfg.getUntrackedParameter<edm::InputTag>("genInfoTag")),
    genInfoToken_           (consumes<GenEventInfoProduct>(genInfoTag_)), 

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

		  //get the PU pt-hat max
		  float pu_pT_hat_max = -1;
   
		  PileupSummaryInfo puSummary_onTime = puInfo.product()->at(0);
		  for(const auto& pu_pT_hat : puSummary_onTime.getPU_pT_hats()) if (pu_pT_hat>pu_pT_hat_max) pu_pT_hat_max = pu_pT_hat;
          event_.max_pt_hats = pu_pT_hat_max;

          continue;
        }
      }
      
    } 
    else  
      edm::LogError("") << "PU collection not found !!!";
    
    edm::Handle<GenEventInfoProduct> genInfo;
    event.getByToken(genInfoToken_, genInfo);
    event_.qScale = genInfo->qScale();
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
      theGen.isPromptDecayed = p.isPromptDecayed();
    
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
      theGen.isPromptDecayed = p.isPromptDecayed();

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

    if       (isL3)      event_.hltmuons   .push_back(theL3Mu);
    else if  (isIterL3)  event_.iterL3muons.push_back(theL3Mu); 
    else if  (isTk)      event_.tkmuons    .push_back(theL3Mu);
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
      theL1Mu.etaAtVtx = muon -> etaAtVtx();
      theL1Mu.phiAtVtx = muon -> phiAtVtx();
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
  
  event_.hltTag.triggers.clear();
  event_.hltTag.objects.clear();

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
  event_.nVtx        = -1;
  event_.trueNI      = -1;
  event_.bxId        = -1;
  event_.instLumi    = -1;
  
  event_.qScale      = -1;
  event_.max_pt_hats = -1;

  
  nGoodVtx = 0; 
}


// define this as a plug-in
DEFINE_FWK_MODULE(MuonNtuples);
