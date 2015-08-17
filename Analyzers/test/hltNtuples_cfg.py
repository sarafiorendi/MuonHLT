import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                      'root://xrootd.unl.edu//store/data/Run2015B/SingleMuon/AOD/PromptReco-v1/000/251/168/00000/02132736-CB26-E511-8127-02163E01386E.root',
                    ),
                    secondaryFileNames = cms.untracked.vstring()
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '74X_dataRun2_HLT_v1'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
# process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
# process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
# process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")


process.muonNtuples =cms.EDAnalyzer("MuonNtuples",
                       offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
                       offlineMuons             = cms.InputTag("muons"),
                       triggerProcess           = cms.string("HLT"),
                       L3Candidates             = cms.untracked.InputTag("hltL3MuonCandidates"), 
					   NeutralDeposit           = cms.untracked.InputTag("hltMuonHcalPFClusterIsoUnseeded"), 
					   PhotonsDeposit           = cms.untracked.InputTag("hltMuonEcalPFClusterIsoUnseeded"), 
                       ChargedDeposit           = cms.untracked.InputTag("hltMuonTkRelIsolationCut0p09Map", "trkIsoDeposits", "HLT"), 
				       RhoCorrectionOnline      = cms.untracked.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"), 
				       RhoCorrectionOffline     = cms.untracked.InputTag("fixedGridRhoFastjetAllCalo"), #to be checked
                       BeamSpotTag              = cms.untracked.InputTag("offlineBeamSpot"),
				       offlineECalPFIso03       = cms.untracked.InputTag("muonEcalPFClusterIsolationProducer03"), 
				       offlineHCalPFIso03       = cms.untracked.InputTag("muonHcalPFClusterIsolationProducer03"), 
				       offlineECalPFIso04       = cms.untracked.InputTag("muonEcalPFClusterIsolationProducer04"), 
				       offlineHCalPFIso04       = cms.untracked.InputTag("muonHcalPFClusterIsolationProducer04"), 
                       )   
process.mypath  = cms.Path(process.muonNtuples)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))   


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
			    	 threshold = cms.untracked.string('ERROR'),
			    	),
#    debugModules  = cms.untracked.vstring('*')
)

