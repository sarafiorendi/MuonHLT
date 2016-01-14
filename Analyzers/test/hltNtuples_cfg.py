import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/04683826-796B-E511-B62E-02163E0118CD.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/04A3BD6D-996B-E511-855D-02163E011D7F.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/06E0F158-7C6B-E511-93AD-02163E014782.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/0C80997A-996B-E511-A8FC-02163E011BB1.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/12210785-756B-E511-808B-02163E013825.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/14330CC0-7B6B-E511-96AD-02163E0142B7.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/16F9EC55-7E6B-E511-B636-02163E012447.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/200EBEA5-9D6B-E511-8DCB-02163E014300.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/20C97268-7A6B-E511-943F-02163E01342D.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/2C02726D-816B-E511-8D72-02163E0145D0.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/2CA6AC11-816B-E511-8A0F-02163E014331.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/2CABF288-7D6B-E511-9B1D-02163E0118BA.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/309203DC-AB6B-E511-AA13-02163E0144D0.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/32C0396D-B46B-E511-88FD-02163E014281.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/32F1C94E-716B-E511-8F6D-02163E0133CD.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/3892B4BB-D16B-E511-8A52-02163E01459B.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/3E7D2ABF-626B-E511-BA69-02163E0127A0.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/4466563B-676B-E511-9A22-02163E0134EE.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/48E41DB2-196C-E511-8591-02163E0143C6.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/48F582B3-846B-E511-90A4-02163E0134BE.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/4ABEAB5C-836B-E511-8991-02163E012013.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/4C218CCD-636B-E511-A538-02163E014387.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/501E5346-826B-E511-ABEC-02163E0137C9.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/50CE7B5C-7F6B-E511-971D-02163E0144B7.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/56BD431C-766B-E511-A5FB-02163E0119F6.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/580A378D-8B6B-E511-A730-02163E011A55.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/5E098905-806B-E511-8302-02163E01195C.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/6069CB35-926B-E511-8476-02163E012466.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/628EC904-626B-E511-B548-02163E013424.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/64D62970-656B-E511-AB98-02163E012AA9.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/6A3173F9-776B-E511-9002-02163E0143D5.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/6A78E02A-836B-E511-B4D0-02163E01195C.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/6E000206-826B-E511-BB2E-02163E013414.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/7467041B-806B-E511-A3FB-02163E012A29.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/76AE2C4B-6B6B-E511-B017-02163E01245B.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/863282FF-776B-E511-B26E-02163E0140D8.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/869F4396-7B6B-E511-BF8B-02163E0119CF.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/8AB52480-646B-E511-9A1B-02163E0139DA.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/8ACB3648-7E6B-E511-BB3F-02163E014587.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/AE57DC41-896B-E511-B707-02163E013557.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/AECD3881-956B-E511-A9DA-02163E011BD5.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/B0094392-046C-E511-AECC-02163E0142B7.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/B2BC55A4-7B6B-E511-A618-02163E0142F4.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/BAE3552B-846B-E511-82CA-02163E0143E1.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/BE34F645-7E6B-E511-9BA0-02163E012A6D.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/C255150B-666B-E511-8586-02163E0143EA.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/CC47BF4B-7A6B-E511-B4D6-02163E01376C.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/CE2DF30A-776B-E511-B1F8-02163E014387.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/CE7776F2-A36B-E511-9FF2-02163E011F30.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/D6718406-BC6B-E511-84CC-02163E0136EA.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/DE7DB122-8C6B-E511-B8B2-02163E0137BC.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/E062CB4C-7F6B-E511-85F8-02163E01428E.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/E2470B55-876B-E511-86FF-02163E013922.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/E63219BC-7D6B-E511-8F67-02163E0133CC.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/E86B2A8E-7D6B-E511-B9A4-02163E0137AE.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/EC86A4D2-686B-E511-AEF6-02163E0138FC.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/F058B8EC-736B-E511-8902-02163E01338D.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/F4DC7DAC-666B-E511-8122-02163E013507.root',
#                       '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/F6220E3C-866B-E511-BC38-02163E01379F.root',
                      '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/FA0B284A-856B-E511-9378-02163E0133E8.root',
                      '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/FAC77F63-8F6B-E511-9911-02163E0145FD.root',
                      '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/FE6D8959-7C6B-E511-B41B-02163E0137BC.root',
                    ),
                    secondaryFileNames = cms.untracked.vstring(),
#                     lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
# process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
# process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
# process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")


process.muonNtuples =cms.EDAnalyzer("MuonNtuples",
                       offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
                       offlineMuons             = cms.InputTag("muons"),
                       triggerProcess           = cms.string("HLT"),
                       tagTriggerProcess        = cms.string("HLT"),
                       L3Candidates             = cms.untracked.InputTag("hltL3MuonCandidates"), 
					   NeutralDeposit           = cms.untracked.InputTag("hltMuonHcalPFClusterIsoForMuons"), 
					   PhotonsDeposit           = cms.untracked.InputTag("hltMuonEcalPFClusterIsoForMuons"), 
					   NeutralDeposit05         = cms.untracked.InputTag("hltMuonHcalPFClusterIsoForMuonsNoEffAreaVeto0p05"), 
					   PhotonsDeposit05         = cms.untracked.InputTag("hltMuonEcalPFClusterIsoForMuonsNoEffAreaVeto0p05"), 
					   NeutralDeposit1          = cms.untracked.InputTag("hltMuonHcalPFClusterIsoForMuonsNoEffAreaVeto0p1"), 
					   PhotonsDeposit1          = cms.untracked.InputTag("hltMuonEcalPFClusterIsoForMuonsNoEffAreaVeto0p1"), 
                       ChargedDeposit           = cms.untracked.InputTag("hltMuonTkRelIsolationCut0p09Map", "trkIsoDeposits", "HLT"), 
				       RhoCorrectionOnline      = cms.untracked.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"), 
				       RhoCorrectionOffline     = cms.untracked.InputTag("fixedGridRhoFastjetAllCalo"), #to be checked
                       BeamSpotTag              = cms.untracked.InputTag("offlineBeamSpot"),
				       offlineECalPFIso03       = cms.untracked.InputTag("muonEcalPFClusterIsolationProducer03"), 
				       offlineHCalPFIso03       = cms.untracked.InputTag("muonHcalPFClusterIsolationProducer03"), 
				       offlineECalPFIso04       = cms.untracked.InputTag("muonEcalPFClusterIsolationProducer04"), 
				       offlineHCalPFIso04       = cms.untracked.InputTag("muonHcalPFClusterIsolationProducer04"), 
                       lumiScalerTag            = cms.InputTag("scalersRawToDigi")
                       )   
process.mypath  = cms.Path(process.muonNtuples)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonNtuple_run258158_test.root"),
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

