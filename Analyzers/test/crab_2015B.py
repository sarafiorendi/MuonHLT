from CRABClient.UserUtilities import config
config = config()

config.General.requestName  = 'muonNtuples_2015B_MuonJsonv4'
config.General.workArea     = '13TeV/'

config.JobType.pluginName   = 'Analysis'
config.JobType.psetName     = 'hltNtuples_cfg.py'
config.JobType.outputFiles  = ['muonNtuple.root']

config.Data.inputDataset    = '/SingleMuon/Run2015B-PromptReco-v1/AOD'
config.Data.inputDBS        = 'global'
config.Data.lumiMask        = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_MuonPhys_v4.txt'
config.Data.splitting       = 'LumiBased'
config.Data.unitsPerJob     = 40
config.Data.outLFNDirBase   = '/store/group/phys_muon/fiorendi/13TeV/2015B/' # or '/store/group/<subdir>'

config.Site.storageSite     = 'T2_CH_CERN'
