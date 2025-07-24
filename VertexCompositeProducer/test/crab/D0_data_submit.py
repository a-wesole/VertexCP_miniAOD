from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.workArea = 'D0_2023PbPb'
config.General.requestName = 'HIPhysicsRawPrime0_D0_16July2025'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../run_edm_and_ttree_DATA.py'
config.JobType.numCores = 2
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD' 
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5 
config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Golden.json'
config.Data.outLFNDirBase = '/store/user/awesolek/D0_2023PbPb' 
config.Data.publication = True
config.Data.outputDatasetTag = 'HIPhysicsRawPrime0_D0_16July2025'

config.Site.storageSite = 'T2_US_Purdue'
#config.Site.ignoreGlobalBlacklist = True
