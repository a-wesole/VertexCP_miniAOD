from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Run3_D0_Data'
config.General.workArea = 'Run3_D0_Data'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../run_edm_and_ttree_DATA.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 4500

config.Data.inputDataset = '/HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD' 
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 2 
config.Data.lumiMask = 'mini2.json'
#for part of dataset0
config.Data.outLFNDirBase = '/store/user/awesolek/Run3_D0_Data' 
config.Data.publication = True
config.Data.outputDatasetTag = 'Run3_D0_Data'

config.Site.storageSite = 'T2_US_Purdue'
#config.Site.ignoreGlobalBlacklist = True
