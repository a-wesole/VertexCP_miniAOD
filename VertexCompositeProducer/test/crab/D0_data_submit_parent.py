from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.workArea = 'D0_2023PbPb_new'
config.General.requestName = 'HIPhysicsRawPrime15_Oct17'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../run_edm_and_ttree_DATA_forD0_withParentFile_andZDC_andEP.py'
#config.JobType.numCores = 2
#config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['../HeavyIonRPRcd_PbPb2023_offline.db']
config.Data.inputDataset = '/HIPhysicsRawPrime15/nsaha-HIPhysicsRawPrime15_D0_10Sept2025-8b581d8886e57c0efd46b0be657273da/USER'
#config.Data.inputDataset = '/HIPhysicsRawPrime8/nsaha-HIPhysicsRawPrime8_D0_7Sept2025-8b581d8886e57c0efd46b0be657273da/USER'
#config.Data.secondaryInputDataset = '/HIPhysicsRawPrime8/HIRun2023A-PromptReco-v2/MINIAOD'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.useParent = True
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1

config.Data.outLFNDirBase = '/store/user/nsaha/D0_2023PbPb_new/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'HIPhysicsRawPrime15_Oct17'

config.Site.storageSite = 'T2_US_Purdue'
#config.Site.ignoreGlobalBlacklist = True
