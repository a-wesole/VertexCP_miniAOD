from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.workArea = 'D0_2023PbPb_ReProduced_Trees'
config.General.requestName = 'HIPhysicsRawPrime21_Oct25'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../run_edm_and_ttree_DATA_forD0_withParentFile_andZDC_andEP.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['../HeavyIonRPRcd_PbPb2023_offline.db']
config.Data.inputDataset = ''
config.Data.splitting = 'FileBased'
config.Data.useParent = True
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1

config.Data.outLFNDirBase = '/store/user/wxie/D0_2023PbPb_ReProduced_Trees/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'HIPhysicsRawPrime21_Oct25'

config.Site.storageSite = 'T2_US_Purdue'

