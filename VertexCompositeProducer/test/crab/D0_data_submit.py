from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = 'test_1_file'
config.General.workArea = 'test_1_file'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/awesole/VertexCP_2/CMSSW_13_2_11/src/VertexCompositeAnalysis/VertexCompositeProducer/test/run_edm_and_ttree_DATA.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 4500

config.Data.inputDataset = '/HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD' 
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1 
config.Data.lumiMask = 'mini.json'
#for part of dataset0
config.Data.outLFNDirBase = '/store/user/awesolek/test_1_file' 
#config.Data.publication = True
config.Data.publication = False
config.Data.outputDatasetTag = 'test_1_file'

config.Site.storageSite = 'T2_US_Purdue'
#config.Site.ignoreGlobalBlacklist = True
