from CRABClient.UserUtilities import config #, getUsernameFromSiteDB
config = config()

config.General.requestName = 'promptD0ToKPi_PT-10'
config.General.workArea = 'Run3_D0_MC'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/awesole/VertexCP_clean/CMSSW_13_2_11/src/VertexCompositeAnalysis/VertexCompositeProducer/test/run_edm_and_ttree_MC.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/promptD0ToKPi_PT-10_TuneCP5_5p36TeV_pythia8-evtgen/HINPbPbSpring23MiniAOD-132X_mcRun3_2023_realistic_HI_v9-v1/MINIAODSIM' 
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1 
config.Data.lumiMask = ''
config.Data.outLFNDirBase = '/store/user/wxie/Run3_D0_MC' 
config.Data.publication = False
config.Data.outputDatasetTag = 'Run3_D0_MC'

config.Site.storageSite = 'T2_US_Purdue'
#config.Site.ignoreGlobalBlacklist = True
