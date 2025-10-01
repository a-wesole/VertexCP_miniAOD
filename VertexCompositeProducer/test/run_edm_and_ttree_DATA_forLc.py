import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023

process = cms.Process('ANASKIM', eras.Run3_2023, Run3_pp_on_PbPb_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')




process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")


# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.TFileService = cms.Service("TFileService",
    fileName =cms.string('TTreeData_Lc_tkpt1_dpt4_maxcand20k_evt100_sept24_v3.root'))


# Define the input source

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    #fileNames = cms.untracked.vstring('file:output_withMC.root'  # Use the EDM output file
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime1/MINIAOD/PromptReco-v2/000/374/681/00000/cd44e980-d445-4a15-9ec1-22b1b7b85253.root'
        #'/store/hidata/HIRun2023A/HIPhysicsRawPrime1/MINIAOD/PromptReco-v2/000/374/668/00000/08a9673f-ec31-4ebf-ab32-7bbd13cbfc1c.root'
    ),
        #eventsToProcess = cms.untracked.VEventRange('1:1430:199505260')  # Replace with your specific run, lumi, event numbers
)



#GT and centrality calibration
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v6', '')
#process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag


# Define centrality binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")



# =============== Import Sequences =====================

# Add PbPb collision event selection
#process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
#process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
#process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')

process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')

#from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
#process.hltobject.triggerNames = trigger_list_data_2023_skimmed

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone(
    HLTPaths = [
        "HLT_HIMinimumBias*"
    ]
)
#process.eventFilter_HM = cms.Sequence(process.hltFilter)


process.event_filters = cms.Sequence(
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter  *
    process.phfCoincFilter2Th4
)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# Define the analysis steps

########## LamC3P candidate rereco ###############################################################


VertexCollection_PAT = "offlineSlimmedPrimaryVertices"
TrackCollection_PAT = "packedPFCandidates"
GenParticleCollection_PAT = "prunedGenParticles"

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalLamC3PCandidates_cff")

process.generalLamC3PCandidatesNew = process.generalLamC3PCandidates.clone()
process.generalLamC3PCandidatesNew.vertexRecoAlgorithm = cms.InputTag(VertexCollection_PAT)
process.generalLamC3PCandidatesNew.trackRecoAlgorithm = cms.InputTag(TrackCollection_PAT)
process.generalLamC3PCandidatesNew.tkEtaDiffCut = cms.double(999.9)
process.generalLamC3PCandidatesNew.tkNhitsCut = cms.int32(11)
process.generalLamC3PCandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalLamC3PCandidatesNew.tkPtCut = cms.double(1.0)#checked
process.generalLamC3PCandidatesNew.alphaCut = cms.double(1.0) #checked
process.generalLamC3PCandidatesNew.alpha2DCut = cms.double(999.9)
process.generalLamC3PCandidatesNew.dPtCut = cms.double(4.0) #checked
#process.generalLamC3PCandidatesNew.lamCMassCut = cms.double(0.125)
process.generalLamC3PCandidatesNew.VtxChiProbCut = cms.double(0.010) #checked


# produce LamC3P trees
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_tree_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3peventinfotree_cff")



process.lamc3pselector = process.lamc3pselector.clone()
process.lamc3pselector.VertexCollection = cms.InputTag(VertexCollection_PAT)
process.lamc3pselector.TrackCollection = cms.InputTag(TrackCollection_PAT)
process.lamc3pselector.GenParticleCollection = cms.InputTag(GenParticleCollection_PAT)
process.lamc3pselector.LamC3P = cms.InputTag("generalLamC3PCandidatesNew:LamC3P")
process.lamc3pselector.DCAValCollection = cms.InputTag("generalLamC3PCandidatesNew:DCAValuesLamC3P")
process.lamc3pselector.DCAErrCollection = cms.InputTag("generalLamC3PCandidatesNew:DCAErrorsLamC3P")
process.lamc3pselector.cand3DDecayLengthSigMin = cms.untracked.double(0.0)
process.lamc3pselector.cand3DPointingAngleMax = cms.untracked.double(1.0)
#process.lamc3pselector.useAnyMVA = cms.bool(False)
#process.lamc3pselector.trkNHitMin = cms.untracked.int32(11)
process.lamc3pselector.trkNHitMin = cms.untracked.int32(-1)


process.LamC3PAna = process.lamc3pana.clone()
process.LamC3PAna.VertexCollection = cms.InputTag(VertexCollection_PAT)
process.LamC3PAna.TrackCollection = cms.InputTag(TrackCollection_PAT)
process.LamC3PAna.GenParticleCollection = cms.InputTag(GenParticleCollection_PAT)
process.LamC3PAna.LamC3P = cms.InputTag("lamc3pselector:LamC3P")
process.LamC3PAna.DCAValCollection = cms.InputTag("lamc3pselector:DCAValuesNewLamC3P")
process.LamC3PAna.DCAErrCollection = cms.InputTag("lamc3pselector:DCAErrorsNewLamC3P")
process.LamC3PAna.isCentrality = cms.bool(True) # Centrality
process.LamC3PAna.centralityBinLabel = cms.InputTag("centralityBin", "HFtowers")#centrality
process.LamC3PAna.centralitySrc = cms.InputTag("hiCentrality") #central
process.LamC3PAna.doGenNtuple = cms.untracked.bool(False) #MConly
process.LamC3PAna.doGenMatching = cms.untracked.bool(False) #MConly
#process.LamC3PAna.useAnyMVA = cms.bool(False); #only set true if you are assigning BDT values +++ change  
#process.LamC3PAna.MVACollection = cms.InputTag("lamc3pselector:MVAValuesNewLamC3P:ANASKIM")


process.LamC3Pana_seq2 = cms.Sequence(process.lamc3pselector * process.LamC3PAna)


process.eventinfoana.selectEvents = cms.untracked.string('EventSelections')
process.eventinfoana.stageL1Trigger = cms.uint32(2)

# Define the process schedule
process.EventSelections = cms.Path(process.centralityBin *
                                   process.hltFilter *
                                   process.event_filters *
                                   process.generalLamC3PCandidatesNew *
                                   process.LamC3Pana_seq2)

process.EventInfoAnalysis = cms.EndPath(process.eventinfoana)
process.schedule = cms.Schedule(process.EventSelections, process.EventInfoAnalysis)




changeToMiniAOD(process)
process.options.numberOfThreads = 1

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('edm.root'),
        outputCommands = cms.untracked.vstring( #which data to include and exclude 
        "drop *", #no data is kept unless explicitly specified
        'keep *_lamc3pselector_LamC3P_*',  # Keep the first collection
        'keep *_generalLamC3PCandidatesNew_LamC3P_*' 
        )
)

process.outputPath = cms.EndPath(process.output)
process.schedule.append(process.outputPath)



