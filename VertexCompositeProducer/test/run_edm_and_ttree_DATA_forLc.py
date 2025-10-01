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
    fileName =cms.string('TTree_Lc_data.root'))


# Define the input source

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime1/MINIAOD/PromptReco-v2/000/374/681/00000/cd44e980-d445-4a15-9ec1-22b1b7b85253.root'
    ),
        #eventsToProcess = cms.untracked.VEventRange('1:1430:199505260')  # Replace with your specific run, lumi, event numbers
)



#GT and centrality calibration
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')


# Define centrality binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")



# =============== Import Sequences =====================

# Add PbPb collision event selection
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone(
    HLTPaths = [
        "HLT_HIMinimumBias*"
    ]
)

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
    fileName = cms.untracked.string('edm_Lc_data.root'),
        outputCommands = cms.untracked.vstring( #which data to include and exclude 
        "drop *", #no data is kept unless explicitly specified
        'keep *_lamc3pselector_LamC3P_*',  #keep output from fitter.cc
        'keep *_generalLamC3PCandidatesNew_LamC3P_*'  #keep output from sellector.cc
        )
)

process.outputPath = cms.EndPath(process.output)
process.schedule.append(process.outputPath)



