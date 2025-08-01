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
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service("TFileService",
    fileName =cms.string('TTree.root'))


# Define the input source

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    #fileNames = cms.untracked.vstring('file:output_withMC.root'  # Use the EDM output file
    fileNames = cms.untracked.vstring(
       #'root://xrootd-cms.infn.it//store/mc/HINPbPbSpring23MiniAOD/promptD0ToKPi_PT-1_TuneCP5_5p36TeV_pythia8-evtgen/MINIAODSIM/132X_mcRun3_2023_realistic_HI_v9-v2/2560000/04335bea-a283-40ea-a050-d71e1b7fac6b.root'
       'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/951/00000/717c1c52-b5a8-4864-89c1-bde6e4184087.root'
       #'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/681/00000/89587dae-c774-4489-a9ea-130849a72872.root'
       #'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/668/00000/06179488-b7e6-44f6-bec9-eb242a290ffd.root'

    #'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/375/055/00000/2d8cd07d-f92f-44df-8e0f-eb28dca3108b.root'
        #'file:2d8cd07d-f92f-44df-8e0f-eb28dca3108b.root'


    ),
        #lumisToProcess = cms.untracked.VLuminosityBlockRange(
        #'374951:30-374951:30'  # run:lumiFirst - run:lumiLast
        #)

            #eventsToProcess = cms.untracked.VEventRange('375055:201:128206575')  # Replace with your specific run, lumi, event numbers    
)



#GT and centrality calibration
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag


# Define centrality binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")


# =============== Import Sequences =====================

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
# process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone(
    HLTPaths = [
        "HLT_HIMinimumBias*"
    ]
)
process.eventFilter_HLT = cms.Sequence(process.hltFilter)

process.event_filters = cms.Sequence(
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter  *
    process.phfCoincFilter2Th4
)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD


VertexCollection_PAT = "offlineSlimmedPrimaryVertices"
TrackCollection_PAT = "packedPFCandidates"
GenParticleCollection_PAT = "prunedGenParticles"
TrkChi2Label = "packedPFCandidateTrackChi2"


# Define the analysis steps

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.trackRecoAlgorithm = cms.InputTag(TrackCollection_PAT)
process.generalD0CandidatesNew.vertexRecoAlgorithm = cms.InputTag(VertexCollection_PAT)
process.generalD0CandidatesNew.GenParticleCollection = cms.untracked.InputTag(GenParticleCollection_PAT)
process.generalD0CandidatesNew.TrackChi2Label = cms.InputTag(TrkChi2Label)
process.generalD0CandidatesNew.tkEtaDiffCut = cms.double(999.9)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(1.0)
process.generalD0CandidatesNew.alphaCut = cms.double(0.30)
process.generalD0CandidatesNew.alpha2DCut = cms.double(999.9)
process.generalD0CandidatesNew.dPtCut = cms.double(0.0)
process.generalD0CandidatesNew.mPiKCutMin = cms.double(1.74)
process.generalD0CandidatesNew.mPiKCutMax = cms.double(2.00)
process.generalD0CandidatesNew.d0MassCut = cms.double(0.125)
process.generalD0CandidatesNew.VtxChiProbCut = cms.double(0.010)


# produce D0 trees
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff")


process.d0Selector = process.d0selector.clone()
process.d0Selector.trackRecoAlgorithm = cms.InputTag(TrackCollection_PAT)
process.d0Selector.vertexRecoAlgorithm = cms.InputTag(VertexCollection_PAT)
process.d0Selector.GenParticleCollection = cms.untracked.InputTag(GenParticleCollection_PAT)
process.d0Selector.D0 = cms.InputTag("generalD0CandidatesNew:D0")
process.d0Selector.DCAValCollection = cms.InputTag("generalD0CandidatesNew:DCAValuesD0")
process.d0Selector.DCAErrCollection = cms.InputTag("generalD0CandidatesNew:DCAErrorsD0")
process.d0Selector.D0 = cms.InputTag("generalD0CandidatesNew:D0")
process.d0Selector.cand3DDecayLengthSigMin = cms.untracked.double(0.)
process.d0Selector.cand3DPointingAngleMax = cms.untracked.double(1.0)
process.d0Selector.input_names = cms.vstring('input')
process.d0Selector.output_names = cms.vstring('probabilities')
process.d0Selector.onnxModelFileName = cms.string("XGBoost_Model_0428_0_OnlyPrompt.onnx")
process.d0Selector.bdtCutsFile = cms.string("bdt_cuts.csv")
process.d0Selector.mvaCut = cms.double(0.9)
process.d0Selector.trkNHitMin = cms.untracked.int32(11)
process.d0Selector.isCentrality = cms.bool(True) # Centrality 
process.d0Selector.useAnyMVA = cms.bool(True)#only set true if you are assigning BDT values  +++change 
process.d0Selector.applyXGB = cms.bool(True)#only set true if you are assigning XGB values


process.d0Analyzer = process.d0ana.clone()
process.d0Analyzer.trackRecoAlgorithm = cms.InputTag(TrackCollection_PAT)
process.d0Analyzer.vertexRecoAlgorithm = cms.InputTag(VertexCollection_PAT)
process.d0Analyzer.GenParticleCollection = cms.untracked.InputTag(GenParticleCollection_PAT)
process.d0Analyzer.D0 = cms.untracked.InputTag("d0Selector:D0")
process.d0Analyzer.DCAValCollection = cms.InputTag("d0Selector:DCAValuesNewD0")
process.d0Analyzer.DCAErrCollection = cms.InputTag("d0Selector:DCAErrorsNewD0")
process.d0Analyzer.isCentrality = cms.bool(True) # Centrality 
process.d0Analyzer.centralityBinLabel = cms.InputTag("centralityBin", "HFtowers")#centrality
process.d0Analyzer.centralitySrc = cms.InputTag("hiCentrality") #central
process.d0Analyzer.doGenNtuple = cms.untracked.bool(False) #MConly
process.d0Analyzer.doGenMatching = cms.untracked.bool(False) #MConly
process.d0Analyzer.useAnyMVA = cms.bool(True); #only set true if you are assigning BDT values +++ change  
process.d0Analyzer.MVACollection = cms.InputTag("d0Selector:MVAValuesNewD0:ANASKIM")
process.d0Analyzer.MVACollection2 = cms.InputTag("d0Selector:MVAValuesNewD02:ANASKIM")


process.d0ana_seq2 = cms.Sequence(process.d0Selector * process.d0Analyzer)

#eventinfoana must be in EndPath, and process.eventinfoana.selectEvents must be the name of eventFilter_HM Path
process.eventinfoana.selectEvents = cms.untracked.string('EventSelections')
process.eventinfoana.stageL1Trigger = cms.uint32(2)

process.EventSelections = cms.Path(
    process.centralityBin *
    process.eventFilter_HLT *
    process.event_filters * 
    process.generalD0CandidatesNew *
    process.d0ana_seq2 
    #process.generalD0CandidatesNew 
)

process.EventInfoAnalysis = cms.EndPath(process.eventinfoana)
process.schedule = cms.Schedule(process.EventSelections, process.EventInfoAnalysis)

changeToMiniAOD(process)
process.options.numberOfThreads = 2

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('edm.root'),
        outputCommands = cms.untracked.vstring( #which data to include and exclude 
        "drop *", #no data is kept unless explicitly specified
        #'keep *_generalD0CandidatesNew_D0_*', 
        'keep *_d0Selector_*_*',  # Keep the MVA collection (adjust the label)

        )
)

process.outputPath = cms.EndPath(process.output)
process.schedule.append(process.outputPath)



