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
    fileName =cms.string('TTree_D0_data_Oct22.root'))


# Define the input source

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
       '/store/hidata/HIRun2023A/HIPhysicsRawPrime15/MINIAOD/PromptReco-v2/000/374/803/00000/3a60ae33-ff73-4889-9a42-6ab394143524.root'
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
process.eventFilter_HLT = cms.Sequence(process.hltFilter)

process.event_filters = cms.Sequence(
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter  *
    process.phfCoincFilter2Th4
)

process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")

process.hiEvtPlane.trackTag = cms.InputTag("packedPFCandidates")
process.hiEvtPlane.vertexTag = cms.InputTag("offlineSlimmedPrimaryVertices")
process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlineSlimmedPrimaryVertices")

process.hiEvtPlane.loadDB = cms.bool(True)
process.hiEvtPlaneFlat.centralityVariable=process.hiEvtPlane.centralityVariable
process.hiEvtPlaneFlat.vertexTag=process.hiEvtPlane.vertexTag
process.hiEvtPlaneFlat.flatminvtx=process.hiEvtPlane.flatminvtx
process.hiEvtPlaneFlat.flatnvtxbins=process.hiEvtPlane.flatnvtxbins
process.hiEvtPlaneFlat.flatdelvtx=process.hiEvtPlane.flatdelvtx
process.hiEvtPlaneFlat.FlatOrder=process.hiEvtPlane.FlatOrder
process.hiEvtPlaneFlat.CentBinCompression=process.hiEvtPlane.CentBinCompression
process.hiEvtPlaneFlat.caloCentRef=process.hiEvtPlane.caloCentRef
process.hiEvtPlaneFlat.caloCentRefWidth=process.hiEvtPlane.caloCentRefWidth

process.CondDB.connect = "sqlite_file:HeavyIonRPRcd_PbPb2023_offline.db"
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                      process.CondDB,
                                      toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                 tag = cms.string('HeavyIonRPRcd')
                                      )
                                    )
)
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')
process.evtplane_seq = cms.Sequence(process.hiEvtPlane * process.hiEvtPlaneFlat)



#from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD


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
process.d0Analyzer.ip_tree = cms.bool(False)

process.d0ana_seq2 = cms.Sequence(process.d0Selector * process.d0Analyzer)


process.eventinfoana = process.eventinfoana.clone()
process.eventinfoana.stageL1Trigger = cms.uint32(2)
process.eventinfoana.VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")

process.EventInfoAnalysis = cms.Sequence(process.eventinfoana)

process.Ana_seq = cms.Path(
    process.centralityBin *
    process.eventFilter_HLT *
    process.event_filters * 
    process.evtplane_seq *
    process.generalD0CandidatesNew *
    process.d0ana_seq2 *
    process.EventInfoAnalysis
)

process.schedule = cms.Schedule(process.Ana_seq)


#changeToMiniAOD(process)
process.options.numberOfThreads = 2

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('edm_D0_data.root'),
        outputCommands = cms.untracked.vstring( #which data to include and exclude 
        "drop *", #no data is kept unless explicitly specified
        #'keep *_generalD0CandidatesNew_D0_*', # Keep everything from D0Fitter.cc 
        'keep *_d0Selector_*_*',  # Keep everything from selector.cc

        )
)

process.outputPath = cms.EndPath(process.output)
process.schedule.append(process.outputPath)



