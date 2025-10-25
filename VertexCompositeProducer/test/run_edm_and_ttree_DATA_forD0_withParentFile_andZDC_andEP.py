import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023


#process = cms.Process('ANASKIM', eras.Run3_2023, Run3_pp_on_PbPb_2023)
process = cms.Process("d0ana", Run3_pp_on_PbPb_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')



process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
#process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")


# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 400
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service("TFileService",
    fileName =cms.string('TTree_fromEDM.root'))


# Define the input source

#miniAOD_files = FileUtils.loadListFromFile('processedFiles_MB11.json')

process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = cms.untracked.vstring(
                                #'root://cms-xrd-global.cern.ch//store/user/nsaha/D0_2023PbPb/HIPhysicsRawPrime11/HIPhysicsRawPrime11_D0_24Aug2025/250825_010454/0000/edm_101.root'
                                '/store/user/nsaha/D0_2023PbPb/HIPhysicsRawPrime15/HIPhysicsRawPrime15_D0_10Sept2025/250910_154532/0000/edm_233.root'
                            ),
                            #secondaryFileNames = fileList.readFiles
                            
                            secondaryFileNames = cms.untracked.vstring(
                                #'/store/hidata/HIRun2023A/HIPhysicsRawPrime11/MINIAOD/PromptReco-v2/000/375/659/00000/a8197d5f-9492-4230-b6c4-f6d0d6249003.root'
                                '/store/hidata/HIRun2023A/HIPhysicsRawPrime15/MINIAOD/PromptReco-v2/000/374/803/00000/3a60ae33-ff73-4889-9a42-6ab394143524.root',
                                '/store/hidata/HIRun2023A/HIPhysicsRawPrime15/MINIAOD/PromptReco-v2/000/374/803/00000/9e7a82d7-b8a6-43f0-8a51-91c9e5c60422.root'
                            ),


                            #firstRun = cms.untracked.uint32(1)
                            #run: 374803 lumi: 522 event: 453201212
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange(
                            #'374803:522'  # run:lumiFirst - run:lumiLast
                            #)
                            
                            
                            #eventsToProcess = cms.untracked.VEventRange('374803:522:453201212')  # Replace with your specific run, lumi, event numbers
                            #eventsToProcess = cms.untracked.VEventRange('374803:522:453813981')  # Replace with your specific run, lumi, event numbers    
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

# ======= ZDC RecHit Producer =========                                                                                                   
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018Producer_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.QWZDC2018RecHit_cfi')
process.load('HeavyIonsAnalysis.ZDCAnalysis.zdcanalyzer_cfi')

process.zdcdigi.SOI = cms.untracked.int32(2)
process.zdcanalyzer.doZDCRecHit = False
process.zdcanalyzer.doZDCDigi = True
process.zdcanalyzer.zdcRecHitSrc = cms.InputTag("QWzdcreco")
process.zdcanalyzer.zdcDigiSrc = cms.InputTag("hcalDigis", "ZDC")
process.zdcanalyzer.calZDCDigi = False
process.zdcanalyzer.verbose = False
process.zdcanalyzer.nZdcTs = cms.int32(6)

process.zdc_seq = cms.Sequence(
    process.zdcdigi +
    process.zdcanalyzer
)

#========== End ZDC =============


# Add PbPb collision event selection
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hffilter_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone(
    HLTPaths = ["HLT_HIMinimumBias*"]
)
process.eventFilter_HLT = cms.Sequence(process.hltFilter)

process.event_filters = cms.Sequence(
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter  *
    process.phfCoincFilter2Th4
)

# Add PbPb event plane
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
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_tree_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff")

process.d0Analyzer = process.d0ana.clone()
process.d0Analyzer.D0 = cms.untracked.InputTag("d0Selector:D0")
process.d0Analyzer.isCentrality = cms.bool(True) # Centrality 
process.d0Analyzer.centralityBinLabel = cms.InputTag("centralityBin", "HFtowers")#centrality
process.d0Analyzer.centralitySrc = cms.InputTag("hiCentrality") #central
process.d0Analyzer.doGenNtuple = cms.untracked.bool(False) #MConly
process.d0Analyzer.doGenMatching = cms.untracked.bool(False) #MConly
process.d0Analyzer.useAnyMVA = cms.bool(True); #only set true if you are assigning BDT values +++ change  
process.d0Analyzer.MVACollection = cms.InputTag("d0Selector:MVAValuesNewD0:ANASKIM")
process.d0Analyzer.MVACollection2 = cms.InputTag("d0Selector:MVAValuesNewD02:ANASKIM")
process.d0Analyzer.ip_tree = cms.bool(True) #True-> Get ip3d,ip3derr from TTree producer!

process.d0ana_seq = cms.Sequence(process.d0Analyzer)

process.eventinfoana = process.eventinfoana.clone()
process.eventinfoana.stageL1Trigger = cms.uint32(2)
process.eventinfoana.VertexCollection = cms.InputTag(VertexCollection_PAT)
process.eventinfoana.TrackCollection = cms.InputTag(TrackCollection_PAT)

process.EventInfoAnalysis = cms.Sequence(process.eventinfoana)

process.Ana_seq = cms.Path(
    process.centralityBin *
    process.eventFilter_HLT *
    process.event_filters *
    process.evtplane_seq *
    process.d0ana_seq *
    process.EventInfoAnalysis *
    process.zdc_seq

)

process.schedule = cms.Schedule(process.Ana_seq)



