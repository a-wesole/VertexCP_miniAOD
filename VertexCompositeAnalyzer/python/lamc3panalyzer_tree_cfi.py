import FWCore.ParameterSet.Config as cms


lamc3pana = cms.EDAnalyzer('VCTreeProducer_LamC3P',
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(True),
  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(4122), #LC
  PID_dau1 = cms.untracked.int32(321),#k-
  PID_dau2 = cms.untracked.int32(2212),#p+
  PID_dau3 = cms.untracked.int32(211),#pi+
  deltaR = cms.untracked.double(0.03),

  VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  TrackCollection = cms.InputTag("packedPFCandidates"),
  #VertexCompositeCollection = cms.InputTag("generalLamC3PCandidatesNew:LamC3P"),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),                           
  centMin = cms.untracked.int32(0),
  centMax = cms.untracked.int32(1000),
                           
  #LamC3P = cms.InputTag("LamC3Pselector:LamC3P"),
  GenParticleCollection = cms.InputTag("prunedGenParticles"),
  BSLabel = cms.InputTag("offlineBeamSpot"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),
  doMuonFull = cms.untracked.bool(False),
  doRecoNtuple = cms.untracked.bool(True),
  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(False),
  saveAllHistogram = cms.untracked.bool(False),
  massHistPeak = cms.untracked.double(1.86),
  massHistWidth = cms.untracked.double(0.2),
  massHistBins = cms.untracked.int32(100),

  pTBins = cms.untracked.vdouble(0,1.2,1.5,2.4,3.0,3.5,4.2,5.0,6.0,7.0,8.0),
  yBins = cms.untracked.vdouble(-2.4,-1.0,0.0,1.0,2.4),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalLamC3PCandidatesNew:MVAValues"),



)


