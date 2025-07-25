import FWCore.ParameterSet.Config as cms


d0ana = cms.EDAnalyzer('VertexCompositeTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  doGenNtuple = cms.untracked.bool(False),
  doGenMatching = cms.untracked.bool(False),
  doGenMatchingTOF = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(True),
  isCentrality = cms.untracked.bool(False),
                       threeProngDecay = cms.untracked.bool(False),
                       #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(421),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(321),
  deltaR = cms.untracked.double(0.02),
  VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  TrackCollection = cms.InputTag("packedPFCandidates"),
  VertexCompositeCollection = cms.untracked.InputTag("generalD0CandidatesNew:D0"),
  GenParticleCollection = cms.untracked.InputTag("prunedGenParticles"),
  BSLabel = cms.InputTag("offlineBeamSpot"),

  saveTree = cms.untracked.bool(True),

  useAnyMVA = cms.bool(False),
  isSkimMVA = cms.untracked.bool(False),
  MVACollection = cms.InputTag("generalD0CandidatesNew:MVAValues")
                              )
