import FWCore.ParameterSet.Config as cms

d0selector = cms.EDProducer('VertexCompositeSelector',
  doGenMatching = cms.untracked.bool(False),
  hasSwap = cms.untracked.bool(True),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(False),

  selectGenMatch = cms.untracked.bool(False),
  selectGenUnMatch = cms.untracked.bool(False),
  selectGenMatchSwap = cms.untracked.bool(False),
  selectGenMatchUnSwap = cms.untracked.bool(False),

  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(421),
  PID_dau1 = cms.untracked.int32(211),
  PID_dau2 = cms.untracked.int32(321),
  VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  TrackCollection = cms.InputTag("packedPFCandidates"),
  VertexCompositeCollection = cms.InputTag("generalD0CandidatesNew:D0"),
  GenParticleCollection = cms.untracked.InputTag("prunedGenParticles"),
  MuonCollection = cms.untracked.InputTag("null"),
  doMuon = cms.untracked.bool(False),
  doMuonFull = cms.untracked.bool(False),
  
  useAnyMVA = cms.bool(False),
  useExistingMVA = cms.bool(False),
  mvaType = cms.string('BDT'),
  GBRForestLabel = cms.string('D0InpPb'),
  GBRForestFileName = cms.string('GBRForestfile_BDT_D0InpPb_1_2.root'),
  MVACollection = cms.InputTag("generalD0CandidatesNew:MVAValues"),
  mvaMax = cms.untracked.double(999.9),
  mvaMin = cms.untracked.double(-999.9),
  mvaCuts = cms.vdouble(-1.,0,0,0,0),

  trkHighPurity = cms.untracked.bool(True),
  trkPMin = cms.untracked.double(0.),
  trkPtMin = cms.untracked.double(0.),
  trkEtaMax = cms.untracked.double(999.),
  trkPSumMin = cms.untracked.double(0.),
  trkPtSumMin = cms.untracked.double(0.),
  trkPtAsymMin = cms.untracked.double(0.),
  trkEtaDiffMax = cms.untracked.double(999.),
  trkPtErrMax = cms.untracked.double(999.),
  trkNHitMin = cms.untracked.int32(0),
  candpTMin = cms.untracked.double(-999.),
  candpTMax = cms.untracked.double(999.),
  candYMin = cms.untracked.double(-999.),
  candYMax = cms.untracked.double(999.),
  cand3DDecayLengthSigMin = cms.untracked.double(0.),
  cand3DPointingAngleMax = cms.untracked.double(999.),
  cand2DDecayLengthSigMin = cms.untracked.double(0.),
  cand2DPointingAngleMax = cms.untracked.double(999.),
  cand3DDCAMin = cms.untracked.double(-999.),
  cand3DDCAMax = cms.untracked.double(999.),
  cand2DDCAMin = cms.untracked.double(-999.),
  cand2DDCAMax = cms.untracked.double(999.),
  candVtxProbMin = cms.untracked.double(0.),

  isCentrality = cms.bool(False),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality")
                              )

