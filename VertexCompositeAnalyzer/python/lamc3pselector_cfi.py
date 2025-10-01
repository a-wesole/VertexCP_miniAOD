#Author: Nihar Saha

# This file contains the default values of parameters passed to VertexCompositeAnalyzer/plugins/VCTreeProducer_LamC3P.cc
#if no other values are specfied in the config file, these will be used.
import FWCore.ParameterSet.Config as cms

lamc3pselector = cms.EDProducer('VCSelector_LamC3P',
  doGenMatching = cms.bool(False),
  hasSwap = cms.untracked.bool(False),
  decayInGen = cms.untracked.bool(False),
  twoLayerDecay = cms.untracked.bool(False),
  threeProngDecay = cms.untracked.bool(True),

  selectGenMatch = cms.untracked.bool(False),
  selectGenUnMatch = cms.untracked.bool(False),
  selectGenMatchSwap = cms.untracked.bool(False),
  selectGenMatchUnSwap = cms.untracked.bool(False),

  #PID used only for GEN and/or GEN match
  PID = cms.untracked.int32(4122),
  PID_dau1 = cms.untracked.int32(321),
  PID_dau2 = cms.untracked.int32(2212),
  PID_dau3 = cms.untracked.int32(211),
  VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  TrackCollection = cms.InputTag("packedPFCandidates"),
  VertexCompositeCollection = cms.InputTag("generalLamC3PCandidatesNew:LamC3P"),
  GenParticleCollection = cms.InputTag("prunedGenParticles"),
  MuonCollection = cms.untracked.InputTag("null"),  
  doMuon = cms.untracked.bool(False),
  doMuonFull = cms.untracked.bool(False),
  
  useAnyMVA = cms.bool(False),
  useExistingMVA = cms.bool(False),
  mvaType = cms.string('BDT'),
  GBRForestLabel = cms.string('LamC3PInpPb'),
  GBRForestFileName = cms.string('GBRForestfile_BDT_LamC3PInpPb_1_2.root'),
  MVACollection = cms.InputTag("generalLamC3PCandidatesNew:MVAValues"),
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

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),
  centMin = cms.untracked.int32(0),
  centMax = cms.untracked.int32(1000)
)

