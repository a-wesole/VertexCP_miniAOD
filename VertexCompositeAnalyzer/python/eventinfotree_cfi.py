import FWCore.ParameterSet.Config as cms

eventinfoana = cms.EDAnalyzer('EventInfoTreeProducer',
  beamSpotSrc = cms.InputTag("offlineBeamSpot"),
  VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
  TrackCollection = cms.InputTag("packedPFcandidates"),

  #Trigger info
  #Trigger info
  #TriggerResultCollection = cms.untracked.InputTag("TriggerResults::HLT"),
  #triggerPathNames = cms.untracked.vstring(
   #   'HLT_HIZeroBias',
      # Minimum Bias trigger
    #  'HLT_HIMinimumBias',
                              #),
  #Filter info
  #FilterResultCollection = cms.untracked.InputTag("TriggerResults::ANASKIM"),
  #eventFilterNames = cms.untracked.vstring(
  #    'Flag_colEvtSel',
  #    'Flag_hfCoincFilter2Th4',
  #    'Flag_primaryVertexFilter',
  #    'Flag_clusterCompatibilityFilter',
                              #),
  selectEvents = cms.untracked.string(""),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),

  #threeProngDecay = cms.untracked.bool(True),
                              
  isEventPlane = cms.bool(True),
  eventplaneSrc = cms.InputTag("hiEvtPlaneFlat")
)

#eventinfoana_mc = eventinfoana.clone()
