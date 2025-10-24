// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// constants, enums and typedefs
//

#define PI 3.1416
#define MAXTRG 1024
#define MAXSEL 100

//
// class decleration
//

class EventInfoTreeProducer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit EventInfoTreeProducer(const edm::ParameterSet&);
  ~EventInfoTreeProducer();

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&) {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initTree();

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* EventInfoNtuple;

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  //float trigPrescale[MAXTRG];
  short centrality;
  int   Ntrkoffline;
  int   NtrkHP;
  int   Npixel;
  short nPV;
  uint candSize;
  //bool  trigHLT[MAXTRG];
  //bool  evtSel[MAXSEL];
  float HFsumETPlus;
  float HFsumETMinus;
  float HFsumET;
  float ZDCPlus;
  float ZDCMinus;
  float ZDCsumE;
  float bestvx;
  float bestvy;
  float bestvz;
  float bestvzError;
  float bestvxError;
  float bestvyError;
  float ephfpAngle[3];
  float ephfmAngle[3];
  float ephfpQ[3];
  float ephfmQ[3];
  //float ephfpQ3[3];
  //float ephfmQ3[3];
  float eptkAngle[2];
  float eptkQ[2];
  float ephfpSumW;
  float ephfmSumW;
  float eptkSumW;
  float ephfpSumW2;
  float ephfmSumW2;
  float eptkSumW2;
  
  bool isCentrality_;
  bool isEventPlane_;

  //token
  edm::EDGetTokenT<reco::BeamSpot> tok_offlineBS_;
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tok_tracks_;
  
};

//
// static data member definitions
//

//
// constructors and destructor
//


EventInfoTreeProducer::EventInfoTreeProducer(const edm::ParameterSet& iConfig)
{
  //input tokens
  tok_offlineBS_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotSrc"));
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
  tok_tracks_ = consumes<std::vector<pat::PackedCandidate>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("TrackCollection")));


  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);

  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  isEventPlane_ = (iConfig.exists("isEventPlane") ? iConfig.getParameter<bool>("isEventPlane") : false);
  if(isEventPlane_) tok_eventplaneSrc_ = consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("eventplaneSrc"));

}


EventInfoTreeProducer::~EventInfoTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EventInfoTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  fillRECO(iEvent, iSetup);
  EventInfoNtuple->Fill();
}


void
EventInfoTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collection
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(tok_offlineBS_, beamspot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  if(!vertices.isValid()) throw cms::Exception("EventInfoTreeProducer") << "Primary vertices  collection not found!" << std::endl;

  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lsNb = iEvent.luminosityBlock();


  centrality = -1;
  if(isCentrality_)
  {
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(tok_centSrc_, cent);
    HFsumETPlus = (cent.isValid() ? cent->EtHFtowerSumPlus() : -1.);
    HFsumETMinus = (cent.isValid() ? cent->EtHFtowerSumMinus() : -1.);
    HFsumET= (cent.isValid() ? cent->EtHFtowerSum() : -1.);
    Npixel = (cent.isValid() ? cent->multiplicityPixel() : -1);
    //ZDCPlus = (cent.isValid() ? cent->zdcSumPlus() : -1.);
    //ZDCMinus = (cent.isValid() ? cent->zdcSumMinus() : -1.);
    //ZDCsumE = (cent.isValid() ? cent->zdcSum() : -1.);

    ZDCPlus = cent->zdcSumPlus();
    ZDCMinus = cent->zdcSumMinus();
    ZDCsumE = cent->zdcSum();
    Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);
      
    edm::Handle<int> cbin;
    iEvent.getByToken(tok_centBinLabel_, cbin);
    centrality = (cbin.isValid() ? *cbin : -1);
  }

  NtrkHP = -1;
  edm::Handle<pat::PackedCandidateCollection> trackColl;
  iEvent.getByToken(tok_tracks_, trackColl);
  if (trackColl.isValid()) {

    NtrkHP = 0;
    
    const reco::TrackBase::TrackQuality highPurity = reco::TrackBase::qualityByName("highPurity");

    for (const auto& trk : *trackColl) {
        if (trk.charge() == 0) continue;
        if (trk.pseudoTrack().quality(highPurity)) {
            NtrkHP++;
        }
    }
  }
  
  if(isEventPlane_)
  {
    edm::Handle<reco::EvtPlaneCollection> eventplanes;
    iEvent.getByToken(tok_eventplaneSrc_, eventplanes);

    //unsigned int collection_size = (*eventplanes).size();
    //std::cout<<"EvtPlaneCollection size="<<collection_size<<std::endl;
    
    const reco::EvtPlane & ephfp1 = (*eventplanes)[15];
    const reco::EvtPlane & ephfm1 = (*eventplanes)[15];
    const reco::EvtPlane & ephfp2 = (*eventplanes)[1];
    const reco::EvtPlane & ephfm2 = (*eventplanes)[0];
    const reco::EvtPlane & ephfp3 = (*eventplanes)[7];
    const reco::EvtPlane & ephfm3 = (*eventplanes)[6];
    const reco::EvtPlane & eptk2 = (*eventplanes)[3];
    const reco::EvtPlane & eptk3 = (*eventplanes)[9];
    
    ephfpAngle[0] = (eventplanes.isValid() ? ephfp1.angle(2) : -99.);
    ephfpAngle[1] = (eventplanes.isValid() ? ephfp2.angle(2) : -99.);
    ephfpAngle[2] = (eventplanes.isValid() ? ephfp3.angle(2) : -99.);

    ephfmAngle[0] = (eventplanes.isValid() ? ephfm1.angle(2) : -99.);
    ephfmAngle[1] = (eventplanes.isValid() ? ephfm2.angle(2) : -99.);
    ephfmAngle[2] = (eventplanes.isValid() ? ephfm3.angle(2) : -99.);

    ephfpQ[0] = (eventplanes.isValid() ? ephfp1.q(2) : -99.);
    ephfpQ[1] = (eventplanes.isValid() ? ephfp2.q(2) : -99.);
    ephfpQ[2] = (eventplanes.isValid() ? ephfp3.q(2) : -99.);

    ephfmQ[0] = (eventplanes.isValid() ? ephfm1.q(2) : -99.);
    ephfmQ[1] = (eventplanes.isValid() ? ephfm2.q(2) : -99.);
    ephfmQ[2] = (eventplanes.isValid() ? ephfm3.q(2) : -99.);
    

    eptkAngle[0] = (eventplanes.isValid() ? eptk2.angle(2): -99.);
    eptkAngle[1] = (eventplanes.isValid() ? eptk3.angle(2): -99.);
    
    eptkQ[0] = (eventplanes.isValid() ? eptk2.q(2): -99.);
    eptkQ[1] = (eventplanes.isValid() ? eptk3.q(2): -99.);

    ephfpSumW = (eventplanes.isValid() ? ephfp2.sumw() : -99.);
    ephfmSumW = (eventplanes.isValid() ? ephfm2.sumw() : -99.);
    eptkSumW  = (eventplanes.isValid() ? eptk2.sumw() : -99.);

    ephfpSumW2 = (eventplanes.isValid() ? ephfp2.sumw2() : -99.);
    ephfmSumW2 = (eventplanes.isValid() ? ephfm2.sumw2() : -99.);
    eptkSumW2  = (eventplanes.isValid() ? eptk2.sumw2() : -99.);
  }

  nPV = vertices->size();
  
  bestvz = -500.9, bestvx = -500.9, bestvy = -500.9;
  bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
  
  const reco::Vertex &vtx = (*vertices)[0]; 
  bestvz = vtx.z();
  bestvx = vtx.x();
  bestvy = vtx.y();
  bestvzError = vtx.zError();
  bestvxError = vtx.xError();
  bestvyError = vtx.yError();

}

// ------------ method called once each job just before starting event
//loop  ------------
void
EventInfoTreeProducer::beginJob()
{
  TH1D::SetDefaultSumw2();

  initTree();
}

void 
EventInfoTreeProducer::initTree()
{ 
  EventInfoNtuple = fs->make<TTree>("EventInfoNtuple","EventInfoNtuple");

  // Event info
  EventInfoNtuple->Branch("RunNb",&runNb,"RunNb/i");
  EventInfoNtuple->Branch("LSNb",&lsNb,"LSNb/i");
  EventInfoNtuple->Branch("EventNb",&eventNb,"EventNb/i");
  EventInfoNtuple->Branch("nPV",&nPV,"nPV/S");
  EventInfoNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
  EventInfoNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
  EventInfoNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
  EventInfoNtuple->Branch("bestvxError",&bestvxError,"bestvxError/F");
  EventInfoNtuple->Branch("bestvyError",&bestvyError,"bestvyError/F");
  EventInfoNtuple->Branch("bestvzError",&bestvzError,"bestvzError/F");

  if(isCentrality_) 
  {
    EventInfoNtuple->Branch("centrality",&centrality,"centrality/S");
    EventInfoNtuple->Branch("Npixel",&Npixel,"Npixel/I");
    EventInfoNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
    EventInfoNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
    EventInfoNtuple->Branch("HFsumET",&HFsumET,"HFsumET/F");
    EventInfoNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
    EventInfoNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
    EventInfoNtuple->Branch("ZDCsumE",&ZDCsumE,"ZDCsumE/F");
    EventInfoNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
    EventInfoNtuple->Branch("NtrkHP",&NtrkHP,"NtrkHP/I");
  }
  if(isEventPlane_) {
    EventInfoNtuple->Branch("ephfpAngle",ephfpAngle,"ephfpAngle[3]/F");
    EventInfoNtuple->Branch("ephfmAngle",ephfmAngle,"ephfmAngle[3]/F");
    EventInfoNtuple->Branch("ephfpQ",ephfpQ,"ephfpQ[3]/F");
    EventInfoNtuple->Branch("ephfmQ",ephfmQ,"ephfmQ[3]/F");
    EventInfoNtuple->Branch("ephfpSumW",&ephfpSumW,"ephfpSumW/F");
    EventInfoNtuple->Branch("ephfmSumW",&ephfmSumW,"ephfmSumW/F");
    EventInfoNtuple->Branch("ephfpSumW2",&ephfpSumW2,"ephfpSumW2/F");
    EventInfoNtuple->Branch("ephfmSumW2",&ephfmSumW2,"ephfmSumW2/F");
    EventInfoNtuple->Branch("eptkAngle",eptkAngle,"eptkAngle[2]/F");
    EventInfoNtuple->Branch("eptkQ",eptkQ,"eptkQ[2]/F");
    EventInfoNtuple->Branch("eptkSumW",&eptkSumW,"eptkSumW/F");
  }
  //EventInfoNtuple->Branch("trigPrescale",trigPrescale,Form("trigPrescale[%d]/F",NTRG_));
  //EventInfoNtuple->Branch("trigHLT",trigHLT,Form("trigHLT[%d]/O",NTRG_));
  //EventInfoNtuple->Branch("evtSel",evtSel,Form("evtSel[%d]/O",NSEL_));
}


//--------------------------------------------------------------------------------------------------
void 
EventInfoTreeProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  /*bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(tok_triggerResults_, triggerResultsLabel);
  hltPrescaleProvider_.init(iRun, iSetup, triggerResultsLabel.process, changed);
  */
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
EventInfoTreeProducer::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventInfoTreeProducer);
