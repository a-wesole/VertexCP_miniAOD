// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <Math/Functions.h>
#include <Math/SVector.h>   
#include <Math/SMatrix.h>

#define PI 3.1416
#define MAXCAN 50000

using namespace std;

class VCTreeProducer_D02kpi : public edm::one::EDAnalyzer<>
{
public:
	explicit VCTreeProducer_D02kpi(const edm::ParameterSet &);
	~VCTreeProducer_D02kpi();

	using MVACollection = std::vector<float>;

	void resize_the_vectors(int newSize);

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void fillRECO(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();
	virtual void initTree();
	void genDecayLength(const uint &, const reco::GenParticle &);

	// ----------member data ---------------------------

	edm::Service<TFileService> fs;

	TTree *VertexCompositeNtuple;
	bool saveTree_;

	// options
	bool doRecoNtuple_;
	bool dogenntuple_;
	bool dogenmatching_;
	bool dogenmatchingtof_;
	bool hasswap_;
	bool decayingen_;
	int PID_;
	int PID_dau1_;
	int PID_dau2_;

	// cut variables
	double multMax_;
	double multMin_;
	double deltaR_; // deltaR for Gen matching

	// tree branches
	// event info
	int centrality;
	int Ntrkoffline;
	int Npixel;
	float HFsumETPlus;
	float HFsumETMinus;
	float ZDCPlus;
	float ZDCMinus;
	float bestvx;
	float bestvy;
	float bestvz;
	float bestvxError;
	float bestvyError;
	float bestvzError;
	float BSx;
	float BSy;
	float BSz;
	float BSxerror;
	float BSyerror;
	float BSzerror;
	int candSize;
	float ephfpAngle[3];
	float ephfmAngle[3];
	float ephfpQ[3];
	float ephfmQ[3];
	float ephfpSumW;
	float ephfmSumW;

	// Composite candidate info
	std::vector<float> mva;
	std::vector<float> mva_xg;
	std::vector<float> pt;
	std::vector<float> eta;
	std::vector<float> phi;
	std::vector<float> flavor;
	std::vector<float> y;
	std::vector<float> mass;
	std::vector<float> VtxProb;
	std::vector<float> dlos;
	std::vector<float> dl;
	std::vector<float> dlerror;
	std::vector<float> DlxyBS;
	std::vector<float> DlxyBSErr;
	std::vector<float> vtxChi2;
	std::vector<float> ndf;
	std::vector<float> agl_abs;
	std::vector<float> ip3d;
	std::vector<float> ip3derr;
	std::vector<float> agl2D_abs;
	std::vector<float> dlos2D;
	std::vector<float> dl2D;
	std::vector<float> dl2Derror;
	std::vector<bool> isSwap;
	std::vector<bool> matchGEN;
	std::vector<int> idmom_reco;
	std::vector<int> idd1_reco;
	std::vector<int> idd2_reco;
	std::vector<float> gen_agl_abs;
	std::vector<float> gen_agl2D_abs;
	std::vector<float> gen_dl;
	std::vector<float> gen_dl2D;
	std::vector<float> twoTrackDCA;

	// dau info
	std::vector<float> dzos1;
	std::vector<float> dzos2;
	std::vector<float> dxyos1;
	std::vector<float> dxyos2;
	std::vector<float> pt1;
	std::vector<float> pt2;
	std::vector<float> ptErr1;
	std::vector<float> ptErr2;
	std::vector<float> p1;
	std::vector<float> p2;
	std::vector<float> Dtrk1Dz1;
	std::vector<float> Dtrk2Dz1;
	std::vector<float> Dtrk1Dxy1;
	std::vector<float> Dtrk2Dxy1;
	std::vector<float> Dtrk1DzError1;
	std::vector<float> Dtrk2DzError1;
	std::vector<float> Dtrk1DxyError1;
	std::vector<float> Dtrk2DxyError1;
	std::vector<float> eta1;
	std::vector<float> eta2;
	std::vector<float> phi1;
	std::vector<float> phi2;
	std::vector<int> charge1;
	std::vector<int> charge2;
	std::vector<int> pid1;
	std::vector<int> pid2;
	std::vector<float> tof1;
	std::vector<float> tof2;
	std::vector<float> H2dedx1;
	std::vector<float> H2dedx2;
	std::vector<float> T4dedx1;
	std::vector<float> T4dedx2;
	std::vector<float> trkChi1;
	std::vector<float> trkChi2;

	// gen info
	std::vector<float> pt_gen;
	std::vector<float> eta_gen;
	std::vector<int> idmom;
	std::vector<float> y_gen;
	std::vector<float> phi_gen;
	std::vector<int> iddau1;
	std::vector<int> iddau2;

	bool useAnyMVA_;
	bool isSkimMVA_;
	bool isCentrality_;
	bool doGenNtuple_;
	bool doGenMatching_;
	bool decayInGen_;

	edm::Handle<int> cbin_;

	// tokens
	edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
	// edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
	edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tok_generalTrk_;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;

	edm::EDGetTokenT<MVACollection> MVAValues_Token_;
	edm::EDGetTokenT<MVACollection> MVAValues_Token_2;

	edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token1_;
	edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token2_;
	edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;

	edm::EDGetTokenT<int> tok_centBinLabel_;
	edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

	edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

	edm::EDGetTokenT<reco::BeamSpot> bsLabel_;
};

void VCTreeProducer_D02kpi::resize_the_vectors(int newSize)
{
	mva.resize(newSize);
	mva_xg.resize(newSize);
	pt.resize(newSize);
	eta.resize(newSize);
	phi.resize(newSize);
	flavor.resize(newSize);
	y.resize(newSize);
	mass.resize(newSize);
	VtxProb.resize(newSize);
	dlos.resize(newSize);
	dl.resize(newSize);
	dlerror.resize(newSize);
	DlxyBS.resize(newSize);
	DlxyBSErr.resize(newSize);
	vtxChi2.resize(newSize);
	ndf.resize(newSize);
	agl_abs.resize(newSize);
	ip3d.resize(newSize);
	ip3derr.resize(newSize);
	agl2D_abs.resize(newSize);
	dlos2D.resize(newSize);
	dl2D.resize(newSize);
	dl2Derror.resize(newSize);
	isSwap.resize(newSize);
	matchGEN.resize(newSize);
	idmom_reco.resize(newSize);
	idd1_reco.resize(newSize);
	idd2_reco.resize(newSize);
	gen_agl_abs.resize(newSize);
	gen_agl2D_abs.resize(newSize);
	gen_dl.resize(newSize);
	gen_dl2D.resize(newSize);
	twoTrackDCA.resize(newSize);

	dzos1.resize(newSize);
	dzos2.resize(newSize);
	dxyos1.resize(newSize);
	dxyos2.resize(newSize);
	pt1.resize(newSize);
	pt2.resize(newSize);
	ptErr1.resize(newSize);
	ptErr2.resize(newSize);
	p1.resize(newSize);
	p2.resize(newSize);
	Dtrk1Dz1.resize(newSize);
	Dtrk2Dz1.resize(newSize);
	Dtrk1Dxy1.resize(newSize);
	Dtrk2Dxy1.resize(newSize);
	Dtrk1DzError1.resize(newSize);
	Dtrk2DzError1.resize(newSize);
	Dtrk1DxyError1.resize(newSize);
	Dtrk2DxyError1.resize(newSize);
	eta1.resize(newSize);
	eta2.resize(newSize);
	phi1.resize(newSize);
	phi2.resize(newSize);
	charge1.resize(newSize);
	charge2.resize(newSize);
	pid1.resize(newSize);
	pid2.resize(newSize);
	tof1.resize(newSize);
	tof2.resize(newSize);
	H2dedx1.resize(newSize);
	H2dedx2.resize(newSize);
	T4dedx1.resize(newSize);
	T4dedx2.resize(newSize);
	trkChi1.resize(newSize);
	trkChi2.resize(newSize);

	pt_gen.resize(newSize);
	eta_gen.resize(newSize);
	idmom.resize(newSize);
	y_gen.resize(newSize);
	phi_gen.resize(newSize);
	iddau1.resize(newSize);
	iddau2.resize(newSize);
}

VCTreeProducer_D02kpi::VCTreeProducer_D02kpi(const edm::ParameterSet &iConfig)
{

	// options
	doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
	doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
	doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
	decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
	PID_ = iConfig.getUntrackedParameter<int>("PID");
	PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
	PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");

	saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");

	useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
	isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA");
	isCentrality_ = iConfig.getParameter<bool>("isCentrality");

	// cut variables
	multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
	multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
	deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.02);

	// input tokens
	tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
	// tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
	tok_generalTrk_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));

	patCompositeCandidateCollection_Token_ = consumes<pat::CompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("D0"));
	MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
	Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxHarmonic2"));
	Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxTruncated40"));
	tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));
	bsLabel_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BSLabel"));

	if (isCentrality_)
	{
		tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
		tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
	}

	if (useAnyMVA_ && iConfig.exists("MVACollection") && iConfig.exists("MVACollection2"))
	{
		MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
		MVAValues_Token_2 = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection2"));
	}
}

VCTreeProducer_D02kpi::~VCTreeProducer_D02kpi()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to for each event  ------------
void VCTreeProducer_D02kpi::analyze(const edm::Event &iEvent, const edm::EventSetup &
																		iSetup)
{
	using std::vector;
	using namespace edm;
	using namespace reco;

	if (doRecoNtuple_) {
		fillRECO(iEvent, iSetup);
	}

	if (saveTree_)
		VertexCompositeNtuple->Fill();
}

void VCTreeProducer_D02kpi::fillRECO(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
#ifdef DEBUG
	using std::cout;
	using std::endl;
#endif
	// get collections
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(tok_offlinePV_, vertices);

	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(bsLabel_, beamSpotHandle);

	// edm::Handle<reco::TrackCollection> tracks;
	edm::Handle<pat::PackedCandidateCollection> tracks;
	iEvent.getByToken(tok_generalTrk_, tracks);

	edm::Handle<pat::CompositeCandidateCollection> D0candidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_, D0candidates);
	const pat::CompositeCandidateCollection *D0candidates_ = D0candidates.product();

	edm::Handle<MVACollection> mvavalues;
	edm::Handle<MVACollection> mvavalues_xg;
	if (useAnyMVA_)
	{
		iEvent.getByToken(MVAValues_Token_, mvavalues);
		iEvent.getByToken(MVAValues_Token_2, mvavalues_xg);
		assert((*mvavalues).size() == D0candidates->size());
		assert((*mvavalues_xg).size() == D0candidates->size());
	}

	edm::Handle<reco::GenParticleCollection> genpars;
	if (doGenMatching_)
		iEvent.getByToken(tok_genParticle_, genpars);

	edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle1;
	iEvent.getByToken(Dedx_Token1_, dEdxHandle1);

	edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle2;
	iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
#ifdef DEBUG
	cout << "Loaded tokens" << endl;
#endif

	centrality = -1;
	if (isCentrality_)
	{
		edm::Handle<reco::Centrality> cent;
		iEvent.getByToken(tok_centSrc_, cent);
		HFsumETPlus = (cent.isValid() ? cent->EtHFtowerSumPlus() : -1.);
		HFsumETMinus = (cent.isValid() ? cent->EtHFtowerSumMinus() : -1.);
		Npixel = (cent.isValid() ? cent->multiplicityPixel() : -1);
		ZDCPlus = (cent.isValid() ? cent->zdcSumPlus() : -1.);
		ZDCMinus = (cent.isValid() ? cent->zdcSumMinus() : -1.);
		Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);

		edm::Handle<int> cbin;
		iEvent.getByToken(tok_centBinLabel_, cbin);
		centrality = (cbin.isValid() ? *cbin : -1);
	}

	// best vertex

	BSx = -999.9;
	BSy = -999.9;
	BSz = -999.9;
	BSxerror = -999.9;
	BSyerror = -999.9;
	BSzerror = -999.9;
	reco::BeamSpot beamSpot;

	float BSdxdz = -999.9;
	float BSdydz = -999.9;
	bestvz = -500.9;
	bestvx = -500.9;
	bestvy = -500.9;
	bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;

	const reco::Vertex &vtx = (*vertices)[0]; // we need ot kwwp this for the vtx covariance, but that is all!
	bestvz = vtx.z();
	bestvx = vtx.x();
	bestvy = vtx.y();
	bestvzError = vtx.zError();
	bestvxError = vtx.xError();
	bestvyError = vtx.yError();
	float xlxyBS = -999.9;
	float ylxyBS = -999.9;
	float vtxYXErr = -999.9;
	float vtxXErr = -999.9;
	float vtxYErr = -999.9;

	beamSpot = *beamSpotHandle;
	reco::Vertex theBeamSpotV(beamSpot.position(), beamSpot.covariance3D());
	BSx = beamSpot.x0();
	BSy = beamSpot.y0();
	BSz = beamSpot.z0();
	BSxerror = beamSpot.x0Error();
	BSyerror = beamSpot.y0Error();
	BSzerror = beamSpot.z0Error();
	BSdxdz = beamSpot.dxdz();
	BSdydz = beamSpot.dydz();

	// Ntrkoffline
	Ntrkoffline = 0;
#ifdef DEBUG
	cout << "Calculated offline Ntrk's" << endl;
#endif

#ifdef DEBUG
	cout << "Gen matching done" << endl;
#endif

	// RECO Candidate info
	candSize = D0candidates_->size();
	resize_the_vectors(candSize);

	//  cout << "for event number " << iEvent.id().event() << " found " << candSize << " D0 candidates in tree producer" << endl;
	for (unsigned it = 0; it < D0candidates_->size(); ++it)
	{

		const pat::CompositeCandidate &trk = (*D0candidates_)[it];
		double secvz = -999.9, secvx = -999.9, secvy = -999.9;
		secvz = trk.userFloat("vtxZ");
		secvx = trk.userFloat("vtxX");
		secvy = trk.userFloat("vtxY");
		bestvz = trk.userFloat("bestvtxZ");
		bestvx = trk.userFloat("bestvtxX");
		bestvy = trk.userFloat("bestvtxY");
		bestvzError = trk.userFloat("zVtxError");
		bestvxError = trk.userFloat("xVtxError");
		bestvyError = trk.userFloat("yVtxError");

		reco::Vertex::CovarianceMatrix sec_covariance;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				sec_covariance(i, j) = trk.userFloat("vertexCovariance_" + std::to_string(i) + "_" + std::to_string(j));
			}
		}

		vtxYXErr = sec_covariance(1, 0);
		vtxXErr = sec_covariance(0, 0);
		vtxYErr = sec_covariance(1, 1);

		eta[it] = trk.eta();
		y[it] = trk.rapidity();
		pt[it] = trk.pt();

		phi[it] = trk.phi();
		flavor[it] = trk.pdgId() / abs(trk.pdgId());
		twoTrackDCA[it] = trk.userFloat("track3DDCA");

		mva[it] = 0.0;
		if (useAnyMVA_)
		{
			mva[it] = (*mvavalues)[it];
			mva_xg[it] = (*mvavalues_xg)[it]; // xgboost
		}

		double px = trk.px();
		double py = trk.py();
		double pz = trk.pz();
		mass[it] = trk.mass();

		if (mass[it] == 0)
		{
			cout << "Error break" << endl;
		}

		// const reco::Candidate * reco_d1 = trk.daughter(0);
		// const reco::Candidate * reco_d2 = trk.daughter(1);

		const reco::Candidate *cand1 = trk.daughter(0);
		const pat::PackedCandidate *reco_d1 = dynamic_cast<const pat::PackedCandidate *>(cand1);

		const reco::Candidate *cand2 = trk.daughter(1);
		const pat::PackedCandidate *reco_d2 = dynamic_cast<const pat::PackedCandidate *>(cand2);

		int reco_d1_charge = TMath::Sign(1, reco_d1->pdgId());
		int reco_d2_charge = TMath::Sign(1, reco_d2->pdgId());

		// const pat::PackedCandidate* reco_d3 = nullptr;

		// Gen match

		matchGEN[it] = false;
		isSwap[it] = false;
		idmom_reco[it] = -77;
		idd1_reco[it] = -77;
		idd2_reco[it] = -77;

		pt_gen[it] = -999.9;
		eta_gen[it] = -999.9;
		idmom[it] = -999;
		y_gen[it] = -999.9;
		phi_gen[it] = -999.9;
		iddau1[it] = -999;
		iddau2[it] = -999;

		if (doGenMatching_)
		{ // doGenMatching

			if (!genpars.isValid())
			{
				cout << "Gen matching cannot be done without Gen collection!!" << endl;
				return;
			}

			for (unsigned genPair = 0; genPair < genpars->size(); ++genPair)
			{ // loop over all gen particles ->known D0 to kPi pairs

				const reco::GenParticle &genD0 = (*genpars)[genPair];

				int id = genD0.pdgId();
				if (fabs(id) != PID_)
					continue; // check to make sure is D0

				const reco::Candidate *gen_d1 = genD0.daughter(0);
				const reco::Candidate *gen_d2 = genD0.daughter(1);

				if (!(fabs(gen_d1->pdgId()) == PID_dau1_ && fabs(gen_d2->pdgId()) == PID_dau2_) && !(fabs(gen_d2->pdgId()) == PID_dau1_ && fabs(gen_d1->pdgId()) == PID_dau2_))
					continue; // make sure k pi pairs
				// cout << "reco d1 pid =" << reco_d1->pdgId() << " and 2 " << reco_d2->pdgId() << endl;
				// cout << "gen d1 pid =" << gen_d1->pdgId() << " and 2 " << gen_d2->pdgId() << endl;

				// if(((reco_d1->charge() == gen_d1->charge() && reco_d2->charge() == gen_d2->charge()) || (reco_d1->charge() == gen_d2->charge() && reco_d2->charge() == gen_d1->charge()))) {
				if (((reco_d1_charge == gen_d1->charge() && reco_d2_charge == gen_d2->charge()) || (reco_d1_charge == gen_d2->charge() && reco_d2_charge == gen_d1->charge())))
				{
					// cout << "okay we passed the charge statement " << endl;

					if (reco_d1_charge == gen_d1->charge())
					{
						double deltaR = sqrt(pow(reco_d1->eta() - gen_d1->eta(), 2) + pow(reco_d1->phi() - gen_d1->phi(), 2));
						if (deltaR > deltaR_)
							continue; // check deltaR matching
						if (fabs((reco_d1->pt() - gen_d1->pt()) / reco_d1->pt()) > 0.2)
							continue; // check deltaPt matching

						deltaR = sqrt(pow(reco_d2->eta() - gen_d2->eta(), 2) + pow(reco_d2->phi() - gen_d2->phi(), 2));
						if (deltaR > deltaR_)
							continue; // check deltaR matching
						if (fabs((reco_d2->pt() - gen_d2->pt()) / reco_d2->pt()) > 0.2)
							continue; // check deltaPt matching

						cout << "SIIIIIIIGNAL -- let's go baby!!" << endl;
						matchGEN[it] = true; // matched gen
						if (reco_d1->pdgId() != gen_d1->pdgId())
							isSwap[it] = true;
						genDecayLength(it, genD0);

						pt_gen[it] = genD0.pt();
						eta_gen[it] = genD0.eta();
						y_gen[it] = genD0.rapidity();
						phi_gen[it] = genD0.phi();

						idmom[it] = genD0.pdgId();

						if (!decayInGen_)
							continue;

						iddau1[it] = gen_d1->pdgId();
						iddau2[it] = gen_d2->pdgId();

						break;
					}

					if (reco_d1->charge() == gen_d2->charge())
					{
						double deltaR = sqrt(pow(reco_d1->eta() - gen_d2->eta(), 2) + pow(reco_d1->phi() - gen_d2->phi(), 2));
						if (deltaR > deltaR_)
							continue; // check deltaR matching
						if (fabs((reco_d1->pt() - gen_d2->pt()) / reco_d1->pt()) > 0.2)
							continue; // check deltaPt matching

						deltaR = sqrt(pow(reco_d2->eta() - gen_d1->eta(), 2) + pow(reco_d2->phi() - gen_d1->phi(), 2));
						if (deltaR > deltaR_)
							continue; // check deltaR matching
						if (fabs((reco_d2->pt() - gen_d1->pt()) / reco_d2->pt()) > 0.2)
							continue; // check deltaPt matching

						cout << "SIIIIIIIGNAL -- let's go baby!!" << endl;
						matchGEN[it] = true; // matched gen
						if (reco_d1->pdgId() != gen_d2->pdgId())
							isSwap[it] = true;
						genDecayLength(it, genD0);

						pt_gen[it] = genD0.pt();
						eta_gen[it] = genD0.eta();
						y_gen[it] = genD0.rapidity();
						phi_gen[it] = genD0.phi();

						idmom[it] = genD0.pdgId();

						if (!decayInGen_)
							continue;

						iddau1[it] = gen_d1->pdgId();
						iddau2[it] = gen_d2->pdgId();

						break;
					}
				}
				else
				{

					//  cout << "did not pass same charge! error! reco d1 charge =" << reco_d1->charge() << " and 2 " << reco_d2->charge() << endl;
					// cout << "gen d1 charge =" << gen_d1->charge() << " and 2 " << gen_d2->charge() << endl;
				}
			} // loop over all gen particles -- to find known D0->kPi pairs
			idmom_reco[it] = trk.pdgId();
			idd1_reco[it] = reco_d1->pdgId();
			idd2_reco[it] = reco_d2->pdgId();

		} // doGenMatching

		double pxd1 = reco_d1->px();
		double pyd1 = reco_d1->py();
		double pzd1 = reco_d1->pz();
		double pxd2 = reco_d2->px();
		double pyd2 = reco_d2->py();
		double pzd2 = reco_d2->pz();

		TVector3 dauvec1(pxd1, pyd1, pzd1);
		TVector3 dauvec2(pxd2, pyd2, pzd2);

		// pt
		pt1[it] = reco_d1->pt();
		pt2[it] = reco_d2->pt();

		//momentum
		p1[it] = reco_d1->p();
		p2[it] = reco_d2->p();

		//eta
		eta1[it] = reco_d1->eta();
		eta2[it] = reco_d2->eta();

		//phi
		phi1[it] = reco_d1->phi();
		phi2[it] = reco_d2->phi();

		//charge
		// charge1[it] = reco_d1->charge();
		// charge2[it] = reco_d2->charge();
		charge1[it] = reco_d1_charge;
		charge2[it] = reco_d2_charge;


		pid1[it] = -99999;
		pid2[it] = -99999;

		//vtxChi2
		vtxChi2[it] = trk.userFloat("vtxChi2");
		ndf[it] = trk.userFloat("vtxNdof");
		VtxProb[it] = TMath::Prob(vtxChi2[it],ndf[it]);

		// PAngle
		TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
		TVector3 secvec(px, py, pz);

		TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
		TVector3 secvec2D(px, py, 0);

        agl_abs[it] = secvec.Angle(ptosvec);

		agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

		float r2lxyBS = (secvx - BSx - (secvz - BSz) * BSdxdz) * (secvx - BSx - (secvz - BSz) * BSdxdz) + (secvy - BSy - (secvz - BSz) * BSdydz) * (secvy - BSy - (secvz - BSz) * BSdydz);
		xlxyBS = secvx - BSx - (secvz - BSz) * BSdxdz;
		ylxyBS = secvy - BSy - (secvz - BSz) * BSdydz;
		DlxyBS[it] = static_cast<float>(TMath::Sqrt(r2lxyBS));
		DlxyBSErr[it] = static_cast<float>(TMath::Sqrt((1. / r2lxyBS) * ((xlxyBS * xlxyBS) * vtxXErr + (2 * xlxyBS * ylxyBS) * vtxYXErr + (ylxyBS * ylxyBS) * vtxYErr)));

		// Decay length 3D
		typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
		typedef ROOT::Math::SVector<double, 3> SVector3;
		typedef ROOT::Math::SVector<double, 6> SVector6;

		SMatrixSym3D totalCov = vtx.covariance() + sec_covariance;
		SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

		dl[it] = ROOT::Math::Mag(distanceVector);
		dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / dl[it];

		dlos[it] = dl[it] / dlerror[it];

		// Decay length 2D
		SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1), vtx.covariance(1, 1), 0, 0, 0);
		SVector6 v2(sec_covariance(0, 0), sec_covariance(0, 1), sec_covariance(1, 1), 0, 0, 0);

		SMatrixSym3D sv1(v1);
		SMatrixSym3D sv2(v2);

		SMatrixSym3D totalCov2D = sv1 + sv2;
		SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

		dl2D[it] = ROOT::Math::Mag(distanceVector2D);
		double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D[it];

		dlos2D[it] = dl2D[it]/dl2Derror;

		//trk info
		//auto dau1 = reco_d1->get<reco::TrackRef>();

		const pat::PackedCandidate* dau1 = dynamic_cast<const pat::PackedCandidate*>(reco_d1);
		const reco::Track& pseudoTrk1 = dau1->pseudoTrack();
		
		trkChi1[it] = pseudoTrk1.normalizedChi2();
		ptErr1[it] = pseudoTrk1.ptError();

		// nhit1 = pseudoTrk1.hitPattern().numberOfValidHits();
		// trkquality1 = pseudoTrk1.quality(reco::TrackBase::highPurity);

		// trk quality

		// trk dEdx
		/*H2dedx1[it] = -999.9;

		if(dEdxHandle1.isValid()){
			const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
			H2dedx1[it] = dEdxTrack[dau1].dEdx();
		}

		T4dedx1[it] = -999.9;

		if(dEdxHandle2.isValid()){
			const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
			T4dedx1[it] = dEdxTrack[dau1].dEdx();
			}*/

		// track Chi2
		// trkChi1[it] = dau1->normalizedChi2();

		// track pT error
		// ptErr1[it] = dau1->ptError();

		// vertexCovariance 00-xError 11-y 22-z
		secvz = trk.vz();
		secvx = trk.vx();
		secvy = trk.vy();

		// trkNHits

		// DCA
		math::XYZPoint bestvtx(bestvx, bestvy, bestvz);
		math::XYZPoint BS_vtx(BSx, BSy, BSz);

		double dzbest1 = dau1->pseudoTrack().dz(bestvtx);	// today change
		double dxybest1 = dau1->pseudoTrack().dxy(bestvtx); // today change
		double dzerror1 = TMath::Sqrt(dau1->pseudoTrack().dzError() * dau1->pseudoTrack().dzError() + bestvzError * bestvzError);
		double dxyerror1 = TMath::Sqrt(dau1->pseudoTrack().dxyError() * dau1->pseudoTrack().dxyError() + bestvxError * bestvyError);
		// double dzbest1 = dau1->dz(bestvtx);
		// double dxybest1 = dau1->dxy(bestvtx);
		// double dzerror1 = std::sqrt(pseudoTrk1.dzError() * pseudoTrk1.dzError() + bestvzError * bestvzError);
		// double dxyerror1 = std::sqrt(pseudoTrk1.d0Error() * pseudoTrk1.d0Error() + bestvxError * bestvyError);

		Dtrk2Dz1[it] = dzbest1;
		Dtrk2Dxy1[it] = dxybest1;
		Dtrk2DzError1[it] = dzerror1;
		Dtrk2DxyError1[it] = dxyerror1;
		dzos1[it] = dzbest1 / dzerror1;
		dxyos1[it] = dxybest1 / dxyerror1;

		// auto dau2 = reco_d2->get<reco::TrackRef>();
		const pat::PackedCandidate *dau2 = dynamic_cast<const pat::PackedCandidate *>(reco_d2);

		const reco::Track &pseudoTrk2 = dau2->pseudoTrack();
		trkChi2[it] = pseudoTrk2.normalizedChi2();
		ptErr2[it] = pseudoTrk2.ptError();

		// nhit2[it] = pseudoTrk2.hitPattern().numberOfValidHits();
		// trkquality2 = pseudoTrk2.quality(reco::TrackBase::highPurity);

		// trk quality

		// trk dEdx
		/*H2dedx2[it] = -999.9;

		if(dEdxHandle1.isValid()){
		  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
		  H2dedx2[it] = dEdxTrack[dau2].dEdx();
		}

		T4dedx2[it] = -999.9;

		if(dEdxHandle2.isValid()){
		  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
		  T4dedx2[it] = dEdxTrack[dau2].dEdx();
		}
		*/

		// track Chi2
		// trkChi2[it] = dau2->normalizedChi2();

		// track pT error
		// ptErr2[it] = dau2->ptError();

		// vertexCovariance 00-xError 11-y 22-z
		secvz = trk.vz();
		secvx = trk.vx();
		secvy = trk.vy();

		// trkNHits

		// DCA
		double dzbest2 = dau2->pseudoTrack().dz(bestvtx);	// today change
		double dxybest2 = dau2->pseudoTrack().dxy(bestvtx); // today change
		double dzerror2 = TMath::Sqrt(dau2->pseudoTrack().dzError() * dau2->pseudoTrack().dzError() + bestvzError * bestvzError);
		double dxyerror2 = TMath::Sqrt(dau2->pseudoTrack().dxyError() * dau2->pseudoTrack().dxyError() + bestvxError * bestvyError);
		// double dzbest2 = dau2->dz(bestvtx);
		// double dxybest2 = dau2->dxy(bestvtx);
		// //double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
		// //double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
		// double dzerror2 = std::sqrt(pseudoTrk2.dzError() * pseudoTrk2.dzError() + bestvzError * bestvzError);
		// double dxyerror2 = std::sqrt(pseudoTrk2.d0Error() * pseudoTrk2.d0Error() + bestvxError * bestvyError);

		dzos2[it] = dzbest2 / dzerror2;
		dxyos2[it] = dxybest2 / dxyerror2;

		Dtrk1Dz1[it] = 1.0 * dzbest2;
		Dtrk1Dxy1[it] = dxybest2;
		Dtrk1DzError1[it] = dzerror2;
		Dtrk1DxyError1[it] = dxyerror2;

#ifdef DEBUG
		cout << "Done reco single iter" << endl;
#endif
	}
#ifdef DEBUG
	cout << "Fill reco done" << endl;
#endif
}

// ------------ method called once each job just before starting event
// loop  ------------
void VCTreeProducer_D02kpi::beginJob()
{
	TH1D::SetDefaultSumw2();

	if (!doRecoNtuple_ && !doGenNtuple_)
	{
		cout << "No output for either RECO or GEN!! Fix config!!" << endl;
		return;
	}

	if (saveTree_)
		initTree();
}

// void VCTreeProducer_D02kpi::initTree()
// {
// 	VertexCompositeNtuple = fs->make<TTree>("VCNtuple_D02kpi", "VCNtuple_D02kpi");

// 	if (doRecoNtuple_)
// 	{

// 		// Event info
// 		VertexCompositeNtuple->Branch("Ntrkoffline", &Ntrkoffline, "Ntrkoffline/I");
// 		VertexCompositeNtuple->Branch("Npixel", &Npixel, "Npixel/I");
// 		VertexCompositeNtuple->Branch("HFsumETPlus", &HFsumETPlus, "HFsumETPlus/F");
// 		VertexCompositeNtuple->Branch("HFsumETMinus", &HFsumETMinus, "HFsumETMinus/F");
// 		VertexCompositeNtuple->Branch("ZDCPlus", &ZDCPlus, "ZDCPlus/F");
// 		VertexCompositeNtuple->Branch("ZDCMinus", &ZDCMinus, "ZDCMinus/F");
// 		VertexCompositeNtuple->Branch("PvtxX", &bestvx, "PvtxX/F");
// 		VertexCompositeNtuple->Branch("PvtxY", &bestvy, "PvtxY/F");
// 		VertexCompositeNtuple->Branch("PvtxZ", &bestvz, "PvtxZ/F");
// 		VertexCompositeNtuple->Branch("BSx", &BSx, "BSx/F");
// 		VertexCompositeNtuple->Branch("BSy", &BSy, "BSy/F");
// 		VertexCompositeNtuple->Branch("BSz", &BSz, "BSz/F");
// 		VertexCompositeNtuple->Branch("PvtxXErr", &bestvxError, "PvtxXErr/F");
// 		VertexCompositeNtuple->Branch("PvtxYErr", &bestvyError, "PvtxYErr/F");
// 		VertexCompositeNtuple->Branch("PvtxZErr", &bestvzError, "PvtxZErr/F");
// 		VertexCompositeNtuple->Branch("BSxErr", &BSxerror, "BSxErr/F");
// 		VertexCompositeNtuple->Branch("BSyErr", &BSyerror, "BSyErr/F");
// 		VertexCompositeNtuple->Branch("BSzErr", &BSzerror, "BSzErr/F");
// 		VertexCompositeNtuple->Branch("candSize", &candSize, "candSize/I");
// 		if (isCentrality_)
// 			VertexCompositeNtuple->Branch("centrality", &centrality, "centrality/I");
// 		// particle info
// 		VertexCompositeNtuple->Branch("pT", &pt, "pT[candSize]/F");
// 		VertexCompositeNtuple->Branch("y", &y, "y[candSize]/F");
// 		VertexCompositeNtuple->Branch("phi", &phi, "phi[candSize]/F");
// 		VertexCompositeNtuple->Branch("mass", &mass, "mass[candSize]/F");
// 		if (useAnyMVA_)
// 		{
// 			VertexCompositeNtuple->Branch("mva", &mva, "mva[candSize]/F");
// 			VertexCompositeNtuple->Branch("mva_xg", &mva_xg, "mva_xg[candSize]/F");
// 		}

// 		if (!isSkimMVA_)
// 		{
// 			// Composite candidate info RECO
// 			VertexCompositeNtuple->Branch("flavor", &flavor, "flavor[candSize]/F");
// 			VertexCompositeNtuple->Branch("eta", &eta, "eta[candSize]/F");
// 			VertexCompositeNtuple->Branch("VtxProb", &VtxProb, "VtxProb[candSize]/F");
// 			VertexCompositeNtuple->Branch("VtxChi2", &vtxChi2, "VtxChi2[candSize]/F");
// 			VertexCompositeNtuple->Branch("3DPointingAngle", &agl_abs, "3DPointingAngle[candSize]/F");
// 			VertexCompositeNtuple->Branch("2DPointingAngle", &agl2D_abs, "2DPointingAngle[candSize]/F");
// 			VertexCompositeNtuple->Branch("3DDecayLengthSignificance", &dlos, "3DDecayLengthSignificance[candSize]/F");
// 			VertexCompositeNtuple->Branch("3DDecayLength", &dl, "3DDecayLength[candSize]/F");
// 			VertexCompositeNtuple->Branch("3DDecayLengthError", &dlerror, "3DDecayLengthError[candSize]/F");
// 			VertexCompositeNtuple->Branch("2DDecayLengthSignificance", &dlos2D, "2DDecayLengthSignificance[candSize]/F");
// 			VertexCompositeNtuple->Branch("2DDecayLength", &dl2D, "2DDecayLength[candSize]/F");
// 			VertexCompositeNtuple->Branch("2DDecayLengthError", &dl2Derror, "2DDecayLengthError[candSize]/F");
// 			VertexCompositeNtuple->Branch("DlxyBS", &DlxyBS, "DlxyBS[candSize]/F");
// 			VertexCompositeNtuple->Branch("DlxyBSErr", &DlxyBSErr, "DlxyBSErr[candSize]/F");
// 			VertexCompositeNtuple->Branch("ip3d", &ip3d, "ip3d[candSize]/F");
// 			VertexCompositeNtuple->Branch("ip3derr", &ip3derr, "ip3derr[candSize]/F");
// 			VertexCompositeNtuple->Branch("zDCASignificanceDaugther1", &dzos1, "zDCASignificanceDaugther1[candSize]/F");
// 			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1", &dxyos1, "xyDCASignificanceDaugther1[candSize]/F");
// 			VertexCompositeNtuple->Branch("pTD1", &pt1, "pTD1[candSize]/F");
// 			VertexCompositeNtuple->Branch("pTerrD1", &ptErr1, "pTerrD1[candSize]/F");
// 			VertexCompositeNtuple->Branch("EtaD1", &eta1, "EtaD1[candSize]/F");
// 			VertexCompositeNtuple->Branch("PhiD1", &phi1, "PhiD1[candSize]/F");
// 			VertexCompositeNtuple->Branch("dedxHarmonic2D1", &H2dedx1, "dedxHarmonic2D1[candSize]/F");
// 			VertexCompositeNtuple->Branch("zDCASignificanceDaugther2", &dzos2, "zDCASignificanceDaugther2[candSize]/F");
// 			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2", &dxyos2, "xyDCASignificanceDaugther2[candSize]/F");
// 			VertexCompositeNtuple->Branch("pTD2", &pt2, "pTD2[candSize]/F");
// 			VertexCompositeNtuple->Branch("pTerrD2", &ptErr2, "pTerrD2[candSize]/F");
// 			VertexCompositeNtuple->Branch("EtaD2", &eta2, "EtaD2[candSize]/F");
// 			VertexCompositeNtuple->Branch("PhiD2", &phi2, "PhiD2[candSize]/F");
// 			VertexCompositeNtuple->Branch("dedxHarmonic2D2", &H2dedx2, "dedxHarmonic2D2[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk1Dz1", &Dtrk1Dz1, "Dtrk1Dz1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk2Dz1", &Dtrk2Dz1, "Dtrk2Dz1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk1Dxy1", &Dtrk1Dxy1, "Dtrk1Dxy1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk2Dxy1", &Dtrk2Dxy1, "Dtrk2Dxy1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk1DzError1", &Dtrk1DzError1, "Dtrk1DzError1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk2DzError1", &Dtrk2DzError1, "Dtrk2DzError1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk1DxyError1", &Dtrk1DxyError1, "Dtrk1DxyError1[candSize]/F");
// 			VertexCompositeNtuple->Branch("Dtrk2DxyError1", &Dtrk2DxyError1, "Dtrk2DxyError1[candSize]/F");
// 			VertexCompositeNtuple->Branch("twoTrackDCA", &twoTrackDCA, "twoTrackDCA[candSize]/F");

// 			if (doGenMatching_)
// 			{
// 				VertexCompositeNtuple->Branch("isSwap", &isSwap, "isSwap[candSize]/O");
// 				VertexCompositeNtuple->Branch("idmom_reco", &idmom_reco, "idmom_reco[candSize]/I");
// 				VertexCompositeNtuple->Branch("idD1_reco", &idd1_reco, "idD1_reco[candSize]/I");
// 				VertexCompositeNtuple->Branch("idD2_reco", &idd2_reco, "idD2_reco[candSize]/I");
// 				VertexCompositeNtuple->Branch("matchGEN", &matchGEN, "matchGEN[candSize]/O");
// 				VertexCompositeNtuple->Branch("gen3DPointingAngle", &gen_agl_abs, "gen3DPointingAngle[candSize]/F");
// 				VertexCompositeNtuple->Branch("gen2DPointingAngle", &gen_agl2D_abs, "gen2DPointingAngle[candSize]/F");
// 				VertexCompositeNtuple->Branch("gen3DDecayLength", &gen_dl, "gen3DDecayLength[candSize]/F");
// 				VertexCompositeNtuple->Branch("gen2DDecayLength", &gen_dl2D, "gen2DDecayLength[candSize]/F");
// 			}
// 		}

// 	} // doRecoNtuple_

// 	if (doGenNtuple_)
// 	{
// 		VertexCompositeNtuple->Branch("pT_gen", &pt_gen, "pT_gen[candSize]/F");
// 		VertexCompositeNtuple->Branch("eta_gen", &eta_gen, "eta_gen[candSize]/F");
// 		VertexCompositeNtuple->Branch("y_gen", &y_gen, "y_gen[candSize]/F");
// 		VertexCompositeNtuple->Branch("phi_gen", &phi_gen, "phi_gen[candSize]/F");
// 		VertexCompositeNtuple->Branch("MotherID_gen", &idmom, "MotherID_gen[candSize]/I");

// 		if (decayInGen_)
// 		{

// 			VertexCompositeNtuple->Branch("DauID1_gen", &iddau1, "DauID1_gen[candSize]/I");
// 			VertexCompositeNtuple->Branch("DauID2_gen", &iddau2, "DauID2_gen[candSize]/I");
// 		}
// 	}
// }

void VCTreeProducer_D02kpi::initTree()
{
	VertexCompositeNtuple = fs->make<TTree>("VCNtuple_D02kpi", "VCNtuple_D02kpi");

	if (doRecoNtuple_)
	{

		// Event info
		VertexCompositeNtuple->Branch("Ntrkoffline", &Ntrkoffline, "Ntrkoffline/I");
		VertexCompositeNtuple->Branch("Npixel", &Npixel, "Npixel/I");
		VertexCompositeNtuple->Branch("HFsumETPlus", &HFsumETPlus, "HFsumETPlus/F");
		VertexCompositeNtuple->Branch("HFsumETMinus", &HFsumETMinus, "HFsumETMinus/F");
		VertexCompositeNtuple->Branch("ZDCPlus", &ZDCPlus, "ZDCPlus/F");
		VertexCompositeNtuple->Branch("ZDCMinus", &ZDCMinus, "ZDCMinus/F");
		VertexCompositeNtuple->Branch("PvtxX", &bestvx, "PvtxX/F");
		VertexCompositeNtuple->Branch("PvtxY", &bestvy, "PvtxY/F");
		VertexCompositeNtuple->Branch("PvtxZ", &bestvz, "PvtxZ/F");
		VertexCompositeNtuple->Branch("BSx", &BSx, "BSx/F");
		VertexCompositeNtuple->Branch("BSy", &BSy, "BSy/F");
		VertexCompositeNtuple->Branch("BSz", &BSz, "BSz/F");
		VertexCompositeNtuple->Branch("PvtxXErr", &bestvxError, "PvtxXErr/F");
		VertexCompositeNtuple->Branch("PvtxYErr", &bestvyError, "PvtxYErr/F");
		VertexCompositeNtuple->Branch("PvtxZErr", &bestvzError, "PvtxZErr/F");
		VertexCompositeNtuple->Branch("BSxErr", &BSxerror, "BSxErr/F");
		VertexCompositeNtuple->Branch("BSyErr", &BSyerror, "BSyErr/F");
		VertexCompositeNtuple->Branch("BSzErr", &BSzerror, "BSzErr/F");
		VertexCompositeNtuple->Branch("candSize", &candSize, "candSize/I");
		if (isCentrality_)
			VertexCompositeNtuple->Branch("centrality", &centrality, "centrality/I");
		// particle info
		VertexCompositeNtuple->Branch("pT", &pt);
		VertexCompositeNtuple->Branch("y", &y);
		VertexCompositeNtuple->Branch("phi", &phi);
		VertexCompositeNtuple->Branch("mass", &mass);
		if (useAnyMVA_)
		{
			VertexCompositeNtuple->Branch("mva", &mva);
			VertexCompositeNtuple->Branch("mva_xg", &mva_xg);
		}

		if (!isSkimMVA_)
		{
			// Composite candidate info RECO
			VertexCompositeNtuple->Branch("flavor", &flavor);
			VertexCompositeNtuple->Branch("eta", &eta);
			VertexCompositeNtuple->Branch("VtxProb", &VtxProb);
			VertexCompositeNtuple->Branch("VtxChi2", &vtxChi2);
			VertexCompositeNtuple->Branch("3DPointingAngle", &agl_abs);
			VertexCompositeNtuple->Branch("2DPointingAngle", &agl2D_abs);
			VertexCompositeNtuple->Branch("3DDecayLengthSignificance", &dlos);
			VertexCompositeNtuple->Branch("3DDecayLength", &dl);
			VertexCompositeNtuple->Branch("3DDecayLengthError", &dlerror);
			VertexCompositeNtuple->Branch("2DDecayLengthSignificance", &dlos2D);
			VertexCompositeNtuple->Branch("2DDecayLength", &dl2D);
			VertexCompositeNtuple->Branch("2DDecayLengthError", &dl2Derror);
			VertexCompositeNtuple->Branch("DlxyBS", &DlxyBS);
			VertexCompositeNtuple->Branch("DlxyBSErr", &DlxyBSErr);
			VertexCompositeNtuple->Branch("ip3d", &ip3d);
			VertexCompositeNtuple->Branch("ip3derr", &ip3derr);
			VertexCompositeNtuple->Branch("zDCASignificanceDaugther1", &dzos1);
			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1", &dxyos1);
			VertexCompositeNtuple->Branch("pTD1", &pt1);
			VertexCompositeNtuple->Branch("pTerrD1", &ptErr1);
			VertexCompositeNtuple->Branch("EtaD1", &eta1);
			VertexCompositeNtuple->Branch("PhiD1", &phi1);
			VertexCompositeNtuple->Branch("dedxHarmonic2D1", &H2dedx1);
			VertexCompositeNtuple->Branch("zDCASignificanceDaugther2", &dzos2);
			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2", &dxyos2);
			VertexCompositeNtuple->Branch("pTD2", &pt2);
			VertexCompositeNtuple->Branch("pTerrD2", &ptErr2);
			VertexCompositeNtuple->Branch("EtaD2", &eta2);
			VertexCompositeNtuple->Branch("PhiD2", &phi2);
			VertexCompositeNtuple->Branch("dedxHarmonic2D2", &H2dedx2);
			VertexCompositeNtuple->Branch("Dtrk1Dz1", &Dtrk1Dz1);
			VertexCompositeNtuple->Branch("Dtrk2Dz1", &Dtrk2Dz1);
			VertexCompositeNtuple->Branch("Dtrk1Dxy1", &Dtrk1Dxy1);
			VertexCompositeNtuple->Branch("Dtrk2Dxy1", &Dtrk2Dxy1);
			VertexCompositeNtuple->Branch("Dtrk1DzError1", &Dtrk1DzError1);
			VertexCompositeNtuple->Branch("Dtrk2DzError1", &Dtrk2DzError1);
			VertexCompositeNtuple->Branch("Dtrk1DxyError1", &Dtrk1DxyError1);
			VertexCompositeNtuple->Branch("Dtrk2DxyError1", &Dtrk2DxyError1);
			VertexCompositeNtuple->Branch("twoTrackDCA", &twoTrackDCA);

			if (doGenMatching_)
			{
				VertexCompositeNtuple->Branch("isSwap", &isSwap);
				VertexCompositeNtuple->Branch("idmom_reco", &idmom_reco);
				VertexCompositeNtuple->Branch("idD1_reco", &idd1_reco);
				VertexCompositeNtuple->Branch("idD2_reco", &idd2_reco);
				VertexCompositeNtuple->Branch("matchGEN", &matchGEN);
				VertexCompositeNtuple->Branch("gen3DPointingAngle", &gen_agl_abs);
				VertexCompositeNtuple->Branch("gen2DPointingAngle", &gen_agl2D_abs);
				VertexCompositeNtuple->Branch("gen3DDecayLength", &gen_dl);
				VertexCompositeNtuple->Branch("gen2DDecayLength", &gen_dl2D);
			}
		}

	} // doRecoNtuple_

	if (doGenNtuple_)
	{
		VertexCompositeNtuple->Branch("pT_gen", &pt_gen);
		VertexCompositeNtuple->Branch("eta_gen", &eta_gen);
		VertexCompositeNtuple->Branch("y_gen", &y_gen);
		VertexCompositeNtuple->Branch("phi_gen", &phi_gen);
		VertexCompositeNtuple->Branch("MotherID_gen", &idmom);

		if (decayInGen_)
		{

			VertexCompositeNtuple->Branch("DauID1_gen", &iddau1);
			VertexCompositeNtuple->Branch("DauID2_gen", &iddau2);
		}
	}
}

// ------------ method called once each job just after ending the event
// loop  ------------
void VCTreeProducer_D02kpi::endJob()
{
}

void VCTreeProducer_D02kpi::genDecayLength(const uint &it, const reco::GenParticle &gCand)
{
   	gen_dl[it] = -99.; gen_agl_abs[it] = -99.; gen_dl2D[it] = -99.; gen_agl2D_abs[it] = -99.;
	if(gCand.numberOfDaughters()==0 || !gCand.daughter(0)) return;
	const auto& dauVtx = gCand.daughter(0)->vertex();
	TVector3 ptosvec(dauVtx.X(), dauVtx.Y(), dauVtx.Z());
	TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
	gen_agl_abs[it] = secvec.Angle(ptosvec);
	gen_dl[it]  = ptosvec.Mag();
	TVector3 ptosvec2D(dauVtx.X(), dauVtx.Y(), 0.0);
	TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
	gen_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
	gen_dl2D[it]  = ptosvec2D.Mag();
}

// define this as a plug-in
DEFINE_FWK_MODULE(VCTreeProducer_D02kpi);
