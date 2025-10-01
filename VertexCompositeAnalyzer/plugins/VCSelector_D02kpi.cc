// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>

#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
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

#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include "BDT_header.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

//
// class decleration
//

#define PI 3.1416

const int y_bins = 2;
const int cent_bins = 4;
const int pT_bins = 10;

using namespace std;
using namespace cms::Ort;

using namespace std;

class BDTHandler
{
	private:
		double bdtCuts[y_bins][cent_bins][pT_bins];

		tuple<int, int, int> get_bins(float y, int centrality, float pT) const
		{
			int y_bin = -1, cent_bin = -1, pT_bin = -1;

			if (centrality >= 0 && centrality < 2 * 10)
				cent_bin = 0;
			else if (centrality >= 2 * 10 && centrality < 2 * 30)
				cent_bin = 1;
			else if (centrality >= 2 * 30 && centrality < 2 * 50)
				cent_bin = 2;
			else if (centrality >= 2 * 50 && centrality < 2 * 90)
				cent_bin = 3;

			if (pT >= 1 && pT < 2)
				pT_bin = 0;
			else if (pT >= 2 && pT < 3)
				pT_bin = 1;
			else if (pT >= 3 && pT < 4)
				pT_bin = 2;
			else if (pT >= 4 && pT < 5)
				pT_bin = 3;
			else if (pT >= 5 && pT < 6)
				pT_bin = 4;
			else if (pT >= 6 && pT < 8)
				pT_bin = 5;
			else if (pT >= 8 && pT < 10)
				pT_bin = 6;
			else if (pT >= 10 && pT < 15)
				pT_bin = 7;
			else if (pT >= 15 && pT < 20)
				pT_bin = 8;
			else if (pT >= 20)
				pT_bin = 9;

			if (abs(y) < 1)
				y_bin = 0;
			else if (abs(y) >= 1 && abs(y) < 3)
				y_bin = 1;

			return make_tuple(y_bin, cent_bin, pT_bin);
		}

	public:
		BDTHandler(bool useAnyMVA_, const edm::ParameterSet &iConfig)
		{
			if (useAnyMVA_) {
				std::string bdt_cuts_path = "";
				if(iConfig.exists("bdtCutsFile")){
					std::string cuts_filename = iConfig.getParameter<std::string>("bdtCutsFile");
					edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeAnalyzer/data/%s", cuts_filename.c_str()));
					bdt_cuts_path = fip.fullPath();
					std::ifstream testFile(bdt_cuts_path);
					if (!testFile.good())
					{
						throw cms::Exception("Configuration") << "cannot find BDT cuts in : " << bdt_cuts_path;
					}
					testFile.close();
				}

				loadCuts(bdt_cuts_path);
			}
		}

		void loadCuts(const string &cuts_filename)
		{
			ifstream file(cuts_filename);
			cout << "the following should be path to csv file" << endl; 
			cout << "cuts filename: " << cuts_filename << endl;
			if (!file.is_open())
				throw cms::Exception("Configuration") << "cannot find BDT cuts in : " << cuts_filename << endl;
			string line;

			while (getline(file, line))
			{
				if (line[0] == '#')
					continue;
				stringstream ss(line);
				int y_bin = -9, cent_bin = -9, pT_bin = -9;
				double cut_value;

				ss >> y_bin >> cent_bin >> pT_bin >> cut_value;

				bdtCuts[y_bin][cent_bin][pT_bin] = cut_value; // writes bdt cut values to the bdCuts array
			}
			file.close();
		}

		inline double getBDTCut(float y, int centrality, float pT) const
		{
			int y_bin, cent_bin, pT_bin;
			tie(y_bin, cent_bin, pT_bin) = get_bins(y, centrality, pT);
			if (y_bin < 0 || cent_bin < 0 || pT_bin < 0)
			{
				return 999;
			}
			return bdtCuts[y_bin][cent_bin][pT_bin];
		}
};

class VCSelector_D02kpi : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>>
{

	public:
		explicit VCSelector_D02kpi(const edm::ParameterSet &, const ONNXRuntime *cache);
		static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
		static void globalEndJob(const ONNXRuntime *);
		~VCSelector_D02kpi();

		using MVACollection = std::vector<float>;

	private:
		virtual void beginJob();
		virtual void produce(edm::Event &, const edm::EventSetup &);
		void fillRECO(edm::Event &iEvent, const edm::EventSetup &iSetup);
		virtual void endJob();

		std::vector<std::string> theInputVars;
		vector<double> inputValues;
		ReadBDT *mva;
		std::unique_ptr<BDTHandler> bdt;
		std::vector<std::string> input_names_;
		std::vector<std::string> output_names_;
		std::vector<std::vector<int64_t>> input_shapes_;
		std::string onnxModelPath_;
		const ONNXRuntime *onnxRuntime_;

		// ----------member data ---------------------------

		// options

		bool doGenMatching_;
		bool hasSwap_;
		bool decayInGen_;
		bool twoLayerDecay_;
		bool threeProngDecay_;
		bool doMuon_;
		bool selectGenMatch_;
		bool selectGenUnMatch_;
		bool selectGenMatchSwap_;
		bool selectGenMatchUnSwap_;

		int PID_;
		int PID_dau1_;
		int PID_dau2_;
		int PID_dau3_;

		// cut variables
		double deltaR_; // deltaR for Gen matching
		bool trkHighPurity_;
		double trkPMin_;
		double trkPtMin_;
		double trkEtaMax_;
		double trkPSumMin_;
		double trkPtSumMin_;
		double trkPtAsymMin_;
		double trkEtaDiffMax_;
		double trkPtErrMax_;
		int trkNHitMin_;
		double candpTMin_;
		double candpTMax_;
		double candYMin_;
		double candYMax_;
		double cand3DDecayLengthSigMin_;
		double cand2DDecayLengthSigMin_;
		double cand3DPointingAngleMax_;
		double cand2DPointingAngleMax_;
		double cand3DDCAMin_;
		double cand3DDCAMax_;
		double cand2DDCAMin_;
		double cand2DDCAMax_;
		double candVtxProbMin_;

		double mvaCut_;

		// tree branches
		// event info
		int centrality;
		int Ntrkoffline;
		float bestvx;
		float bestvy;
		float bestvz;

		// Composite candidate info
		float mva_value;
		float onnxVal;
		float mva_old;
		float pt;
		float eta;
		float flavor;
		float y;
		float mass;
		float VtxProb;
		float dlos;
		float dl;
		float dlerror;
		float agl;
		float vtxChi2;
		float ndf;
		float agl_abs;
		float agl2D;
		float agl2D_abs;
		float dlos2D;
		bool isSwap;
		bool matchGEN;
		int idmom_reco;

		// dau info
		float dzos1;
		float dzos2;
		float dzos3;
		float dxyos1;
		float dxyos2;
		float dxyos3;

		float nhit1;
		float nhit2;
		bool trkquality1;
		bool trkquality2;
		float pt1;
		float pt2;
		float ptErr1;
		float ptErr2;
		float p1;
		float p2;
		float eta1;
		float eta2;
		float phi1;
		float phi2;
		float H2dedx1;
		float H2dedx2;
		float T4dedx1;
		float T4dedx2;
		bool isPionD1;
		bool isPionD2;
		bool isKaonD1;
		bool isKaonD2;

		int selectFlavor_;
		bool usePID_;
		bool useAnyMVA_;
		bool applyXGB_;

		bool assignBDT = true;

		std::vector<float> mvaVals_;

		TF2 *func_mva;
		std::vector<double> mvaCuts_;

		TH2D *hist_bdtcut;

		float mvaMin_;
		float mvaMax_;

		int centMin_;
		int centMax_;
		bool isCentrality_;

		edm::Handle<int> cbin_;

		// tokens
		edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tok_generalTrk_;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;
		edm::EDGetTokenT<MVACollection> MVAValues_Token_;
		edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token1_;
		edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> Dedx_Token2_;
		edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
		edm::EDGetTokenT<int> tok_centBinLabel_;
		edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

		std::string d0IDName_;

		pat::CompositeCandidateCollection theGoodCandidates;
		MVACollection theMVANew;
		MVACollection theMVANew_xg;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

VCSelector_D02kpi::VCSelector_D02kpi(const edm::ParameterSet &iConfig, const ONNXRuntime *cache)
	:  input_shapes_(), onnxRuntime_(cache)
{
	string a1 = "log3ddls";
	string a2 = "nVtxProb";
	string a3 = "n3DPointingAngle";
	string a4 = "nDtrk1Pt";
	string a5 = "nDtrk2Pt";
	string a6 = "nxyDCASigD1";
	string a7 = "nxyDCASigD2";
	string a8 = "nzDCASigD1";
	string a9 = "nzDCASigD2";
	string a10 = "npT";
	string a11 = "ny";
	string a12 = "ncent";

	theInputVars.push_back(a1);
	theInputVars.push_back(a2);
	theInputVars.push_back(a3);
	theInputVars.push_back(a4);
	theInputVars.push_back(a5);
	theInputVars.push_back(a6);
	theInputVars.push_back(a7);
	theInputVars.push_back(a8);
	theInputVars.push_back(a9);
	theInputVars.push_back(a10);
	theInputVars.push_back(a12);
	theInputVars.push_back(a11);
	mva = new ReadBDT(theInputVars);
	// options

	threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");

	PID_ = iConfig.getUntrackedParameter<int>("PID");
	PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
	PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
	if (threeProngDecay_)
		PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");

	// cut variables
	centMin_ = iConfig.getUntrackedParameter<int>("centMin", 0);
	centMax_ = iConfig.getUntrackedParameter<int>("centMax", 10000);
	deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);
	mvaMax_ = iConfig.getUntrackedParameter<double>("mvaMax", 999.9);
	mvaMin_ = iConfig.getUntrackedParameter<double>("mvaMin", -999.9);

	trkHighPurity_ = iConfig.getUntrackedParameter<bool>("trkHighPurity");
	trkPMin_ = iConfig.getUntrackedParameter<double>("trkPMin");
	trkPtMin_ = iConfig.getUntrackedParameter<double>("trkPtMin");
	trkEtaMax_ = iConfig.getUntrackedParameter<double>("trkEtaMax");
	trkPSumMin_ = iConfig.getUntrackedParameter<double>("trkPSumMin");
	trkPtSumMin_ = iConfig.getUntrackedParameter<double>("trkPtSumMin");
	trkPtAsymMin_ = iConfig.getUntrackedParameter<double>("trkPtAsymMin");
	trkEtaDiffMax_ = iConfig.getUntrackedParameter<double>("trkEtaDiffMax");
	trkPtErrMax_ = iConfig.getUntrackedParameter<double>("trkPtErrMax");
	trkNHitMin_ = iConfig.getUntrackedParameter<int>("trkNHitMin");
	candpTMin_ = iConfig.getUntrackedParameter<double>("candpTMin");
	candpTMax_ = iConfig.getUntrackedParameter<double>("candpTMax");
	candYMin_ = iConfig.getUntrackedParameter<double>("candYMin");
	candYMax_ = iConfig.getUntrackedParameter<double>("candYMax");
	cand3DDecayLengthSigMin_ = iConfig.getUntrackedParameter<double>("cand3DDecayLengthSigMin");
	cand2DDecayLengthSigMin_ = iConfig.getUntrackedParameter<double>("cand2DDecayLengthSigMin");
	cand3DPointingAngleMax_ = iConfig.getUntrackedParameter<double>("cand3DPointingAngleMax");
	cand2DPointingAngleMax_ = iConfig.getUntrackedParameter<double>("cand2DPointingAngleMax");
	cand3DDCAMin_ = iConfig.getUntrackedParameter<double>("cand3DDCAMin");
	cand3DDCAMax_ = iConfig.getUntrackedParameter<double>("cand3DDCAMax");
	cand2DDCAMin_ = iConfig.getUntrackedParameter<double>("cand2DDCAMin");
	cand2DDCAMax_ = iConfig.getUntrackedParameter<double>("cand2DDCAMax");
	candVtxProbMin_ = iConfig.getUntrackedParameter<double>("candVtxProbMin");
	mvaCut_ = iConfig.getParameter<double>(string("mvaCut"));

	// input tokens

	tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));

	patCompositeCandidateCollection_Token_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("VertexCompositeCollection"));

	tok_generalTrk_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));
	Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxHarmonic2"));
	Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxTruncated40"));
	tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

	usePID_ = false;
	selectFlavor_ = 0;
	if (iConfig.exists("usePID"))
		usePID_ = iConfig.getParameter<bool>("usePID");
	if (iConfig.exists("useFlavor"))
		selectFlavor_ = iConfig.getUntrackedParameter<int>("selectFlavor");

	// Loading TMVA
	useAnyMVA_ = false;

	isCentrality_ = false;
	if (iConfig.exists("isCentrality"))
		isCentrality_ = iConfig.getParameter<bool>("isCentrality");
	if (isCentrality_)
	{
		tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
		tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
	}

	if (iConfig.exists("applyXGB"))
		applyXGB_ = iConfig.getParameter<bool>("applyXGB");
	if (iConfig.exists("useAnyMVA"))
		useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");

	if (useAnyMVA_)
	{
		if (iConfig.exists("input_names") || iConfig.exists("output_names"))
		{
			input_names_ = iConfig.getParameter<std::vector<std::string>>("input_names");
			output_names_ = iConfig.getParameter<std::vector<std::string>>("output_names");
		}
		else
		{
			throw cms::Exception("Configuration") << "onnxModelName not provided in ParameterSet";
		}
	}
	bdt = std::make_unique<BDTHandler>(useAnyMVA_, iConfig);


	d0IDName_ = (iConfig.getParameter<edm::InputTag>("VertexCompositeCollection")).instance();

	produces<pat::CompositeCandidateCollection>(d0IDName_);
	produces<MVACollection>(Form("MVAValuesNew%s", d0IDName_.c_str()));

	produces<MVACollection>(Form("MVAValuesNew%s2", d0IDName_.c_str()));

	isPionD1 = true;
	isPionD2 = true;
	isKaonD1 = false;
	isKaonD2 = false;
}
std::unique_ptr<ONNXRuntime> VCSelector_D02kpi::initializeGlobalCache(const edm::ParameterSet &iConfig)
{
	bool useAnyMVA = iConfig.exists("useAnyMVA") ? iConfig.getParameter<bool>("useAnyMVA") : false;

	if (!useAnyMVA)
		return nullptr;

	if (iConfig.exists("onnxModelFileName"))
	{
		std::string onnxModelPath = iConfig.getParameter<std::string>("onnxModelFileName");

		edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeProducer/data/%s", onnxModelPath.c_str()));
		std::string fullPath = fip.fullPath();
		std::cout << fullPath << std::endl;

		std::ifstream testFile(fullPath);
		if (!testFile.good())
		{
			throw cms::Exception("Configuration") << "cannot find ONNX Model in : " << fullPath;
		}
		testFile.close();

		return std::make_unique<ONNXRuntime>(fip.fullPath());
	}

	return nullptr;
}
void VCSelector_D02kpi::globalEndJob(const ONNXRuntime *cache) {}

VCSelector_D02kpi::~VCSelector_D02kpi()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to for each event  ------------
void VCSelector_D02kpi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{

	using std::vector;
	using namespace edm;

	edm::Handle<pat::CompositeCandidateCollection> patCandidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_, patCandidates);

	if (!patCandidates.isValid())
	{
		edm::LogError("VCSelector_D02kpi") << "Error: patCandidates collection not found!";
		return;
	}

	fillRECO(iEvent, iSetup);

	auto theNewD0Cands = std::make_unique<pat::CompositeCandidateCollection>();
	theNewD0Cands->reserve(theGoodCandidates.size());

	std::copy(theGoodCandidates.begin(),
			theGoodCandidates.end(),
			std::back_inserter(*theNewD0Cands));

	// Store final EDM output
	iEvent.put(std::move(theNewD0Cands), d0IDName_);

	theGoodCandidates.clear();

	if (useAnyMVA_)
	{

		auto mvas = std::make_unique<MVACollection>(theMVANew.begin(), theMVANew.end());
		auto mvas_xg = std::make_unique<MVACollection>(theMVANew_xg.begin(), theMVANew_xg.end());

		iEvent.put(std::move(mvas), Form("MVAValuesNew%s", d0IDName_.c_str()));
		iEvent.put(std::move(mvas_xg), Form("MVAValuesNew%s2", d0IDName_.c_str()));

		theMVANew.clear();
		theMVANew_xg.clear();
	}
}

void VCSelector_D02kpi::fillRECO(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
	// get collections
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(tok_offlinePV_, vertices);

	edm::Handle<pat::PackedCandidateCollection> tracks;
	iEvent.getByToken(tok_generalTrk_, tracks);

	edm::Handle<pat::CompositeCandidateCollection> d0candidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_, d0candidates);
	const pat::CompositeCandidateCollection *d0candidates_ = d0candidates.product();

	edm::Handle<MVACollection> mvavalues;

	edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle1;
	if (usePID_)
		iEvent.getByToken(Dedx_Token1_, dEdxHandle1);

	edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle2;
	if (usePID_)
		iEvent.getByToken(Dedx_Token2_, dEdxHandle2);

	centrality = -1;
	if (isCentrality_)
	{
		edm::Handle<reco::Centrality> cent;
		iEvent.getByToken(tok_centSrc_, cent);

		iEvent.getByToken(tok_centBinLabel_, cbin_);
		centrality = *cbin_;
	}
	if (centrality != -1 && (centrality >= centMax_ || centrality < centMin_))
		return;

	// best vertex
	bestvz = -999.9;
	bestvx = -999.9;
	bestvy = -999.9;
	double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
	const reco::Vertex &vtx = (*vertices)[0];
	bestvz = vtx.z();
	bestvx = vtx.x();
	bestvy = vtx.y();
	bestvzError = vtx.zError();
	bestvxError = vtx.xError();
	bestvyError = vtx.yError();

	// Ntrkoffline
	Ntrkoffline = 0;

	// RECO Candidate info
	for (unsigned it = 0; it < d0candidates_->size(); ++it)
	{

		const pat::CompositeCandidate &trk = (*d0candidates_)[it];
		bestvz = trk.userFloat("bestvtxZ");
		bestvx = trk.userFloat("bestvtxX");
		bestvy = trk.userFloat("bestvtxY"); 
		bestvzError = trk.userFloat("zVtxError");
		bestvxError = trk.userFloat("xVtxError");
		bestvyError = trk.userFloat("yVtxError");

		double bdt_cut_value = -999.9;
		double secvz = -999.9, secvx = -999.9, secvy = -999.9;
		secvz = trk.userFloat("vtxZ");
		secvx = trk.userFloat("vtxX");
		secvy = trk.userFloat("vtxY");

		eta = trk.eta();
		y = trk.rapidity();
		pt = trk.pt();
		flavor = trk.pdgId() / 421;

		double px = trk.px();
		double py = trk.py();
		double pz = trk.pz();
		mass = trk.mass();


		auto const *d1 = dynamic_cast<const pat::PackedCandidate *>(trk.daughter(0));
		auto const *d2 = dynamic_cast<const pat::PackedCandidate *>(trk.daughter(1));

		auto pseudoTrk1 = d1->pseudoTrack(); // today change
		auto pseudoTrk2 = d2->pseudoTrack(); // today change

		if (!d1)
		{
			edm::LogError("VertexCompositeTreeProducer") << "Failed to cast daughter 0 to pat::PackedCandidate";
		}

		if (!d2)
		{
			edm::LogError("VertexCompositeTreeProducer") << "Failed to cast daughter 1 to pat::PackedCandidate";
		}

		// Gen match
		//  select particle vs antiparticle
		if (usePID_ && selectFlavor_ && (int)flavor != selectFlavor_)
			continue;
		if (pt < candpTMin_ || pt > candpTMax_)
			continue;
		if (abs(y) < candYMin_ || abs(y) > candYMax_)
			continue;

		pt1 = d1->pt();
		pt2 = d2->pt();

		if (pt1 < trkPtMin_ || pt2 < trkPtMin_)
			continue;
		if ((pt1 + pt2) < trkPtSumMin_)
			continue;

		if (pt2 / pt1 < trkPtAsymMin_ || pt1 / pt2 < trkPtAsymMin_)
			continue;

		// momentum
		p1 = d1->p();
		p2 = d2->p();

		if (p1 < trkPMin_ || p2 < trkPMin_)
			continue;
		if ((p1 + p2) < trkPSumMin_)
			continue;

		// eta
		eta1 = d1->eta();
		eta2 = d2->eta();

		if (fabs(eta1) > trkEtaMax_ || fabs(eta2) > trkEtaMax_)
			continue;
		if (fabs(eta1 - eta2) > trkEtaDiffMax_)
			continue;

		// phi
		phi1 = d1->phi();
		phi2 = d2->phi();

		vtxChi2 = trk.userFloat("vtxChi2");
		ndf = trk.userFloat("vtxNdof");
		VtxProb = TMath::Prob(vtxChi2, ndf);

		if (VtxProb < candVtxProbMin_)
			continue;

		// PAngle
		TVector3 ptosvec(secvx - bestvx, secvy - bestvy, secvz - bestvz);
		TVector3 secvec(px, py, pz);

		TVector3 ptosvec2D(secvx - bestvx, secvy - bestvy, 0);
		TVector3 secvec2D(px, py, 0);

		agl = cos(secvec.Angle(ptosvec));
		agl_abs = secvec.Angle(ptosvec);
		if (agl_abs > cand3DPointingAngleMax_)
			continue;

		agl2D = cos(secvec2D.Angle(ptosvec2D));
		agl2D_abs = secvec2D.Angle(ptosvec2D);
		if (agl2D_abs > cand2DPointingAngleMax_)
			continue;

		// Decay length 3D
		typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
		typedef ROOT::Math::SVector<double, 3> SVector3;
		typedef ROOT::Math::SVector<double, 6> SVector6;

		reco::Vertex::CovarianceMatrix sec_covariance;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				sec_covariance(i, j) = trk.userFloat("vertexCovariance_" + std::to_string(i) + "_" + std::to_string(j));
			}
		}

		SMatrixSym3D totalCov = vtx.covariance() + sec_covariance;
		SVector3 distanceVector(secvx - bestvx, secvy - bestvy, secvz - bestvz);

		dl = ROOT::Math::Mag(distanceVector);
		dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / dl;

		dlos = dl / dlerror;

		if (dlos < cand3DDecayLengthSigMin_ || dlos > 1000.)
			continue;

		// Decay length 2D
		SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1), vtx.covariance(1, 1), 0, 0, 0);
		SVector6 v2(sec_covariance(0, 0), sec_covariance(0, 1), sec_covariance(1, 1), 0, 0, 0);

		SMatrixSym3D sv1(v1);
		SMatrixSym3D sv2(v2);

		SMatrixSym3D totalCov2D = sv1 + sv2;
		SVector3 distanceVector2D(secvx - bestvx, secvy - bestvy, 0);

		double dl2D = ROOT::Math::Mag(distanceVector2D);
		double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D)) / dl2D;

		dlos2D = dl2D / dl2Derror;
		if (dlos2D < cand2DDecayLengthSigMin_ || dlos2D > 1000.)
			continue;

		double dca3D = dl * sin(agl_abs);
		if (dca3D < cand3DDCAMin_ || dca3D > cand3DDCAMax_)
			continue;

		double dca2D = dl2D * sin(agl2D_abs);
		if (dca2D < cand2DDCAMin_ || dca2D > cand2DDCAMax_)
			continue;

		// trk info

		float ptErr1 = pseudoTrk1.ptError(); // need this  for XGBoost
		if (ptErr1 / pt1 > trkPtErrMax_)
			continue;
		int nhit1 = pseudoTrk1.hitPattern().numberOfValidHits(); // cuts
		if (nhit1 < trkNHitMin_) continue;

		// trk dEdx
		H2dedx1 = -999.9;
		T4dedx1 = -999.9;

		H2dedx1 = -999.9;
		T4dedx1 = -999.9;
		math::XYZPoint bestvtx(bestvx, bestvy, bestvz);


		double dzbest1 = d1->pseudoTrack().dz(bestvtx); // today change
		double dxybest1 = d1->pseudoTrack().dxy(bestvtx); // today change
		double dzerror1 = TMath::Sqrt(d1->pseudoTrack().dzError() * d1->pseudoTrack().dzError() + bestvzError * bestvzError);
		double dxyerror1 = TMath::Sqrt(d1->pseudoTrack().dxyError()*d1->pseudoTrack().dxyError() + bestvxError*bestvyError);

		dzos1 = dzbest1 / dzerror1;
		dxyos1 = dxybest1 / dxyerror1;

		float ptErr2 = pseudoTrk2.ptError();
		if (ptErr2 / pt2 > trkPtErrMax_)
			continue;
		int nhit2 = pseudoTrk2.hitPattern().numberOfValidHits();
		if (nhit2 < trkNHitMin_)
			continue;

		// trk dEdx
		H2dedx2 = -999.9;
		T4dedx2 = -999.9;


		double dzbest2 = d2->pseudoTrack().dz(bestvtx);    // today change
		double dxybest2 = d2->pseudoTrack().dxy(bestvtx); // today change
		double dzerror2 = TMath::Sqrt(d2->pseudoTrack().dzError() * d2->pseudoTrack().dzError() + bestvzError * bestvzError);
		double dxyerror2 = TMath::Sqrt(d2->pseudoTrack().dxyError()*d2->pseudoTrack().dxyError() + bestvxError*bestvyError);

		dzos2 = dzbest2 / dzerror2;
		dxyos2 = dxybest2 / dxyerror2;

		mva_value = -999.9;
		onnxVal = -999.9;
		// if (useAnyMVA_ && onnxRuntime_)
		if (useAnyMVA_ )
		{

			if (applyXGB_ && pt > 2 && abs(y) < 1.5){
				cms::Ort::FloatArrays data_;
				data_.emplace_back(19, 0);
				std::vector<float> &onnxVals_ = data_[0];
				onnxVals_[0] = pt;
				onnxVals_[1] = y;
				onnxVals_[2] = VtxProb;
				onnxVals_[3] = centrality;
				onnxVals_[4] = agl;
				onnxVals_[5] = agl_abs;
				onnxVals_[6] = agl2D;
				onnxVals_[7] = agl2D_abs;
				onnxVals_[8] = dl;
				onnxVals_[9] = dlos;
				onnxVals_[10] = dl2D;
				onnxVals_[11] = dlos2D;
				onnxVals_[12] = pt1;
				onnxVals_[13] = eta1;
				onnxVals_[14] = pt2;
				onnxVals_[15] = eta2;
				onnxVals_[16] = ptErr1;
				onnxVals_[17] = ptErr2;
				onnxVals_[18] = trk.userFloat("track3DDCA");

				std::vector<float> outputs = onnxRuntime_->run(input_names_, data_, input_shapes_, output_names_)[0];
				onnxVal = outputs[1];
			}

			inputValues.clear();
			inputValues.push_back(log10(dlos)); // 00
			inputValues.push_back(VtxProb);     // 01
			inputValues.push_back(agl_abs);     // 02
			inputValues.push_back(pt1);         // 03
			inputValues.push_back(pt2);         // 04
			inputValues.push_back(dxyos1);      // 04
			inputValues.push_back(dxyos2);      // 04
			inputValues.push_back(dzos1);       // 04
			inputValues.push_back(dzos2);
			inputValues.push_back(pt);
			inputValues.push_back(centrality);
			inputValues.push_back(y);
			mva_value = mva->GetMvaValue(inputValues);
			bdt_cut_value = bdt->getBDTCut(y, centrality, pt);


			if (mva_value <= bdt_cut_value && onnxVal <= mvaCut_)
			//if (mva_value <= bdt_cut_value )
		        //if (bdt_cut_value < -2)
			continue;

			// cout << "candidate passed the cuts no problem" << endl;
			theMVANew.push_back(mva_value);
			theMVANew_xg.push_back(onnxVal);
		}

		// select MVA value
		theGoodCandidates.push_back(trk);
	}
}

// ------------ method called once each job just before starting event
// loop  ------------
void VCSelector_D02kpi::beginJob()
{
}

// ------------ method called once each job just after ending the event
// loop  ------------
void VCSelector_D02kpi::endJob()
{
}

// define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(VCSelector_D02kpi);
