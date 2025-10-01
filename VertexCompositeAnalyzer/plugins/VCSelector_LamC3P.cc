//+++++++++++++++++++++++++++++++++++                                                                                                            

//Code owner: Nihar Ranjan saha                                                                                                                  
//contact mail: nihar.ranjan.saha@cern.ch                                                                                                        
//Date: 24 Sept 2025                                                                                                                             

/*Note: This code, called selector,where
  we apply all necessary selections/cuts to
  Lc candidates and daughters coming from
  Fitter code before saving into TTree.
  */
//+++++++++++++++++++++++++++++++++++


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

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
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#define MAXCAN 20000

//
// class decleration
//

#define PI 3.1416

using namespace std;


class VCSelector_LamC3P : public edm::one::EDProducer<> {
	public:
		explicit VCSelector_LamC3P(const edm::ParameterSet&);
		~VCSelector_LamC3P();

		using MVACollection = std::vector<float>;

	private:
		virtual void beginJob() ;
		virtual void produce(edm::Event&, const edm::EventSetup&);
		virtual void fillRECO(edm::Event&, const edm::EventSetup&) ;
		virtual void endJob() ;

		// ----------member data ---------------------------

		//options
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

		//cut variables
		double multMax_;
		double multMin_;
		double deltaR_; //deltaR for Gen matching
		bool   trkHighPurity_;
		double trkPMin_;
		double trkPtMin_;
		double trkEtaMax_;
		double trkPSumMin_;
		double trkPtSumMin_;
		double trkPtAsymMin_;
		double trkEtaDiffMax_;
		double trkPtErrMax_;
		int    trkNHitMin_;
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

		//tree branches
		//event info
		int centrality;
		int Ntrkoffline;
		float bestvx;
		float bestvy;
		float bestvz;

		//Composite candidate info
		float mva;
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

		//dau candidate info
		float grand_mass;
		float grand_VtxProb;
		float grand_dlos;
		float grand_dl;
		float grand_dlerror;
		float grand_agl;
		float grand_vtxChi2;
		float grand_ndf;
		float grand_agl_abs;
		float grand_agl2D;
		float grand_agl2D_abs;
		float grand_dlos2D;


		//dau info
		float dzos1;
		float dzos2;
		float dzos3;
		float dxyos1;
		float dxyos2;
		float dxyos3;
		float nhit1;
		float nhit2;
		float nhit3;
		bool trkquality1;
		bool trkquality2;
		bool trkquality3;
		float pt1;
		float pt2;
		float pt3;
		float ptErr1;
		float ptErr2;
		float ptErr3;
		float p1;
		float p2;
		float p3;
		float eta1;
		float eta2;
		float eta3;
		float phi1;
		float phi2;
		float phi3;
		float H2dedx1;
		float H2dedx2;
		float H2dedx3;
		float T4dedx1;
		float T4dedx2;
		float T4dedx3;
		float trkChi1;
		float trkChi2;
		float trkChi3;
		bool  isPionD1;
		bool  isPionD2;
		bool  isPionD3;
		bool  isKaonD1;
		bool  isKaonD2;
		bool  isKaonD3;

		//grand-dau info
		float grand_dzos1;
		float grand_dzos2;
		float grand_dxyos1;
		float grand_dxyos2;
		float grand_nhit1;
		float grand_nhit2;
		bool grand_trkquality1;
		bool grand_trkquality2;
		float grand_pt1;
		float grand_pt2;
		float grand_ptErr1;
		float grand_ptErr2;
		float grand_p1;
		float grand_p2;
		float grand_eta1;
		float grand_eta2;
		float grand_H2dedx1;
		float grand_H2dedx2;
		float grand_T4dedx1;
		float grand_T4dedx2;
		float grand_trkChi1;
		float grand_trkChi2;

		//dau muon info
		float nmatchedst1;
		float nmatchedch1;
		float ntrackerlayer1;
		float npixellayer1;
		float matchedenergy1;
		float nmatchedst2;
		float nmatchedch2;
		float ntrackerlayer2;
		float npixellayer2;
		float matchedenergy2;
		float dx1_seg_;
		float dy1_seg_;
		float dxSig1_seg_;
		float dySig1_seg_;
		float ddxdz1_seg_;
		float ddydz1_seg_;
		float ddxdzSig1_seg_;
		float ddydzSig1_seg_;
		float dx2_seg_;
		float dy2_seg_;
		float dxSig2_seg_;
		float dySig2_seg_;
		float ddxdz2_seg_;
		float ddydz2_seg_;
		float ddxdzSig2_seg_;
		float ddydzSig2_seg_;

		//vector for gen match
		vector< vector<double> > *pVect;
		vector<double> *Dvector1;
		vector<double> *Dvector2;
		vector<double> *Dvector3;
		vector<int> *pVectIDmom;

		int  selectFlavor_;
		bool usePID_;
		bool useAnyMVA_;
		bool useExistingMVA_;

		std::string mvaType_;
		std::string forestLabel_;
		GBRForest * forest_;
		bool useForestFromDB_;
		std::vector<float> mvaVals_;
		std::string dbFileName_;

		TF2* func_mva;
		std::vector<double> mvaCuts_;

		TH2D* hist_bdtcut;

		float mvaMin_;
		float mvaMax_;

		int   centMin_;
		int   centMax_;
		bool isCentrality_;

		edm::Handle<int> cbin_;

		//tokens
		edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;

		edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;
		edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
		edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
		edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
		edm::EDGetTokenT<reco::MuonCollection> tok_muon_;
		edm::EDGetTokenT<int> tok_centBinLabel_;
		edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

		std::string CandIDName_;

		pat::CompositeCandidateCollection theVertexComps;
};



VCSelector_LamC3P::VCSelector_LamC3P(const edm::ParameterSet& iConfig)
{
	//options
	twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
	threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
	hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
	decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
	doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
	selectGenMatch_ = iConfig.getUntrackedParameter<bool>("selectGenMatch");
	selectGenUnMatch_ = iConfig.getUntrackedParameter<bool>("selectGenUnMatch");
	selectGenMatchSwap_ = iConfig.getUntrackedParameter<bool>("selectGenMatchSwap");
	selectGenMatchUnSwap_ = iConfig.getUntrackedParameter<bool>("selectGenMatchUnSwap");

	PID_ = iConfig.getUntrackedParameter<int>("PID");
	PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
	PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
	if(threeProngDecay_) PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");

	//cut variables
	centMin_ = iConfig.getUntrackedParameter<int>("centMin");
	centMax_ = iConfig.getUntrackedParameter<int>("centMax");
	multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
	multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
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

	//input tokens
	tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));

	patCompositeCandidateCollection_Token_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("VertexCompositeCollection"));
	tok_muon_ = consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection"));
	Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
	Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));

	usePID_ = false;
	selectFlavor_ = 0;
	if(iConfig.exists("usePID")) usePID_ = iConfig.getParameter<bool>("usePID");
	if(iConfig.exists("useFlavor")) selectFlavor_ = iConfig.getUntrackedParameter<int>("selectFlavor");

	// Loading TMVA
	useExistingMVA_ = false;

	forestLabel_ = "DStarInPbPb";
	std::string type = "BDT";
	useForestFromDB_ = false;
	dbFileName_ = "";

	forest_ = nullptr;

	isCentrality_ = iConfig.getParameter<bool>("isCentrality");

	if(isCentrality_)
	{
		tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
		tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
	}


	CandIDName_ = (iConfig.getParameter<edm::InputTag>("VertexCompositeCollection")).instance();


	produces< pat::CompositeCandidateCollection >(CandIDName_);

}


VCSelector_LamC3P::~VCSelector_LamC3P()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
	void
VCSelector_LamC3P::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using std::vector;
	using namespace edm;

	fillRECO(iEvent,iSetup);

	edm::Handle<pat::CompositeCandidateCollection> patCandidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_, patCandidates);

	if (!patCandidates.isValid())
	{
		edm::LogError("VCSelector_LamC3P") << "Error: patCandidates collection not found!";
		return;
	}


	auto theNewLcCands = std::make_unique<pat::CompositeCandidateCollection>();

	theNewLcCands->reserve( theVertexComps.size() );

	std::copy( theVertexComps.begin(),
			theVertexComps.end(),
			std::back_inserter(*theNewLcCands) );

	iEvent.put(std::move(theNewLcCands), CandIDName_);
	theVertexComps.clear();


}

	void
VCSelector_LamC3P::fillRECO(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//get collections
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(tok_offlinePV_,vertices);

	edm::Handle<pat::CompositeCandidateCollection> LamC3PCandidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_,LamC3PCandidates);
	const pat::CompositeCandidateCollection * LamC3PCandidates_ = LamC3PCandidates.product();

	edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
	if(usePID_) iEvent.getByToken(Dedx_Token1_, dEdxHandle1);

	edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
	if(usePID_) iEvent.getByToken(Dedx_Token2_, dEdxHandle2);

	centrality=-1;
	if(isCentrality_)
	{
		edm::Handle<reco::Centrality> cent;
		iEvent.getByToken(tok_centSrc_, cent);

		iEvent.getByToken(tok_centBinLabel_,cbin_);
		centrality = *cbin_;

	}
	if(centrality!=-1 && (centrality >= centMax_ || centrality < centMin_)) return;

	//best vertex
	bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
	double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
	const reco::Vertex & vtx = (*vertices)[0];
	bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
	bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();


	//RECO Candidate info

	int candSize = LamC3PCandidates_->size();


	for(int it=0; it<candSize; ++it){

		const pat::CompositeCandidate & trk = (*LamC3PCandidates_)[it];


		bestvz = trk.userFloat("bestvtxZ");
		bestvx = trk.userFloat("bestvtxX");
		bestvy = trk.userFloat("bestvtxY"); 
		bestvzError = trk.userFloat("zVtxError");
		bestvxError = trk.userFloat("xVtxError");
		bestvyError = trk.userFloat("yVtxError");

		double secvz = -999.9, secvx = -999.9, secvy = -999.9;
		secvz = trk.userFloat("vtxZ");
		secvx = trk.userFloat("vtxX");
		secvy = trk.userFloat("vtxY");



		eta = trk.eta();
		y = trk.rapidity();
		pt = trk.pt();
		flavor = trk.pdgId()/4122;
		double px = trk.px();
		double py = trk.py();
		double pz = trk.pz();
		mass = trk.mass();



		if(pt<candpTMin_ || pt>candpTMax_) continue;
		if(y<candYMin_ || y>candYMax_) continue;

		vtxChi2 = trk.userFloat("vtxChi2");
		ndf = trk.userFloat("vtxNdof");
		VtxProb = TMath::Prob(vtxChi2, ndf);
		if(VtxProb < candVtxProbMin_ || VtxProb > 1.0) continue;

		const pat::PackedCandidate *d1_hypo = dynamic_cast<const pat::PackedCandidate *>(trk.daughter(0));
		const pat::PackedCandidate *d2_hypo = dynamic_cast<const pat::PackedCandidate *>(trk.daughter(1));
		const pat::PackedCandidate *d3_hypo = nullptr;
		if (threeProngDecay_) {d3_hypo = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(2));}



		if (!d1_hypo){edm::LogError("VertexCompositeTreeProducer") << "Failed to cast daughter 0 to pat::PackedCandidate";}
		if (!d2_hypo){edm::LogError("VertexCompositeTreeProducer") << "Failed to cast daughter 1 to pat::PackedCandidate";}
		if (!d3_hypo){edm::LogError("VertexCompositeTreeProducer") << "Failed to cast daughter 2 to pat::PackedCandidate";}



		//pt
		pt1 = d1_hypo->pt();
		p1 = d1_hypo->p();
		eta1 = d1_hypo->eta();
		phi1 = d1_hypo->phi();

		pt2 = d2_hypo->pt();
		p2 = d2_hypo->p();
		eta2 = d2_hypo->eta();
		phi2 = d2_hypo->phi();



		if(pt1 < trkPtMin_ || pt2 < trkPtMin_) continue;
		if(p1 < trkPMin_ || p2 < trkPMin_) continue;
		if(fabs(eta1) > trkEtaMax_ || fabs(eta2) > trkEtaMax_) continue;

		if((pt1+pt2) < trkPtSumMin_) continue;
		if(pt2/pt1 < trkPtAsymMin_ || pt1/pt2 < trkPtAsymMin_) continue;
		if((p1+p2) < trkPSumMin_) continue;
		if(fabs(eta1-eta2) > trkEtaDiffMax_) continue;



		if(threeProngDecay_)
		{
			pt3 = d3_hypo->pt();
			p3 = d3_hypo->p();
			eta3 = d3_hypo->eta();
			phi3 = d3_hypo->phi();

			if(pt3 < trkPtMin_) continue;
			if(p3 < trkPMin_) continue;
			if(fabs(eta3) > trkEtaMax_) continue;

		}



		//PAngle
		TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
		TVector3 secvec(px,py,pz);

		TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
		TVector3 secvec2D(px,py,0);

		agl = cos(secvec.Angle(ptosvec));
		agl_abs = secvec.Angle(ptosvec);
		if(agl_abs > cand3DPointingAngleMax_) continue;



		agl2D = cos(secvec2D.Angle(ptosvec2D));
		agl2D_abs = secvec2D.Angle(ptosvec2D);
		if(agl2D_abs > cand2DPointingAngleMax_) continue;

		//Decay length 3D
		typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
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

		if (dlos < cand3DDecayLengthSigMin_ || dlos > 1000.) continue;


		//Decay length 2D
		SVector6 v1(vtx.covariance(0, 0), vtx.covariance(0, 1), vtx.covariance(1, 1), 0, 0, 0);
		SVector6 v2(sec_covariance(0, 0), sec_covariance(0, 1), sec_covariance(1, 1), 0, 0, 0);


		SMatrixSym3D sv1(v1);
		SMatrixSym3D sv2(v2);

		SMatrixSym3D totalCov2D = sv1 + sv2;
		SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);

		double dl2D = ROOT::Math::Mag(distanceVector2D);
		double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D;

		dlos2D = dl2D/dl2Derror;
		if(dlos2D < cand2DDecayLengthSigMin_ || dlos2D > 1000.) continue;

		double dca3D = dl*sin(agl_abs);
		if(dca3D < cand3DDCAMin_ || dca3D > cand3DDCAMax_) continue;


		double dca2D = dl2D*sin(agl2D_abs);
		if(dca2D < cand2DDCAMin_ || dca2D > cand2DDCAMax_) continue;



		const reco::Track& pseudoTrk1 = d1_hypo->pseudoTrack();

		float ptErr1 = pseudoTrk1.ptError();
		int nhit1 = pseudoTrk1.hitPattern().numberOfValidHits();


		if(ptErr1/d1_hypo->pt() > trkPtErrMax_) continue;
		if(nhit1 < trkNHitMin_) continue;

		math::XYZPoint bestvtx(bestvx, bestvy, bestvz);	
		double dzbest1 = pseudoTrk1.dz(bestvtx);
		double dxybest1 = pseudoTrk1.dxy(bestvtx);
		double dzerror1 = std::sqrt(pseudoTrk1.dzError() * pseudoTrk1.dzError() + bestvzError * bestvzError);
		double dxyerror1 = std::sqrt(pseudoTrk1.dxyError() * pseudoTrk1.dxyError() + bestvxError * bestvyError);

		dzos1 = dzbest1/dzerror1;
		dxyos1 = dxybest1/dxyerror1;



		//Get Daughter-2 info!
		const reco::Track& pseudoTrk2 = d2_hypo->pseudoTrack();
		float ptErr2 = pseudoTrk2.ptError();
		int nhit2 = pseudoTrk2.hitPattern().numberOfValidHits();

		if(ptErr2/d2_hypo->pt() > trkPtErrMax_) continue;
		if(nhit2 < trkNHitMin_) continue;

		double dzbest2 = pseudoTrk2.dz(bestvtx);
		double dxybest2 = pseudoTrk2.dxy(bestvtx);
		double dzerror2 = std::sqrt(pseudoTrk2.dzError() * pseudoTrk2.dzError() + bestvzError * bestvzError);
		double dxyerror2 = std::sqrt(pseudoTrk2.dxyError() * pseudoTrk2.dxyError() + bestvxError * bestvyError);

		dzos2 = dzbest2/dzerror2;
		dxyos2 = dxybest2/dxyerror2;


		if(threeProngDecay_)
		{
			const reco::Track& pseudoTrk3 = d3_hypo->pseudoTrack();
			float ptErr3 = pseudoTrk3.ptError();
			int nhit3 = pseudoTrk3.hitPattern().numberOfValidHits();

			if(ptErr3/d3_hypo->pt() > trkPtErrMax_) continue;
			if(nhit3 < trkNHitMin_) continue;

			double dzbest3 = pseudoTrk3.dz(bestvtx);
			double dxybest3 = pseudoTrk3.dxy(bestvtx);
			double dzerror3 = std::sqrt(pseudoTrk3.dzError() * pseudoTrk3.dzError() + bestvzError * bestvzError);
			double dxyerror3 = std::sqrt(pseudoTrk3.dxyError() * pseudoTrk3.dxyError() + bestvxError * bestvyError);

			dzos3 = dzbest3/dzerror3;
			dxyos3 = dxybest3/dxyerror3;


		}


		theVertexComps.push_back( trk );


	}
}




// ------------ method called once each job just before starting event

	void
VCSelector_LamC3P::beginJob()
{
}

// ------------ method called once each job just after ending the event

void 
VCSelector_LamC3P::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(VCSelector_LamC3P);
