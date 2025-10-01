//+++++++++++++++++++++++++++++++++++                                                                                                            

//Code owner: Nihar Ranjan saha                                                                                                                  
//contact mail: nihar.ranjan.saha@cern.ch                                                                                                        
//Date: 24 Sept 2025                                                                                                                             

/*
Note: This code, called TTree producer, where
we fill all the required variables into output
TTree. You can add more branches here according
to your analysis. Here, all the candidates are
already passed selector conditions.
*/ 
//+++++++++++++++++++++++++++++++++++

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
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>
#include <vector>


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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//#include "Utils.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <iomanip>


#define PI 3.1416
#define MAXCAN 20000

using namespace std;


struct CandidateData {
    // Candidate info
    std::vector<float> mva;
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
    std::vector<float> Ddca;
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
    std::vector<int> idd3_reco;

    // Gen-level info (parent)
    std::vector<float> gen_agl_abs;
    std::vector<float> gen_agl2D_abs;
    std::vector<float> gen_dl;
    std::vector<float> gen_dl2D;

    // Daughter info
    std::vector<float> dzos1;
    std::vector<float> dzos2;
    std::vector<float> dzos3;
    std::vector<float> dxyos1;
    std::vector<float> dxyos2;
    std::vector<float> dxyos3;
    std::vector<float> pt1;
    std::vector<float> pt2;
    std::vector<float> pt3;
    std::vector<float> ptErr1;
    std::vector<float> ptErr2;
    std::vector<float> ptErr3;
    std::vector<float> p1;
    std::vector<float> p2;
    std::vector<float> p3;
    std::vector<float> Dtrk1Dz1;
    std::vector<float> Dtrk2Dz1;
    std::vector<float> Dtrk3Dz1;
    std::vector<float> Dtrk1Dxy1;
    std::vector<float> Dtrk2Dxy1;
    std::vector<float> Dtrk3Dxy1;
    std::vector<float> Dtrk1DzError1;
    std::vector<float> Dtrk2DzError1;
    std::vector<float> Dtrk3DzError1;
    std::vector<float> Dtrk1DxyError1;
    std::vector<float> Dtrk2DxyError1;
    std::vector<float> Dtrk3DxyError1;
    std::vector<float> eta1;
    std::vector<float> eta2;
    std::vector<float> eta3;
    std::vector<float> phi1;
    std::vector<float> phi2;
    std::vector<float> phi3;
    std::vector<int> charge1;
    std::vector<int> charge2;
    std::vector<int> charge3;
    std::vector<int> pid1;
    std::vector<int> pid2;
    std::vector<int> pid3;
    std::vector<float> tof1;
    std::vector<float> tof2;
    std::vector<float> tof3;
    std::vector<float> H2dedx1;
    std::vector<float> H2dedx2;
    std::vector<float> H2dedx3;
    std::vector<float> T4dedx1;
    std::vector<float> T4dedx2;
    std::vector<float> T4dedx3;
    std::vector<float> trkChi1;
    std::vector<float> trkChi2;
    std::vector<float> trkChi3;

    // Gen-level info (daughters)
    std::vector<float> pt_gen;
    std::vector<float> eta_gen;
    std::vector<int> idmom;
    std::vector<float> y_gen;
    std::vector<float> phi_gen;
    std::vector<int> iddau1;
    std::vector<int> iddau2;
    std::vector<int> iddau3;
};



class VCTreeProducer_LamC3P : public edm::one::EDAnalyzer<> {
	public:
		explicit VCTreeProducer_LamC3P(const edm::ParameterSet&);
		~VCTreeProducer_LamC3P();

		using MVACollection = std::vector<float>;

	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void fillRECO(const edm::Event&, const edm::EventSetup&) ;
		virtual void endJob() ;
		virtual void initTree();

		void genDecayLength(const uint&, const reco::GenParticle&, CandidateData&);

		// ----------member data ---------------------------

		edm::Service<TFileService> fs;


		TTree* VertexCompositeNtuple;

		CandidateData candInfo; 
		bool saveTree_;

		//options
		bool doRecoNtuple_;
		bool dogenntuple_;   

		bool dogenmatchingtof_;
		bool hasswap_;
		bool decayingen_;
		bool threeProngDecay_;
		int PID_;
		int PID_dau1_;
		int PID_dau2_;
		int PID_dau3_;

		//cut variables
		double multMax_;
		double multMin_;
		double deltaR_; //deltaR for Gen matching


		//tree branches
		//event info
		int   centMin_;
		int centMax_;
		int centrality;
		int Ntrkoffline;
		int Npixel;
		float HFsumET;
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





		void manageVectorsSize(CandidateData& data, size_t newSize) {

			data.mva.resize(newSize);
			data.pt.resize(newSize);
			data.eta.resize(newSize);
			data.phi.resize(newSize);
			data.flavor.resize(newSize);
			data.y.resize(newSize);
			data.mass.resize(newSize);
			data.VtxProb.resize(newSize);
			data.dlos.resize(newSize);
			data.dl.resize(newSize);
			data.dlerror.resize(newSize);
			data.DlxyBS.resize(newSize);
			data.DlxyBSErr.resize(newSize);
			data.vtxChi2.resize(newSize);
			data.ndf.resize(newSize);
			data.agl_abs.resize(newSize);
			data.Ddca.resize(newSize);
			data.ip3d.resize(newSize);
			data.ip3derr.resize(newSize);
			data.agl2D_abs.resize(newSize);
			data.dlos2D.resize(newSize);
			data.dl2D.resize(newSize);
			data.dl2Derror.resize(newSize);
			data.isSwap.resize(newSize);
			data.matchGEN.resize(newSize);
			data.idmom_reco.resize(newSize);
			data.idd1_reco.resize(newSize);
			data.idd2_reco.resize(newSize);
			data.idd3_reco.resize(newSize);

			data.gen_agl_abs.resize(newSize);
			data.gen_agl2D_abs.resize(newSize);
			data.gen_dl.resize(newSize);
			data.gen_dl2D.resize(newSize);

			data.dzos1.resize(newSize);
			data.dzos2.resize(newSize);
			data.dzos3.resize(newSize);
			data.dxyos1.resize(newSize);
			data.dxyos2.resize(newSize);
			data.dxyos3.resize(newSize);
			data.pt1.resize(newSize);
			data.pt2.resize(newSize);
			data.pt3.resize(newSize);
			data.ptErr1.resize(newSize);
			data.ptErr2.resize(newSize);
			data.ptErr3.resize(newSize);
			data.p1.resize(newSize);
			data.p2.resize(newSize);
			data.p3.resize(newSize);
			data.Dtrk1Dz1.resize(newSize);
			data.Dtrk2Dz1.resize(newSize);
			data.Dtrk3Dz1.resize(newSize);
			data.Dtrk1Dxy1.resize(newSize);
			data.Dtrk2Dxy1.resize(newSize);
			data.Dtrk3Dxy1.resize(newSize);
			data.Dtrk1DzError1.resize(newSize);
			data.Dtrk2DzError1.resize(newSize);
			data.Dtrk3DzError1.resize(newSize);
			data.Dtrk1DxyError1.resize(newSize);
			data.Dtrk2DxyError1.resize(newSize);
			data.Dtrk3DxyError1.resize(newSize);
			data.eta1.resize(newSize);
			data.eta2.resize(newSize);
			data.eta3.resize(newSize);
			data.phi1.resize(newSize);
			data.phi2.resize(newSize);
			data.phi3.resize(newSize);
			data.charge1.resize(newSize);
			data.charge2.resize(newSize);
			data.charge3.resize(newSize);
			data.pid1.resize(newSize);
			data.pid2.resize(newSize);
			data.pid3.resize(newSize);
			data.tof1.resize(newSize);
			data.tof2.resize(newSize);
			data.tof3.resize(newSize);
			data.H2dedx1.resize(newSize);
			data.H2dedx2.resize(newSize);
			data.H2dedx3.resize(newSize);
			data.T4dedx1.resize(newSize);
			data.T4dedx2.resize(newSize);
			data.T4dedx3.resize(newSize);
			data.trkChi1.resize(newSize);
			data.trkChi2.resize(newSize);
			data.trkChi3.resize(newSize);

			data.pt_gen.resize(newSize);
			data.eta_gen.resize(newSize);
			data.idmom.resize(newSize);
			data.y_gen.resize(newSize);
			data.phi_gen.resize(newSize);
			data.iddau1.resize(newSize);
			data.iddau2.resize(newSize);
			data.iddau3.resize(newSize);
		}


		vector< vector<double> > *pVect;
		vector<double> *Dvector1;
		vector<double> *Dvector2;
		vector<double> *Dvector3;
		vector<int> *pVectIDmom;



		bool isSkimMVA_;
		bool isCentrality_;
		bool doGenNtuple_;
		bool useAnyMVA_;
		bool doGenMatching_;
		bool decayInGen_;

		edm::Handle<int> cbin_;

		//tokens
		edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tok_generalTrk_;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;

		edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;

		edm::EDGetTokenT<int> tok_centBinLabel_;
		edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

		edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

		edm::EDGetTokenT< reco::BeamSpot > bsLabel_;

};


VCTreeProducer_LamC3P::VCTreeProducer_LamC3P(const edm::ParameterSet& iConfig)

{
	//options
	doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
	doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
	doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
	decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
	PID_ = iConfig.getUntrackedParameter<int>("PID");
	PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
	PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
	PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");

	saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
	threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
	isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA"); 
	isCentrality_ = iConfig.getParameter<bool>("isCentrality");

	centMin_ = iConfig.getUntrackedParameter<int>("centMin");
	centMax_ = iConfig.getUntrackedParameter<int>("centMax");

	multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
	multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
	deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.02);


	//input tokens
	tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
	tok_generalTrk_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));
	patCompositeCandidateCollection_Token_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LamC3P"));
	tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getParameter<edm::InputTag>("GenParticleCollection")));
	bsLabel_        = consumes< reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("BSLabel"));

	if(isCentrality_)
	{
		tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
		tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
	}

}


VCTreeProducer_LamC3P::~VCTreeProducer_LamC3P()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VCTreeProducer_LamC3P::analyze(const edm::Event& iEvent, const edm::EventSetup&
		iSetup)
{
	using std::vector;
	using namespace edm;
	using namespace reco;

	if(doRecoNtuple_) fillRECO(iEvent,iSetup);

	if(saveTree_) VertexCompositeNtuple->Fill();
}

void VCTreeProducer_LamC3P::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


#ifdef DEBUG
	using std::cout;
	using std::endl;
#endif
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(tok_offlinePV_,vertices);

	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(bsLabel_, beamSpotHandle);

	edm::Handle<pat::CompositeCandidateCollection> lamC3Pcandidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_,lamC3Pcandidates);

	const pat::CompositeCandidateCollection * lamC3Pcandidates_ = lamC3Pcandidates.product();

	edm::Handle<reco::GenParticleCollection> genpars;	
	if (doGenMatching_) {
		iEvent.getByToken(tok_genParticle_, genpars);
	}

#ifdef DEBUG
	cout << "Loaded tokens" << endl;
#endif

	centrality=-1;
	if(isCentrality_)
	{
		edm::Handle<reco::Centrality> cent;
		iEvent.getByToken(tok_centSrc_, cent);
		HFsumET = (cent.isValid() ? cent->EtHFtowerSum() : -1.); //Nihar
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



	//best vertex
	BSx=-999.9; BSy=-999.9; BSz=-999.9;
	BSxerror=-999.9; BSyerror=-999.9; BSzerror=-999.9;

	float BSdxdz = -999.9;
	float BSdydz = -999.9;
	bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
	bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;

	const reco::Vertex & vtx = (*vertices)[0];
	bestvx = vtx.x();
	bestvy = vtx.y();
	bestvz = vtx.z();

	bestvxError = vtx.xError();
	bestvyError = vtx.yError();
	bestvzError = vtx.zError();

	float xlxyBS = -999.9;
	float ylxyBS = -999.9;
	float vtxYXErr = -999.9;
	float vtxXErr = -999.9;
	float vtxYErr = -999.9;


	reco::BeamSpot beamSpot = *beamSpotHandle;
	reco::Vertex theBeamSpotV(beamSpot.position(), beamSpot.covariance3D());

	BSx      = beamSpot.x0();
	BSy      = beamSpot.y0();
	BSz      = beamSpot.z0();
	BSxerror = beamSpot.x0Error();
	BSyerror = beamSpot.y0Error();
	BSzerror = beamSpot.z0Error();
	BSdxdz   = beamSpot.dxdz();
	BSdydz   = beamSpot.dydz();


	//RECO Candidate info
	candSize = lamC3Pcandidates_->size();	
	manageVectorsSize(candInfo, candSize);


	for(int it=0; it<candSize; ++it){


		const pat::CompositeCandidate & trk = (*lamC3Pcandidates_)[it];


		double secvz=-999.9, secvx=-999.9, secvy=-999.9;
		secvz = trk.userFloat("vtxZ"); secvx = trk.userFloat("vtxX"); secvy = trk.userFloat("vtxY");
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


		candInfo.eta[it] = trk.eta();
		candInfo.y[it] = trk.rapidity();
		candInfo.pt[it] = trk.pt();
		candInfo.phi[it] = trk.phi();
		candInfo.flavor[it] = trk.pdgId()/abs(trk.pdgId());

		double px = trk.px();
		double py = trk.py();
		double pz = trk.pz();
		candInfo.mass[it] = trk.mass();

		if (candInfo.mass[it] == 0) {
			cout << "Error break" << endl;
		}


		const pat::PackedCandidate* reco_d1 = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(0));
		const pat::PackedCandidate* reco_d2 = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(1));
		const pat::PackedCandidate* reco_d3 = nullptr;

		if (threeProngDecay_) {
			reco_d3 = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(2)); 
		}

		//Note: This is very important to assign the charge based on pdgId!!!
		//Note: Without this we can't get correct charge assignment based on permutation!!
		int reco_d1_charge = TMath::Sign(1, reco_d1->pdgId());
		int reco_d2_charge = TMath::Sign(1, reco_d2->pdgId());
		int reco_d3_charge = TMath::Sign(1, reco_d3->pdgId());




		//Gen-matching!!	  
		candInfo.matchGEN[it] = false;
		candInfo.isSwap[it] = false;
		candInfo.idmom_reco[it] = -77;
		candInfo.idd1_reco[it] = -77;
		candInfo.idd2_reco[it] = -77;
		candInfo.idd3_reco[it] = -77;

		candInfo.pt_gen[it] = -999.9;
		candInfo.eta_gen[it] = -999.9;
		candInfo.idmom[it] = -999;
		candInfo.y_gen[it] = -999.9;
		candInfo.phi_gen[it] = -999.9;
		candInfo.iddau1[it] = -999;
		candInfo.iddau2[it] = -999;
		candInfo.iddau3[it] = -999;


		if(doGenMatching_){

			if(!genpars.isValid()){
				cout<<"Gen matching cannot be done without Gen collection!!"<<endl;
				return;
			}

			for(unsigned genPair=0; genPair<genpars->size(); ++genPair){

				const reco::GenParticle & genLamC = (*genpars)[genPair];

				if(fabs(genLamC.pdgId()) != PID_) continue;


				const reco::Candidate * gen_d1 = genLamC.daughter(0);
				const reco::Candidate * gen_d2 = genLamC.daughter(1);
				const reco::Candidate * gen_d3 = genLamC.daughter(2);


				bool passPID = (
						(fabs(gen_d1->pdgId())==PID_dau1_ && fabs(gen_d2->pdgId())==PID_dau2_ && fabs(gen_d3->pdgId())==PID_dau3_) ||
						(fabs(gen_d1->pdgId())==PID_dau1_ && fabs(gen_d2->pdgId())==PID_dau3_ && fabs(gen_d3->pdgId())==PID_dau2_) ||
						(fabs(gen_d1->pdgId())==PID_dau2_ && fabs(gen_d2->pdgId())==PID_dau1_ && fabs(gen_d3->pdgId())==PID_dau3_) ||
						(fabs(gen_d1->pdgId())==PID_dau2_ && fabs(gen_d2->pdgId())==PID_dau3_ && fabs(gen_d3->pdgId())==PID_dau1_) ||
						(fabs(gen_d1->pdgId())==PID_dau3_ && fabs(gen_d2->pdgId())==PID_dau1_ && fabs(gen_d3->pdgId())==PID_dau2_) ||
						(fabs(gen_d1->pdgId())==PID_dau3_ && fabs(gen_d2->pdgId())==PID_dau2_ && fabs(gen_d3->pdgId())==PID_dau1_)
					       );
				if(!passPID) continue;


				// To check all 6 charge permutations
				struct CandidateCombo { 
					const reco::Candidate* g1; 
					const reco::Candidate* g2; 
					const reco::Candidate* g3; 
				};


				std::vector<CandidateCombo> gen_combos = {
					{gen_d1, gen_d2, gen_d3},
					{gen_d1, gen_d3, gen_d2},
					{gen_d2, gen_d1, gen_d3},
					{gen_d2, gen_d3, gen_d1},
					{gen_d3, gen_d1, gen_d2},
					{gen_d3, gen_d2, gen_d1}
				};


				for(auto &c : gen_combos) {


					if(reco_d1_charge!=c.g1->charge() || reco_d2_charge!=c.g2->charge() || reco_d3_charge!=c.g3->charge()) continue; 

					double dR1 = reco::deltaR(*reco_d1, *c.g1);
					if(dR1 >= deltaR_) continue;
					if(fabs((reco_d1->pt() - c.g1->pt()) / reco_d1->pt()) >=0.2) continue;
					double dR2 = reco::deltaR(*reco_d2, *c.g2);
					if(dR2 >= deltaR_) continue;
					if(fabs((reco_d2->pt() - c.g2->pt()) / reco_d2->pt()) >= 0.2) continue;
					double dR3 = reco::deltaR(*reco_d3, *c.g3);
					if(dR3 >= deltaR_) continue;
					if(fabs((reco_d3->pt() - c.g3->pt()) / reco_d3->pt()) >= 0.2) continue;


					candInfo.matchGEN[it] = true;


					if(reco_d1->pdgId() != c.g1->pdgId() ||
							reco_d2->pdgId() != c.g2->pdgId() ||
							reco_d3->pdgId() != c.g3->pdgId()) {
						candInfo.isSwap[it] = true;
					}

					if(decayInGen_){
						candInfo.iddau1[it] = c.g1->pdgId();
						candInfo.iddau2[it] = c.g2->pdgId();
						candInfo.iddau3[it] = c.g3->pdgId();
					}


					genDecayLength(it, genLamC, candInfo);

					candInfo.pt_gen[it] = genLamC.pt();
					candInfo.eta_gen[it] = genLamC.eta();
					candInfo.y_gen[it] = genLamC.rapidity();
					candInfo.phi_gen[it] = genLamC.phi();
					candInfo.idmom[it] = genLamC.pdgId();

					candInfo.idmom_reco[it] = trk.pdgId();                                                                
					candInfo.idd1_reco[it] = reco_d1->pdgId();
					candInfo.idd2_reco[it] = reco_d2->pdgId();
					candInfo.idd3_reco[it] = reco_d3->pdgId();

					break;
				}// End gen combo loop

				if(candInfo.matchGEN[it]) break;


			}// End genPair

		}//End doGenMatching!



		double pxd1 = reco_d1->px();
		double pyd1 = reco_d1->py();
		double pzd1 = reco_d1->pz();
		double pxd2 = reco_d2->px();
		double pyd2 = reco_d2->py();
		double pzd2 = reco_d2->pz();


		TVector3 dauvec1(pxd1,pyd1,pzd1);
		TVector3 dauvec2(pxd2,pyd2,pzd2);


		//pt
		candInfo.pt1[it] = reco_d1->pt();
		candInfo.pt2[it] = reco_d2->pt();

		//momentum
		candInfo.p1[it] = reco_d1->p();
		candInfo.p2[it] = reco_d2->p();

		//eta
		candInfo.eta1[it] = reco_d1->eta();
		candInfo.eta2[it] = reco_d2->eta();

		//phi
		candInfo.phi1[it] = reco_d1->phi();
		candInfo.phi2[it] = reco_d2->phi();

		//charge
		candInfo.charge1[it] = reco_d1->charge();
		candInfo.charge2[it] = reco_d2->charge();

		double pxd3 = -999.9;
		double pyd3 = -999.9;
		double pzd3 = -999.9;


		if(threeProngDecay_ )
		{

			pxd3 = reco_d3->px();
			pyd3 = reco_d3->py();
			pzd3 = reco_d3->pz();
			candInfo.pt3[it] = reco_d3->pt();
			candInfo.p3[it] = reco_d3->p();
			candInfo.eta3[it] = reco_d3->eta();
			candInfo.phi3[it] = reco_d3->phi();
			candInfo.charge3[it] = reco_d3->charge();
			TVector3 dauvec3(pxd3,pyd3,pzd3);
		}



		candInfo.pid1[it] = -99999;
		candInfo.pid2[it] = -99999;
		candInfo.pid3[it] = -99999;


		//vtxChi2
		candInfo.vtxChi2[it] = trk.userFloat("vtxChi2");
		candInfo.ndf[it] = trk.userFloat("vtxNdof");
		candInfo.VtxProb[it] = TMath::Prob(candInfo.vtxChi2[it],candInfo.ndf[it]);


		//PAngle
		TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
		TVector3 secvec(px,py,pz);

		TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
		TVector3 secvec2D(px,py,0);


		candInfo.agl_abs[it] = secvec.Angle(ptosvec);
		candInfo.Ddca[it] = ptosvec.Mag() * TMath::Sin(candInfo.agl_abs[it]);
		candInfo.ip3d[it] = trk.userFloat("ip3d");
		candInfo.ip3derr[it] = trk.userFloat("ip3derr");

		candInfo.agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

		float r2lxyBS = (secvx-BSx-(secvz-BSz)*BSdxdz) * (secvx-BSx-(secvz-BSz)*BSdxdz) + (secvy-BSy-(secvz-BSz)*BSdydz) * (secvy-BSy-(secvz-BSz)*BSdydz);
		xlxyBS = secvx-BSx - (secvz-BSz)*BSdxdz;
		ylxyBS = secvy-BSy - (secvz-BSz)*BSdydz;
		candInfo.DlxyBS[it] = static_cast<float>(TMath::Sqrt(r2lxyBS));
		candInfo.DlxyBSErr[it] = static_cast<float>(TMath::Sqrt ((1./r2lxyBS) * ((xlxyBS*xlxyBS)*vtxXErr + (2*xlxyBS*ylxyBS)*vtxYXErr + (ylxyBS*ylxyBS)*vtxYErr)));


		//Decay length 3D
		typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
		typedef ROOT::Math::SVector<double, 3> SVector3;
		typedef ROOT::Math::SVector<double, 6> SVector6;


		SMatrixSym3D totalCov = vtx.covariance() + sec_covariance;
		SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

		candInfo.dl[it] = ROOT::Math::Mag(distanceVector);
		candInfo.dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/candInfo.dl[it];

		candInfo.dlos[it] = candInfo.dl[it]/candInfo.dlerror[it];


		//Decay length 2D

		SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
		SVector6 v2(sec_covariance(0,0), sec_covariance(0,1),sec_covariance(1,1),0,0,0);

		SMatrixSym3D sv1(v1);
		SMatrixSym3D sv2(v2);

		SMatrixSym3D totalCov2D = sv1 + sv2;
		SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);

		candInfo.dl2D[it] = ROOT::Math::Mag(distanceVector2D);
		candInfo.dl2Derror[it] = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/candInfo.dl2D[it];

		candInfo.dlos2D[it] = candInfo.dl2D[it]/candInfo.dl2Derror[it];


		const reco::Track& pseudoTrk1 = reco_d1->pseudoTrack();


		candInfo.trkChi1[it] = pseudoTrk1.normalizedChi2();
		candInfo.ptErr1[it] = pseudoTrk1.ptError();

		//DCA
		math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
		math::XYZPoint BS_vtx(BSx,BSy,BSz);


		double dzbest1 = pseudoTrk1.dz(bestvtx);	  
		double dxybest1 = pseudoTrk1.dxy(bestvtx); 
		double dzerror1 = std::sqrt(pseudoTrk1.dzError() * pseudoTrk1.dzError() + bestvzError * bestvzError);
		double dxyerror1 = std::sqrt(pseudoTrk1.dxyError() * pseudoTrk1.dxyError() + bestvxError * bestvyError);

		candInfo.dzos1[it] = dzbest1/dzerror1;
		candInfo.dxyos1[it] = dxybest1/dxyerror1;
		candInfo.Dtrk1Dz1[it] = 1.0*dzbest1;
		candInfo.Dtrk1Dxy1[it] = dxybest1;
		candInfo.Dtrk1DzError1[it] = dzerror1;
		candInfo.Dtrk1DxyError1[it] = dxyerror1;


		const reco::Track& pseudoTrk2 = reco_d2->pseudoTrack();
		candInfo.trkChi2[it] = pseudoTrk2.normalizedChi2();
		candInfo.ptErr2[it] = pseudoTrk2.ptError();

		double dzbest2 = pseudoTrk2.dz(bestvtx);	  
		double dxybest2 = pseudoTrk2.dxy(bestvtx); 
		double dzerror2 = TMath::Sqrt(pseudoTrk2.dzError() * pseudoTrk2.dzError() + bestvzError * bestvzError);
		double dxyerror2 = TMath::Sqrt(pseudoTrk2.dxyError() * pseudoTrk2.dxyError() + bestvxError * bestvyError);

		candInfo.dzos2[it] = dzbest2/dzerror2;
		candInfo.dxyos2[it] = dxybest2/dxyerror2;
		candInfo.Dtrk2Dz1[it] = 1.0*dzbest2;
		candInfo.Dtrk2Dxy1[it] = dxybest2;
		candInfo.Dtrk2DzError1[it] = dzerror2;
		candInfo.Dtrk2DxyError1[it] = dxyerror2;



		if (threeProngDecay_){

			const reco::Track& pseudoTrk3 = reco_d3->pseudoTrack();
			candInfo.trkChi3[it] = pseudoTrk3.normalizedChi2();
			candInfo.ptErr3[it] = pseudoTrk3.ptError();

			double dzbest3 = pseudoTrk3.dz(bestvtx);	  
			double dxybest3 = pseudoTrk3.dxy(bestvtx); 
			double dzerror3 = TMath::Sqrt(pseudoTrk3.dzError() * pseudoTrk3.dzError() + bestvzError * bestvzError);
			double dxyerror3 = TMath::Sqrt(pseudoTrk3.dxyError() * pseudoTrk3.dxyError() + bestvxError * bestvyError);

			candInfo.dzos3[it] = dzbest3/dzerror3;
			candInfo.dxyos3[it] = dxybest3/dxyerror3;
			candInfo.Dtrk3Dz1[it] = 1.0*dzbest3;
			candInfo.Dtrk3Dxy1[it] = dxybest3;
			candInfo.Dtrk3DzError1[it] = dzerror3;
			candInfo.Dtrk3DxyError1[it] = dxyerror3;
		}



#ifdef DEBUG
		cout << "Done reco single iter" << endl;
#endif

	}// Candidate loop

#ifdef DEBUG
	cout << "Fill reco done" << endl;
#endif
}




// ------------ method called once each job just before starting event loop  ------------
	void
VCTreeProducer_LamC3P::beginJob()
{
	TH1D::SetDefaultSumw2();

	if(!doRecoNtuple_ && !doGenNtuple_)
	{
		cout<<"No output for either RECO or GEN!! Fix config!!"<<endl; return;
	}


	if(saveTree_) initTree();
}

	void 
VCTreeProducer_LamC3P::initTree()
{ 
	VertexCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");

	if(doRecoNtuple_) 
	{ 

		// Event info
		VertexCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
		VertexCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
		VertexCompositeNtuple->Branch("HFsumET",&HFsumET,"HFsumET/F");
		VertexCompositeNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
		VertexCompositeNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
		VertexCompositeNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
		VertexCompositeNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
		VertexCompositeNtuple->Branch("PvtxX",&bestvx,"PvtxX/F");
		VertexCompositeNtuple->Branch("PvtxY",&bestvy,"PvtxY/F");
		VertexCompositeNtuple->Branch("PvtxZ",&bestvz,"PvtxZ/F");
		VertexCompositeNtuple->Branch("BSx",&BSx,"BSx/F");
		VertexCompositeNtuple->Branch("BSy",&BSy,"BSy/F");
		VertexCompositeNtuple->Branch("BSz",&BSz,"BSz/F");
		VertexCompositeNtuple->Branch("PvtxXErr",&bestvxError,"PvtxXErr/F");
		VertexCompositeNtuple->Branch("PvtxYErr",&bestvyError,"PvtxYErr/F");
		VertexCompositeNtuple->Branch("PvtxZErr",&bestvzError,"PvtxZErr/F");
		VertexCompositeNtuple->Branch("BSxErr",&BSxerror,"BSxErr/F");
		VertexCompositeNtuple->Branch("BSyErr",&BSyerror,"BSyErr/F");
		VertexCompositeNtuple->Branch("BSzErr",&BSzerror,"BSzErr/F");
		VertexCompositeNtuple->Branch("candSize",&candSize,"candSize/I");
		if(isCentrality_) VertexCompositeNtuple->Branch("centrality",&centrality,"centrality/I");
		VertexCompositeNtuple->Branch("pT", &candInfo.pt);
		VertexCompositeNtuple->Branch("y", &candInfo.y);
		VertexCompositeNtuple->Branch("phi", &candInfo.phi);
		VertexCompositeNtuple->Branch("mass", &candInfo.mass);

		if(!isSkimMVA_)  
		{
			//Composite candidate info RECO
			VertexCompositeNtuple->Branch("flavor", &candInfo.flavor);
			VertexCompositeNtuple->Branch("eta", &candInfo.eta);
			VertexCompositeNtuple->Branch("VtxProb", &candInfo.VtxProb);
			VertexCompositeNtuple->Branch("VtxChi2", &candInfo.vtxChi2);
			VertexCompositeNtuple->Branch("3DPointingAngle", &candInfo.agl_abs);
			VertexCompositeNtuple->Branch("Ddca", &candInfo.Ddca);
			VertexCompositeNtuple->Branch("ip3d", &candInfo.ip3d);
			VertexCompositeNtuple->Branch("ip3derr", &candInfo.ip3derr);
			VertexCompositeNtuple->Branch("2DPointingAngle", &candInfo.agl2D_abs);
			VertexCompositeNtuple->Branch("3DDecayLengthSignificance", &candInfo.dlos);
			VertexCompositeNtuple->Branch("3DDecayLength", &candInfo.dl);
			VertexCompositeNtuple->Branch("3DDecayLengthError", &candInfo.dlerror);
			VertexCompositeNtuple->Branch("2DDecayLengthSignificance", &candInfo.dlos2D);
			VertexCompositeNtuple->Branch("2DDecayLength", &candInfo.dl2D);
			VertexCompositeNtuple->Branch("2DDecayLengthError", &candInfo.dl2Derror);
			VertexCompositeNtuple->Branch("DlxyBS", &candInfo.DlxyBS);
			VertexCompositeNtuple->Branch("DlxyBSErr", &candInfo.DlxyBSErr);
			VertexCompositeNtuple->Branch("zDCASignificanceDaugther1", &candInfo.dzos1);
			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1", &candInfo.dxyos1);
			VertexCompositeNtuple->Branch("pTD1", &candInfo.pt1);
			VertexCompositeNtuple->Branch("pTerrD1", &candInfo.ptErr1);
			VertexCompositeNtuple->Branch("EtaD1", &candInfo.eta1);
			VertexCompositeNtuple->Branch("PhiD1", &candInfo.phi1);
			VertexCompositeNtuple->Branch("dedxHarmonic2D1", &candInfo.H2dedx1);
			VertexCompositeNtuple->Branch("zDCASignificanceDaugther2", &candInfo.dzos2);
			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2", &candInfo.dxyos2);
			VertexCompositeNtuple->Branch("pTD2", &candInfo.pt2);
			VertexCompositeNtuple->Branch("pTerrD2", &candInfo.ptErr2);
			VertexCompositeNtuple->Branch("EtaD2", &candInfo.eta2);
			VertexCompositeNtuple->Branch("PhiD2", &candInfo.phi2);
			VertexCompositeNtuple->Branch("dedxHarmonic2D2", &candInfo.H2dedx2);
			VertexCompositeNtuple->Branch("Dtrk1Dz1", &candInfo.Dtrk1Dz1);
			VertexCompositeNtuple->Branch("Dtrk2Dz1", &candInfo.Dtrk2Dz1);
			VertexCompositeNtuple->Branch("Dtrk1Dxy1", &candInfo.Dtrk1Dxy1);
			VertexCompositeNtuple->Branch("Dtrk2Dxy1", &candInfo.Dtrk2Dxy1);
			VertexCompositeNtuple->Branch("Dtrk1DzError1", &candInfo.Dtrk1DzError1);
			VertexCompositeNtuple->Branch("Dtrk2DzError1", &candInfo.Dtrk2DzError1);
			VertexCompositeNtuple->Branch("Dtrk1DxyError1", &candInfo.Dtrk1DxyError1);
			VertexCompositeNtuple->Branch("Dtrk2DxyError1", &candInfo.Dtrk2DxyError1);

			if(threeProngDecay_){
				VertexCompositeNtuple->Branch("pTD3", &candInfo.pt3);
				VertexCompositeNtuple->Branch("pTerrD3", &candInfo.ptErr3);
				VertexCompositeNtuple->Branch("EtaD3", &candInfo.eta3);
				VertexCompositeNtuple->Branch("PhiD3", &candInfo.phi3);

				VertexCompositeNtuple->Branch("Dtrk3Dz1", &candInfo.Dtrk3Dz1);
				VertexCompositeNtuple->Branch("Dtrk3Dxy1", &candInfo.Dtrk3Dxy1);
				VertexCompositeNtuple->Branch("Dtrk3DzError1", &candInfo.Dtrk3DzError1);
				VertexCompositeNtuple->Branch("Dtrk3DxyError1", &candInfo.Dtrk3DxyError1);

			}

			if(doGenMatching_)
			{
				VertexCompositeNtuple->Branch("isSwap", &candInfo.isSwap);
				VertexCompositeNtuple->Branch("idmom_reco", &candInfo.idmom_reco);
				VertexCompositeNtuple->Branch("idD1_reco", &candInfo.idd1_reco);
				VertexCompositeNtuple->Branch("idD2_reco", &candInfo.idd2_reco);
				VertexCompositeNtuple->Branch("idD3_reco", &candInfo.idd3_reco);
				VertexCompositeNtuple->Branch("matchGEN", &candInfo.matchGEN);
				VertexCompositeNtuple->Branch("gen3DPointingAngle", &candInfo.gen_agl_abs);
				VertexCompositeNtuple->Branch("gen2DPointingAngle", &candInfo.gen_agl2D_abs);
				VertexCompositeNtuple->Branch("gen3DDecayLength", &candInfo.gen_dl);
				VertexCompositeNtuple->Branch("gen2DDecayLength", &candInfo.gen_dl2D);
			}

		}

	} // doRecoNtuple_

	if(doGenNtuple_)
	{
		VertexCompositeNtuple->Branch("pT_gen", &candInfo.pt_gen);
		VertexCompositeNtuple->Branch("eta_gen", &candInfo.eta_gen);
		VertexCompositeNtuple->Branch("y_gen", &candInfo.y_gen);
		VertexCompositeNtuple->Branch("phi_gen", &candInfo.phi_gen);
		VertexCompositeNtuple->Branch("MotherID_gen", &candInfo.idmom);	  
		if(decayInGen_)
		{ 
			VertexCompositeNtuple->Branch("DauID1_gen", &candInfo.iddau1);
			VertexCompositeNtuple->Branch("DauID2_gen", &candInfo.iddau2);
			VertexCompositeNtuple->Branch("DauID3_gen", &candInfo.iddau3);
		}
	}
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
VCTreeProducer_LamC3P::endJob() {

}


void VCTreeProducer_LamC3P::genDecayLength(const uint& it, const reco::GenParticle& gCand, CandidateData& candInfo) {

	candInfo.gen_dl[it] = -99.9;
	candInfo.gen_agl_abs[it] = -99.9;
	candInfo.gen_dl2D[it] = -99.9;
	candInfo.gen_agl2D_abs[it] = -99.9;

	if (gCand.numberOfDaughters() == 0 || !gCand.daughter(0)) return;

	const auto& dauVtx = gCand.daughter(0)->vertex();
	TVector3 ptosvec(dauVtx.X(), dauVtx.Y(), dauVtx.Z());
	TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());

	// Use the struct object to access and modify the vectors
	candInfo.gen_agl_abs[it] = secvec.Angle(ptosvec);
	candInfo.gen_dl[it] = ptosvec.Mag();

	TVector3 ptosvec2D(dauVtx.X(), dauVtx.Y(), 0.0);
	TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);

	// Use the struct object to access and modify the vectors
	candInfo.gen_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
	candInfo.gen_dl2D[it] = ptosvec2D.Mag();
}

//define this as a plug-in
DEFINE_FWK_MODULE(VCTreeProducer_LamC3P);
