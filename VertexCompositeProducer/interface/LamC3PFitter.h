
#ifndef VertexCompositeAnalysis__LAMC3P_FITTER_H
#define VertexCompositeAnalysis__LAMC3P_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


#include <string>
#include <fstream>
#include <typeinfo>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include "TLorentzVector.h"

#define ELECTRON_MASS 0.0005
#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916
#define PSI2S_MASS  3.686109
#define PROTON_MASS 0.9383
#define D0_MASS 1.8648
#define DSTAR_MASS 2.01028
#define DPLUS_MASS 1.8696
#define DSUBS_MASS 1.9683
#define BPLUS_MASS 5.27931
#define BZERO_MASS 5.27962
#define BSUBS_MASS 5.36682
#define F0980_MASS 0.99
#define LAMBDAC_MASS 2.286
#define KSTAR892_MASS 0.89581
#define DELTA1232PLUSPLUS_MASS 1.209
#define LAMBDA1520_MASS 1.5195

#define MUON_PDGID 13
#define PION_PDGID 211
#define KAON_PDGID 321
#define KSTAR_PDGID 313
#define PHI_PDGID 333
#define JPSI_PDGID 443
#define BZERO_PDGID 511
#define BPLUS_PDGID 521
#define BSUBS_PDGID 531
#define DZERO_PDGID 421
#define DPLUS_PDGID 411
#define DSUBS_PDGID 431
#define DSTAR_PDGID 413
#define F0980_PDGID 9010221
#define PROTON_PDGID 2212
#define LAMBDAC_PDGID 4122
#define KSTAR892_PDGID 313
#define DELTA1232PLUSPLUS_PDGID 2224
#define LAMBDA1520_PDGID  3124
#define CHIC1_PDGID 20443
#define X_PDGID 9920443
#define PSI2S_PDGID 100443



class CommonFuncts {
	public:
		void test() {}

		float getParticleSigma(double mass) {
			if (mass == ELECTRON_MASS) return 0.013E-9f;
			if (mass == MUON_MASS) return 4E-9f;
			if (mass == PION_MASS) return 3.5E-7f;
			if (mass == KAON_MASS) return 1.6E-5f;
			if (mass == PROTON_MASS) return 8E-8f;
			return 1E-6f;
		}

		const double mp2 = 0.9383 * 0.9383;
		const double mK2 = 0.493677 * 0.493677;
		const double mpi2 = 0.13957 * 0.13957;

		double MQtbl[16][3] = {
			{0., 0., 0.}, {0., 0., 0.},
			{mp2, mpi2, mK2}, {mpi2, mp2, mK2},
			{mp2, mK2, mpi2}, {mpi2, mK2, mp2},
			{mK2, mp2, mpi2}, {mK2, mpi2, mp2},
			{mK2, mp2, mpi2}, {mK2, mpi2, mp2},
			{mp2, mK2, mpi2}, {mpi2, mK2, mp2},
			{mp2, mpi2, mK2}, {mpi2, mp2, mK2},
			{0., 0., 0.}, {0., 0., 0.}
		};

		double totE(int p, double p1sq, double p2sq, double p3sq) {
			double m1 = MQtbl[p][0];
			double m2 = MQtbl[p][1];
			double m3 = MQtbl[p][2];
			return sqrt(p1sq + m1) + sqrt(p2sq + m2) + sqrt(p3sq + m3);
		}
};









class LamC3PFitter {
	public:
		LamC3PFitter(const edm::ParameterSet& theParams, edm::ConsumesCollector&& iC, std::vector<std::vector<int>>& selectedTkhidxSetIn);
		~LamC3PFitter();


		struct Track {
			double pt, eta, phi;
			int q;
			int index;

			Track(double ptp, double etap, double phip, int qp, int indexp) : pt(ptp), eta(etap), phi(phip), q(qp), index(indexp) { }

		};


		struct TrackXYZ {
			double px, py, pz;
			int q;
			int index;

			TrackXYZ(double pt, double eta, double phi, int qp, int indexp) : q(qp), index(indexp) { setup(pt, eta, phi); }

			TrackXYZ(const Track& tr) : q(tr.q) , index(tr.index) { setup(tr.pt, tr.eta, tr.phi); }

			void setup(double pt, double eta, double phi) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
			}

		};


		struct TrackXYZP2 {
			double px, py, pz, p2;
			int q;
			int index;

			TrackXYZP2(double pt, double eta, double phi, int qp, int indexp) : q(qp),index(indexp) { setup(pt, eta, phi); }

			TrackXYZP2(const Track& tr) : q(tr.q), index(tr.index) { setup(tr.pt, tr.eta, tr.phi); }

			TrackXYZP2(const TrackXYZ& tr) : px(tr.px), py(tr.py), pz(tr.pz), q(tr.q), index(tr.index) {
				p2 = px * px + py * py + pz * pz;
			}

			void setup(double pt, double eta, double phi) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
				p2 = pt * pt + pz * pz;
			}

		};

		struct Triplet {
			int i1, i2, i3;

			Triplet(int i1p, int i2p, int i3p) : i1(i1p), i2(i2p), i3(i3p) { }
		};

		struct TrackSelected{
			int index_tk1, index_tk2, index_tk3, permutation_number;

			TrackSelected(int index_tk1p, int index_tk2p, int index_tk3p, int permutation_numberp) : index_tk1(index_tk1p), index_tk2(index_tk2p), index_tk3(index_tk3p), permutation_number(permutation_numberp) { }
		};

		struct P3 {
			double px, py, pz;

			P3(double pt, double eta, double phi) {
				setup(pt, eta, phi);
			}

			P3(const Track& tr) {
				setup(tr.pt, tr.eta, tr.phi);
			}

			P3(const TrackXYZ& tr) {
				px = tr.px;
				py = tr.py;
				pz = tr.pz;
			}

			P3(const TrackXYZP2& tr) {
				px = tr.px;
				py = tr.py;
				pz = tr.pz;
			}

			void setup(double pt, double eta, double phi) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
			}

			P3& operator += (const P3& a) {
				px += a.px;
				py += a.py;
				pz += a.pz;
				return *this;
			}

		};//P3

		struct P4 {
			double px, py, pz, p0;

			P4(double pt, double eta, double phi, double m) {
				setup(pt, eta, phi, m);
			}

			P4(const Track& tr, double m) {
				setup(tr.pt, tr.eta, tr.phi, m);
			}

			P4(const TrackXYZ& tr, double m) {
				px = tr.px;
				py = tr.py;
				pz = tr.pz;
				double p2 = px * px + py * py + pz * pz;
				p0 = sqrt(p2 + m * m);
			}

			void setup(double pt, double eta, double phi, double m) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
				double p2 = pt * pt + pz * pz;
				p0 = sqrt(p2 + m * m);
			}

			P4& operator += (const P4& a) {
				px += a.px;
				py += a.py;
				pz += a.pz;
				p0 += a.p0;
				return *this;
			}

		};//P4  

		struct ImpactParameters {
			double dz;
			double dxy;
			double dzError;
			double dxyError;
		};



		std::vector<Track> lst;
		std::vector<TrackXYZP2> lstXYZP2;

		std::vector<int> selectedTkhidx;
		std::vector< std::vector<int> > &selectedTkhidxSet;
		TLorentzVector temp_vec;


		void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);

		CommonFuncts Functs;

		void TkCombinationPermutation_Lc_v3(
				const std::vector<pat::PackedCandidate> input_tracks,
				std::vector<LamC3PFitter::Track> lst,
				std::vector<LamC3PFitter::TrackXYZP2> lstXYZP2,
				float* mass_window,
				std::vector<LamC3PFitter::TrackSelected>& selectedTkhidxSet
				);



		void fitLamCCandidates(
				const std::vector<pat::PackedCandidate> input_tracks,
				const math::XYZPoint bestvtx,
				const math::XYZPoint bestvtxErr,
				const ImpactParameters ip_params,
				std::vector<Track> lst,
				std::vector<TrackXYZP2> lstXYZP2,
				float *mass_window,
				const edm::Event& iEvent,
				const edm::EventSetup &iSetup
				);



		// Switching to L. Lista's reco::Candidate infrastructure for LamC3P storage
		const pat::CompositeCandidateCollection& getLamC3P() const;
		const std::vector<float>& getMVAVals() const; 

		void resetAll();

	private:
		// STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
		pat::CompositeCandidateCollection theLamC3Ps;

		// Tracker geometry for discerning hit positions
		const TrackerGeometry* trackerGeom;

		const MagneticField* magField;

		edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bField_esToken_;

		edm::InputTag recoAlg;
		edm::InputTag vtxAlg;
		edm::EDGetTokenT<reco::TrackCollection> token_tracks;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> token_packedCandidates;

		edm::EDGetTokenT<reco::VertexCollection> token_vertices;
		edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > token_dedx;
		edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
		edm::EDGetTokenT< edm::ValueMap< float > > chi2Map_;

		// Cuts
		double mKPCutMin;
		double mKPCutMax;
		double mPiKPCutMin;
		double mPiKPCutMax;
		double tkDCACut;
		double tkChi2Cut;
		int    tkNhitsCut;
		double tkPtErrCut;
		double tkPtCut;
		double tkEtaCut;
		double tkPtSumCut;
		double tkEtaDiffCut;
		double chi2Cut;
		double rVtxCut;
		double rVtxSigCut;
		double lVtxCut;
		double lVtxSigCut;
		double collinCut2D;
		double collinCut3D;
		double lamCMassCut;
		double dauTransImpactSigCut;
		double dauLongImpactSigCut;
		double VtxChiProbCut;
		double dPt3Cut;
		double alphaCut;
		double alpha2DCut;
		bool   isWrongSign;

		//Nihar
		double dPtCut_;
		double dRapidityCut_;


		std::vector<reco::TrackBase::TrackQuality> qualities;

		bool useAnyMVA_;
		std::vector<bool> useMVA_;
		std::vector<double> min_MVA_;
		std::string mvaType_;
		std::string forestLabel_;
		GBRForest * forest_;
		bool useForestFromDB_;

		std::vector<float> mvaVals_;
		edm::ESGetToken<GBRForest, GBRWrapperRcd> mvaToken_;

		std::string dbFileName_;

};




#endif
