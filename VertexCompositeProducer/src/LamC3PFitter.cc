
//+++++++++++++++++++++++++++++++++++

//Code owner: Nihar Ranjan saha
//contact mail: nihar.ranjan.saha@cern.ch
//Date: 24 Sept 2025

/*
Note: This code named as Fitter, where
the reconstruction of Lc candidates is
done out of all charged track from miniAOD data.
*/

//+++++++++++++++++++++++++++++++++++




#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include <vector>
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

const float piMassLamC3P = 0.13957018;
const float piMassLamC3PSquared = piMassLamC3P*piMassLamC3P;
const float kaonMassLamC3P = 0.493677;
const float kaonMassLamC3PSquared = kaonMassLamC3P*kaonMassLamC3P;
const float protonMassLamC3P = 0.938272013; 
const float protonMassLamC3PSquared = protonMassLamC3P*protonMassLamC3P;
const float lamCMassLamC3P = 2.28646;
float piMassLamC3P_sigma = 3.5E-7f;
float kaonMassLamC3P_sigma = 1.6E-5f;
float protonMassLamC3P_sigma = 1.6E-5f;
float lamCMassLamC3P_sigma = lamCMassLamC3P*1.e-6;

float cand1Mass[2] = {piMassLamC3P, protonMassLamC3P};
float cand2Mass[2] = {protonMassLamC3P, piMassLamC3P};
float cand1Mass_sigma[2] = {piMassLamC3P_sigma, protonMassLamC3P_sigma};
float cand2Mass_sigma[2] = {protonMassLamC3P_sigma, piMassLamC3P_sigma};



#define PROTON_MASS 0.9383
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677

//These are defined to limit the array size, to be consistent with Dfinder!
#define MAX_CAN       20000
#define MAX_Vertices 4000
#define MAX_TRACK    6000
#define MAX_MUON     10000
#define MAX_GEN      6000
#define MAX_BX       150
#define MAX_TRIGGER  30


using std::vector;
using std::string;
using namespace reco;
using namespace edm;
using namespace std;



// Constructor and (empty) destructor

LamC3PFitter::LamC3PFitter(const edm::ParameterSet& theParameters, edm::ConsumesCollector && iC, std::vector<std::vector<int>>& selectedTkhidxSetIn): selectedTkhidxSet(selectedTkhidxSetIn) , bField_esToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>()) 
{

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));

  token_packedCandidates = iC.consumes<std::vector<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));

  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
  chi2Map_  = iC.consumes< edm::ValueMap< float > >(theParameters.getParameter< edm::InputTag >( "TrackChi2Label" ) );


  // Second, initialize post-fit cuts
  mKPCutMin = theParameters.getParameter<double>(string("mKPCutMin"));
  mKPCutMax = theParameters.getParameter<double>(string("mKPCutMax"));
  mPiKPCutMin = theParameters.getParameter<double>(string("mPiKPCutMin"));
  mPiKPCutMax = theParameters.getParameter<double>(string("mPiKPCutMax"));
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  tkPtErrCut = theParameters.getParameter<double>(string("tkPtErrCut"));
  tkEtaCut = theParameters.getParameter<double>(string("tkEtaCut"));
//  tkPtSumCut = theParameters.getParameter<double>(string("tkPtSumCut"));
//  tkEtaDiffCut = theParameters.getParameter<double>(string("tkEtaDiffCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut2D = theParameters.getParameter<double>(string("collinearityCut2D"));
  collinCut3D = theParameters.getParameter<double>(string("collinearityCut3D"));
  lamCMassCut = theParameters.getParameter<double>(string("lamCMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPt3Cut = theParameters.getParameter<double>(string("dPt3Cut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  alpha2DCut = theParameters.getParameter<double>(string("alpha2DCut"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));

  //Nihar
  dPtCut_= theParameters.getParameter<double>(string("dPtCut"));
  dRapidityCut_= theParameters.getParameter<double>(string("dRapidityCut"));
  
  useAnyMVA_ = false;
  forestLabel_ = "LamC3PInpPb";
  std::string type = "BDT";
  useForestFromDB_ = true;
  dbFileName_ = "";

  forest_ = nullptr;

  if(theParameters.exists("useAnyMVA")) useAnyMVA_ = theParameters.getParameter<bool>("useAnyMVA");

  if(useAnyMVA_){
    if(theParameters.exists("mvaType"))type = theParameters.getParameter<std::string>("mvaType");
    if(theParameters.exists("GBRForestLabel"))forestLabel_ = theParameters.getParameter<std::string>("GBRForestLabel");
    if(theParameters.exists("GBRForestFileName")){
      dbFileName_ = theParameters.getParameter<std::string>("GBRForestFileName");
      useForestFromDB_ = false;
    }

    if(!useForestFromDB_){
      edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeProducer/data/%s",dbFileName_.c_str()));
      TFile gbrfile(fip.fullPath().c_str(),"READ");
      forest_ = (GBRForest*)gbrfile.Get(forestLabel_.c_str());
      gbrfile.Close();
    }

    mvaType_ = type;
    mvaToken_ = iC.esConsumes<GBRForest, GBRWrapperRcd>(edm::ESInputTag("", forestLabel_));
  }

  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }
}

LamC3PFitter::~LamC3PFitter() {
  delete forest_;
}






// Method containing the algorithm for vertex reconstruction


static const double Mass_in_permutation[16][3] = { 
  {  0.,              0.,                   0.   }, // 0 (- - -)
  {  0.,              0.,                   0.   },
  {  PROTON_MASS,     PION_MASS,            KAON_MASS  }, // 2 (- - +)  perm2+p=2
  {  PION_MASS,       PROTON_MASS,          KAON_MASS  },//perm2+p=3
  {  PROTON_MASS,     KAON_MASS,            PION_MASS }, // 4 (- + -)
  {  PION_MASS,       KAON_MASS,            PROTON_MASS  },
  {  KAON_MASS,       PROTON_MASS,          PION_MASS }, // 6 (- + +)
  {  KAON_MASS,       PION_MASS,            PROTON_MASS  },
  {  KAON_MASS,       PROTON_MASS,          PION_MASS }, // 8 (+ - -)
  {  KAON_MASS,       PION_MASS,            PROTON_MASS  },
  {  PROTON_MASS,     KAON_MASS,            PION_MASS }, // 10 (+ - +)
  {  PION_MASS,       KAON_MASS,            PROTON_MASS  },
  {  PROTON_MASS,     PION_MASS,            KAON_MASS  }, // 12 (+ + -)
  {  PION_MASS,       PROTON_MASS,          KAON_MASS  },
  {  0.,              0.,                   0.   }, // 14 (+ + +)
  {  0.,              0.,                   0.   }
};

static const int PDGID_in_permutation[16][3] = {
    { 0,            0,             0            }, // 0 (- - -)
    { 0,            0,             0            },
    { -PROTON_PDGID,  -PION_PDGID,    KAON_PDGID   }, // 2 (- - +)
    { -PION_PDGID,    -PROTON_PDGID,  KAON_PDGID   }, // 3
    { -PROTON_PDGID,  KAON_PDGID,    -PION_PDGID   }, // 4 (- + -)
    { -PION_PDGID,    KAON_PDGID,    -PROTON_PDGID }, // 5
    { -KAON_PDGID,    PROTON_PDGID,  PION_PDGID   }, // 6 (- + +)
    { -KAON_PDGID,    PION_PDGID,    PROTON_PDGID }, // 7
    { KAON_PDGID,    -PROTON_PDGID,  -PION_PDGID   }, // 8 (+ - -)
    { KAON_PDGID,    -PION_PDGID,    -PROTON_PDGID }, // 9
    { PROTON_PDGID,  -KAON_PDGID,    PION_PDGID   }, // 10 (+ - +)
    { PION_PDGID,    -KAON_PDGID,    PROTON_PDGID }, // 11
    { PROTON_PDGID,  PION_PDGID,    -KAON_PDGID   }, // 12 (+ + -)
    { PION_PDGID,    PROTON_PDGID,  -KAON_PDGID   }, // 13
    { 0,            0,             0            }, // 14 (+ + +)
    { 0,            0,             0            }  // 15
};






void LamC3PFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  using namespace std;
  using namespace edm;
  using namespace reco;


  // Prepare custom track vectors
  std::vector<Track> lst;
  std::vector<TrackXYZP2> lstXYZP2;

  // Handles
  Handle<pat::PackedCandidateCollection> packedHandle;
  Handle<VertexCollection> vertexHandle;
  Handle<BeamSpot> beamSpotHandle;
  edm::Handle<edm::ValueMap<float>> chi2Handle;
  
  // Get data
  iEvent.getByToken(token_packedCandidates, packedHandle);
  iEvent.getByToken(token_vertices, vertexHandle);
  iEvent.getByToken(token_beamSpot, beamSpotHandle);
  iEvent.getByToken(chi2Map_, chi2Handle);

  //if (!packedHandle.isValid() || packedHandle->empty()) return;

  // Magnetic field (needed for other things in your module)
  auto bFieldHandle = iSetup.getHandle(bField_esToken_);
  magField = bFieldHandle.product();


  //bool isVtxPV = 0;
  double xVtx = -99999.0;
  double yVtx = -99999.0;
  double zVtx = -99999.0;
  double xVtxError = -999.0;
  double yVtxError = -999.0;
  double zVtxError = -999.0;
  double dzvtx = -999, dxyvtx = -999;
  double dzerror = -999, dxyerror = -999;
  double trk_chi2 = -999;



  
  const reco::VertexCollection vtxCollection = *(vertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();

  if (vtxCollection.size() > 0 && !vtxPrimary->isFake() )
  {
    if (vtxCollection.size() >=MAX_Vertices){
      std::cout<<"ERROR: number of  Vertices exceeds the size of array!"<<std::endl;
      return;
    }
    //isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else
  {
    //isVtxPV = 0;
    xVtx = beamSpotHandle->position().x();
    yVtx = beamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = beamSpotHandle->BeamWidthX();
    yVtxError = beamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }
  math::XYZPoint bestvtx(xVtx, yVtx, zVtx);
  math::XYZPoint bestvtxErr(xVtxError, yVtxError, zVtxError);


  std::vector<pat::PackedCandidate> input_tracks;
  std::vector<pat::PackedCandidate> input_tracks_passed;


  input_tracks = *packedHandle;
  

    int i_tk=0;
    int passedTrk=0;

    std::vector<bool> isNeededTrack;
    ImpactParameters ip_params;
    
    //for (size_t i = 0; i < packedHandle->size(); ++i) {
    for (std::vector<pat::PackedCandidate>::const_iterator tk_it=input_tracks.begin(); tk_it != input_tracks.end(); tk_it++) {

      i_tk++;


      if (passedTrk >= MAX_TRACK){
	std::cout<<"ERROR: number of tracks exceeds the size of array!"<<std::endl;
	break;
      }
    
      isNeededTrack.push_back(false);

      // Make sure the candidate is charged and has tracking information    
      if (!tk_it->hasTrackDetails()) continue;
      if (abs(tk_it->charge()) != 1) continue;
      if (tk_it->pt() < tkPtCut) continue;
      if (fabs(tk_it->eta()) > tkEtaCut) continue;

      
      //if(tk_it->pseudoTrack().normalizedChi2() > tkChi2Cut) continue;
      //if(tk_it->pseudoTrack().numberOfValidHits() < tkNhitsCut) continue;

      if( !(tk_it->pseudoTrack().quality(reco::TrackBase::qualityByName("highPurity")))) continue;
      if( tk_it->pseudoTrack().ptError()/tk_it->pt()>tkPtErrCut ) continue;

      //LATER
    /*if (chi2Handle.isValid() && !chi2Handle.failedToGet()) {
      trk_chi2 = (*chi2Handle)[ptr];
    }
    else {trk_chi2 = cand.pseudoTrack().normalizedChi2();}
    */
    
    //trk_chi2 = tk_it->pseudoTrack().normalizedChi2();
    //cout<<"trk_chi2="<<trk_chi2<<endl;
      
    // Optional quality filters (you can expand this)
    //if (trk_chi2 < tkChi2Cut || cand.pseudoTrack().numberOfValidHits() < tkNhitsCut || (cand.bestTrack() && (cand.bestTrack()->ptError() / cand.pt()) >= tkPtErrCut)) continue;



    dzvtx = tk_it->pseudoTrack().dz(bestvtx);
    dxyvtx = tk_it->pseudoTrack().dxy(bestvtx);
    dzerror = TMath::Sqrt(tk_it->pseudoTrack().dzError()*tk_it->pseudoTrack().dzError() + zVtxError*zVtxError);
    dxyerror = TMath::Sqrt(tk_it->pseudoTrack().dxyError()*tk_it->pseudoTrack().dxyError() + xVtxError*yVtxError);

    ip_params = {dzvtx, dxyvtx, dzerror, dxyerror};

    
    if (fabs(dzvtx / dzerror) <= dauLongImpactSigCut) continue;
    if (fabs(dxyvtx / dxyerror) <= dauTransImpactSigCut) continue;

    
    //input_tracks.push_back(cand);
    //const auto& cand = *tk_it;
    input_tracks_passed.push_back(*tk_it);

    
    
    //lst.push_back(Track(cand.pt(), cand.eta(), cand.phi(), cand.charge()+1, input_tracks.begin()));

    isNeededTrack[tk_it-input_tracks.begin()] = true;
    
    lst.push_back(
		  Track(input_tracks[tk_it-input_tracks.begin()].pt(),
			input_tracks[tk_it-input_tracks.begin()].eta(),
			input_tracks[tk_it-input_tracks.begin()].phi(),
			input_tracks[tk_it-input_tracks.begin()].charge()+1,
			tk_it-input_tracks.begin()
			)
		  );
    

    passedTrk++;
			

    }//-----End track loops

    //int n = (int) lst.size();    
    for (int i = 0; i < (int)lst.size(); i ++) {
      TrackXYZP2 tr = lst[i];
      lstXYZP2.push_back(tr);
    }

    
    //std::cout<<"Input tracks="<<input_tracks.size()<<"  "<<"input tracks passed="<<input_tracks_passed.size()<<"  "<<"lst size="<<lst.size()<<"  "<<"lstXYZP2 size="<<lstXYZP2.size()<<std::endl;
  

  float mass_window[2] = {2.1, 2.5};  
  fitLamCCandidates(input_tracks, bestvtx, bestvtxErr, ip_params, lst, lstXYZP2, mass_window, iEvent, iSetup);

}




void LamC3PFitter::TkCombinationPermutation_Lc_v3(
						  const std::vector<pat::PackedCandidate> input_tracks, 
						  //const std::vector<edm::Ptr<pat::PackedCandidate>> input_tracks,
						  std::vector<Track> lst,
						  std::vector<TrackXYZP2> lstXYZP2,
						  float *mass_window,
						  std::vector<TrackSelected>& selectedTkhidxSet
						  
						  ){
  
  int tk1_hindex = -1;
  int tk2_hindex = -1;
  int tk3_hindex = -1;

	
	//auto t0 = std::chrono::high_resolution_clock::now();	
	int number_NeededTrack = (int) lst.size();

	//std::cout<<"number_NeededTrack="<<number_NeededTrack<<std::endl;
	
	for (int tk1idx = 0; tk1idx < number_NeededTrack; tk1idx++) {
	  const TrackXYZP2& tr1 = lstXYZP2[tk1idx];
	  tk1_hindex = tr1.index;
	  int perm1 = tr1.q;

		  
	  for (int tk2idx = tk1idx + 1; tk2idx < number_NeededTrack; tk2idx++) {
	    const TrackXYZP2& tr2 = lstXYZP2[tk2idx];
	    tk2_hindex = tr2.index;

	    int perm2 = (perm1 << 1) + tr2.q;
	    P3 p12(tr2);
	    p12 += tr1;

	
	    
	    for (int tk3idx = tk2idx + 1; tk3idx < number_NeededTrack; tk3idx++) {

	      const TrackXYZP2& tr3 = lstXYZP2[tk3idx];
	      tk3_hindex = tr3.index;

	      int perm3 = (perm2 << 1) + tr3.q;

	      if (perm3 == 0 || perm3 == 14) continue;
	      
	      P3 pD(tr3);
	      pD += p12;

		      
	      double ptD2 = pD.px * pD.px + pD.py * pD.py;


	      if (ptD2 < dPtCut_ * dPtCut_) continue;
	      
	      double pzD2 = pD.pz * pD.pz;
	      for (int p = 0; p < 2; p++) {
                double p0 = Functs.totE(perm3 + p, tr1.p2, tr2.p2, tr3.p2);

		double mD2 = p0 * p0 - ptD2 - pzD2;

		
			
		if (mD2 < (mass_window[0] * mass_window[0]) || mD2 > (mass_window[1] * mass_window[1])) continue;

		
                double mtD2 = mD2 + ptD2;
		

		if (pzD2 > sinh(dRapidityCut_) * sinh(dRapidityCut_) * mtD2) continue;
		
                selectedTkhidxSet.push_back(TrackSelected(tk1_hindex,tk2_hindex,tk3_hindex, perm3 + p));

		continue;

	      }
	    }
	  }
	}
	//std::cout<<"TkCombinationPermutation, selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
	return;
}


//Define all 


  void LamC3PFitter::fitLamCCandidates(
				       const std::vector<pat::PackedCandidate> input_tracks, 
				       const math::XYZPoint bestvtx,
				       const math::XYZPoint bestvtxErr,
				       const ImpactParameters ip,
				       vector<Track> lst,
				       vector<TrackXYZP2> lstXYZP2,
				       float * mass_window,
				       const edm::Event& iEvent,
				       const edm::EventSetup &iSetup
				       ){



	
    vector<TrackSelected> selectedTkhidxSet;
    TkCombinationPermutation_Lc_v3( input_tracks,lst,lstXYZP2, mass_window, selectedTkhidxSet);


    float chi = 0.;
    float ndf = 0.;
    
    //particle factory: produce transient tracks
    KinematicParticleFactoryFromTransientTrack pFactory;
    VirtualKinematicParticleFactory vFactory;
    //fitter for D
    KinematicParticleVertexFitter   tktk_fitter;
    RefCountedKinematicTree         tktk_VFT;
    RefCountedKinematicParticle     tktk_VFP;
    RefCountedKinematicVertex       tktk_VFPvtx;
    //constrain fit fitter
    KinematicConstrainedVertexFitter kcv_tktk_fitter;


    const MagneticField& bField = iSetup.getData(bField_esToken_);
    const MagneticField* field = &bField;
  

    
    TLorentzVector v4_tk;
    std::vector<TLorentzVector> tktk_4vecs;//fitted tks
    TLorentzVector tktk_4vec;//fitted D
    TLorentzVector unfitted_tktk_4vec;//unfitted D
    std::vector<TLorentzVector> tktkRes_4vecs;//fitted Res tks
    TLorentzVector tktkRes_4vec;//fitted Res
    TLorentzVector unfitted_tktkRes_4vec;//unfitted Res
    std::vector<RefCountedKinematicParticle> tktk_candidate;//input tracks to D fitter
    std::vector<RefCountedKinematicParticle> tktkRes_candidate;//input tracks to Res fitter
    std::vector<RefCountedKinematicParticle> tktkCands;//output tracks from D fitter
    std::vector<RefCountedKinematicParticle> tktkResCands;//output tracks from Res fitter
    TLorentzVector temp_vec;//for temporary usage

    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
    typedef ROOT::Math::SVector<double, 3> SVector3;




    int Trk_idx = 0;
    
    for(int i = 0; i < int(selectedTkhidxSet.size()); i++){
    
      //if(theLamC3Ps.size() >= MAX_CAN) break;      
      if(Trk_idx >= MAX_CAN) break; //To match with Dfinder!!!
      
      //clear before using
      v4_tk.Clear();
      tktk_4vecs.clear();
      tktk_4vec.Clear();
      unfitted_tktk_4vec.Clear();
      tktkRes_4vecs.clear();
      tktkRes_4vec.Clear();
      unfitted_tktkRes_4vec.Clear();
      tktk_candidate.clear();
      tktkRes_candidate.clear();
      tktkCands.clear();
      tktkResCands.clear();
      
      unfitted_tktk_4vec.SetPxPyPzE(0., 0., 0., 0.);
      unfitted_tktkRes_4vec.SetPxPyPzE(0., 0., 0., 0.);
      

      ParticleMass tk_mass;
      std::vector<int> pushbackTrkIdx;
      std::vector<int> triplet_charges;
      std::vector<int> triplet_pdgids;

      float tk_sigma;
      //int tk1_charge=0, tk2_charge=0, tk3_charge=0;
      int tk1_pdgid=0, tk2_pdgid=0, tk3_pdgid=0;
      




      //Here index_tk1,2,3 are the member names of the struct
      int idx1 = selectedTkhidxSet[i].index_tk1;
      int idx2 = selectedTkhidxSet[i].index_tk2;
      int idx3 = selectedTkhidxSet[i].index_tk3;
      int permIdx = selectedTkhidxSet[i].permutation_number;


      
      const pat::PackedCandidate* d_cand1 = &input_tracks[idx1];
      const pat::PackedCandidate* d_cand2 = &input_tracks[idx2];
      const pat::PackedCandidate* d_cand3 = &input_tracks[idx3];
      

      pat::PackedCandidate cand1 = *d_cand1;
      pat::PackedCandidate cand2 = *d_cand2;
      pat::PackedCandidate cand3 = *d_cand3;




      reco::TransientTrack tkTT1(input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack(), &(*field) );
      //reco::TransientTrack tkTT1(cand1.pseudoTrack(), &(*field) );
      if (tkTT1.isValid())
	{
	  tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0];
	  tk1_pdgid = PDGID_in_permutation[selectedTkhidxSet[i].permutation_number][0];
	  tk_sigma = Functs.getParticleSigma(tk_mass);
	  tktk_candidate.push_back(pFactory.particle(tkTT1,tk_mass,chi,ndf,tk_sigma));
	}
	

      reco::TransientTrack tkTT2(input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack(), &(*field) );
      //reco::TransientTrack tkTT2(cand2.pseudoTrack(), &(*field) );
      if (tkTT2.isValid())
	{

	  tk_mass   = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1];
	  tk2_pdgid = PDGID_in_permutation[selectedTkhidxSet[i].permutation_number][1];
	  tk_sigma = Functs.getParticleSigma(tk_mass);
	  tktk_candidate.push_back(pFactory.particle(tkTT2, tk_mass, chi, ndf, tk_sigma));
	  
	  }
	

      reco::TransientTrack tkTT3(input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack(), &(*field) );
      //reco::TransientTrack tkTT3(cand3.pseudoTrack(), &(*field) );
      if (tkTT3.isValid())
	{
	  tk_mass   = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2];
	  tk3_pdgid = PDGID_in_permutation[selectedTkhidxSet[i].permutation_number][2];
	  tk_sigma = Functs.getParticleSigma(tk_mass);
	  tktk_candidate.push_back(pFactory.particle(tkTT3, tk_mass, chi, ndf, tk_sigma));
	  }

      
      double px1 = input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack().px();
      double py1 = input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack().py();
      double pz1 = input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack().pz();
      double m1 = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0];
      double E1 = sqrt(m1*m1 + px1*px1 + py1*py1 + pz1*pz1);
      
      double px2 = input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack().px();
      double py2 = input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack().py();
      double pz2 = input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack().pz();
      double m2 = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1];
      double E2 = sqrt(m2*m2 + px2*px2 + py2*py2 + pz2*pz2);
      
      double px3 = input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack().px();
      double py3 = input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack().py();
      double pz3 = input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack().pz();
      double m3 = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2];
      double E3 = sqrt(m3*m3 + px3*px3 + py3*py3 + pz3*pz3);
      
      double invmasssquare = (E1+E2+E3)*(E1+E2+E3) - (px1+px2+px3)*(px1+px2+px3) - (py1+py2+py3)*(py1+py2+py3) - (pz1+pz2+pz3)*(pz1+pz2+pz3);
      double invmass = sqrt(invmasssquare);
      //std::cout<<"The inv mass : "<<invmass<<std::endl;
    


      //double MaximumDoca = Functs.getMaxDoca(tktk_candidate);
      tktk_VFT = tktk_fitter.fit(tktk_candidate);
      
      if (!tktk_VFT || !tktk_VFT->isValid()) continue;
      

      
      
      tktk_VFT->movePointerToTheTop();
      RefCountedKinematicParticle LamC3Pcand = tktk_VFT->currentParticle();
      RefCountedKinematicVertex LamC3P_DecayVertex = tktk_VFT->currentDecayVertex();

      if (!LamC3Pcand->currentState().isValid())continue;
      
      if (!LamC3P_DecayVertex || !LamC3P_DecayVertex->vertexIsValid()) continue;
      
      double chi2_prob_tktk = TMath::Prob(LamC3P_DecayVertex->chiSquared(),LamC3P_DecayVertex->degreesOfFreedom());
      if( chi2_prob_tktk < VtxChiProbCut || chi2_prob_tktk > 1.0) continue; 

      

      // Convert GlobalPoint (float-based) to reco::Vertex::Point (double-based)
      reco::Vertex::Point position(LamC3P_DecayVertex->vertexState().position().x(),
				   LamC3P_DecayVertex->vertexState().position().y(),
				   LamC3P_DecayVertex->vertexState().position().z()
				   );

      GlobalVector LamC3P_TotalP = GlobalVector(LamC3Pcand->currentState().globalMomentum().x(),
						LamC3Pcand->currentState().globalMomentum().y(),
						LamC3Pcand->currentState().globalMomentum().z()
						);
      
      //float LamC3P_TotalE = sqrt(LamC3P_TotalP.mag2() + pow(invmass, 2));
        
      /*const reco::Particle::LorentzVector lamCP4(LamC3Pcand->currentState().globalMomentum().x(),
                                           LamC3Pcand->currentState().globalMomentum().y(),
                                           LamC3Pcand->currentState().globalMomentum().z(),
                                           LamC3P_TotalE);
      */
      const reco::Particle::LorentzVector lamCP4(LamC3Pcand->currentState().kinematicParameters().momentum().x(),
                                           LamC3Pcand->currentState().kinematicParameters().momentum().y(),
                                           LamC3Pcand->currentState().kinematicParameters().momentum().z(),
                                           LamC3Pcand->currentState().kinematicParameters().energy());


      Particle::Point LamC3P_Vtx((*LamC3P_DecayVertex).position().x(), (*LamC3P_DecayVertex).position().y(), (*LamC3P_DecayVertex).position().z());

        std::vector<double> LamC3P_VtxEVec;
        LamC3P_VtxEVec.push_back(LamC3P_DecayVertex->error().cxx());
        LamC3P_VtxEVec.push_back(LamC3P_DecayVertex->error().cyx());
        LamC3P_VtxEVec.push_back(LamC3P_DecayVertex->error().cyy());
        LamC3P_VtxEVec.push_back(LamC3P_DecayVertex->error().czx());
        LamC3P_VtxEVec.push_back(LamC3P_DecayVertex->error().czy());
        LamC3P_VtxEVec.push_back(LamC3P_DecayVertex->error().czz());
        SMatrixSym3D LamC3P_VtxCovMatrix(LamC3P_VtxEVec.begin(), LamC3P_VtxEVec.end());
        const Vertex::CovarianceMatrix LamC3P_VtxCov(LamC3P_VtxCovMatrix);
        double LamC3P_VtxChi2(LamC3P_DecayVertex->chiSquared());
        double LamC3P_VtxNdof(LamC3P_DecayVertex->degreesOfFreedom());
        double LamC3P_NormalizedChi2 = LamC3P_VtxChi2 / LamC3P_VtxNdof;

        double rVtxMag = 99999.0;
        double lVtxMag = 99999.0;
        double sigmaRvtxMag = 999.0;
        double sigmaLvtxMag = 999.0;
        double LamC3P_Angle3D = -100.0;
        double LamC3P_Angle2D = -100.0;

        GlobalVector LamC3P_LineOfFlight = GlobalVector(LamC3P_Vtx.x() - bestvtx.x(), LamC3P_Vtx.y() - bestvtx.y(), LamC3P_Vtx.z() - bestvtx.z());

	edm::Handle<reco::VertexCollection> vertexHandle;
	iEvent.getByToken(token_vertices, vertexHandle);
	const reco::VertexCollection vtxCollection = *(vertexHandle.product());
	reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();

	Handle<BeamSpot> beamSpotHandle;
	iEvent.getByToken(token_beamSpot, beamSpotHandle);
	//bool isVtxPV= True;
	
	SMatrixSym3D LamC3P_TotalCov;
        if (vtxCollection.size() > 0 && !vtxPrimary->isFake())
          LamC3P_TotalCov = LamC3P_VtxCovMatrix + vtxPrimary->covariance();
        else
          LamC3P_TotalCov = LamC3P_VtxCovMatrix + beamSpotHandle->rotatedCovariance3D();

        SVector3 distanceVector3D(LamC3P_LineOfFlight.x(), LamC3P_LineOfFlight.y(), LamC3P_LineOfFlight.z());
        SVector3 distanceVector2D(LamC3P_LineOfFlight.x(), LamC3P_LineOfFlight.y(), 0.0);

        LamC3P_Angle3D = angle(LamC3P_LineOfFlight.x(), LamC3P_LineOfFlight.y(), LamC3P_LineOfFlight.z(),
                          LamC3P_TotalP.x(), LamC3P_TotalP.y(), LamC3P_TotalP.z());
        LamC3P_Angle2D = angle(LamC3P_LineOfFlight.x(), LamC3P_LineOfFlight.y(), (float)0.0,
                          LamC3P_TotalP.x(), LamC3P_TotalP.y(), (float)0.0);

        lVtxMag = LamC3P_LineOfFlight.mag();
        rVtxMag = LamC3P_LineOfFlight.perp();
        sigmaLvtxMag = sqrt(ROOT::Math::Similarity(LamC3P_TotalCov, distanceVector3D)) / lVtxMag;
        sigmaRvtxMag = sqrt(ROOT::Math::Similarity(LamC3P_TotalCov, distanceVector2D)) / rVtxMag;

        if (LamC3P_NormalizedChi2 > chi2Cut ||
            rVtxMag < rVtxCut ||
            rVtxMag / sigmaRvtxMag < rVtxSigCut ||
            lVtxMag < lVtxCut ||
            lVtxMag / sigmaLvtxMag < lVtxSigCut ||
            cos(LamC3P_Angle3D) < collinCut3D || cos(LamC3P_Angle2D) < collinCut2D || LamC3P_Angle3D > alphaCut || LamC3P_Angle2D > alpha2DCut)
          continue;


	AnalyticalImpactPointExtrapolator extrapolator(magField);
	//TransverseImpactPointExtrapolator extrapolator(magField); 
	TrajectoryStateOnSurface tsos = extrapolator.extrapolate(LamC3Pcand->currentState().freeTrajectoryState(), RecoVertex::convertPos(vtxPrimary->position()));
        if (!tsos.isValid())continue;
	
        Measurement1D cur3DIP;
        VertexDistance3D a3d;
        GlobalPoint refPoint = tsos.globalPosition();
        GlobalError refPointErr = tsos.cartesianError().position();
        GlobalPoint vertexPosition = RecoVertex::convertPos(vtxPrimary->position());
        GlobalError vertexPositionErr = RecoVertex::convertError(vtxPrimary->error());
        cur3DIP = (a3d.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));



	int lamCCharge = cand1.charge() + cand2.charge() + cand3.charge();
	if (abs(lamCCharge) != 1) continue;  
	int LamC3P_pdgID = (lamCCharge == 1) ? 4122 : -4122;
      
      
      // Build the Lambda_c candidate
      pat::CompositeCandidate theLamC3P;


      theLamC3P.setPdgId(LamC3P_pdgID);
      theLamC3P.setCharge(lamCCharge);
      theLamC3P.setP4(lamCP4);



      cand1.setPdgId(tk1_pdgid);
      cand2.setPdgId(tk2_pdgid);
      cand3.setPdgId(tk3_pdgid);
      
      //cand1.setPdgId();
      theLamC3P.addDaughter(cand1);
      theLamC3P.addDaughter(cand2);
      theLamC3P.addDaughter(cand3);

      //Add user variables!      
      theLamC3P.addUserFloat("vtxX", LamC3P_Vtx.x());
      theLamC3P.addUserFloat("vtxY", LamC3P_Vtx.y());
      theLamC3P.addUserFloat("vtxZ", LamC3P_Vtx.z());
      theLamC3P.addUserFloat("bestvtxX", bestvtx.x());
      theLamC3P.addUserFloat("bestvtxY", bestvtx.y());
      theLamC3P.addUserFloat("bestvtxZ", bestvtx.z());
      theLamC3P.addUserFloat("xVtxError", bestvtxErr.x());
      theLamC3P.addUserFloat("yVtxError", bestvtxErr.y());
      theLamC3P.addUserFloat("zVtxError", bestvtxErr.z());
      theLamC3P.addUserFloat("vtxChi2", LamC3P_VtxChi2);
      theLamC3P.addUserFloat("vtxNdof", LamC3P_VtxNdof);
      theLamC3P.addUserFloat("dzvtx", ip.dz);
      theLamC3P.addUserFloat("dxyvtx", ip.dxy);
      theLamC3P.addUserFloat("dzError", ip.dzError);
      theLamC3P.addUserFloat("dxyError", ip.dxyError);
      //Newly added
      theLamC3P.addUserFloat("ip3d", cur3DIP.value()); 
      theLamC3P.addUserFloat("ip3derr", cur3DIP.error()); 
      
      
      
	
      for (int ii = 0; ii < 3; ++ii)
        {
          for (int jj = 0; jj < 3; ++jj)
	    {
	      theLamC3P.addUserFloat("vertexCovariance_" + std::to_string(ii) + "_" + std::to_string(jj), LamC3P_VtxCov(ii, jj));
	    }
        }
      
      //std::cout<<"LamC mass="<<theLamC3P.mass()<<std::endl;
      theLamC3Ps.push_back(theLamC3P);


      
      Trk_idx++;
    }//---End of Track index loop ----
    
    
  }// End function
  


const pat::CompositeCandidateCollection& LamC3PFitter::getLamC3P() const {
  return theLamC3Ps;
}


const std::vector<float>& LamC3PFitter::getMVAVals() const {
  return mvaVals_;
}


void LamC3PFitter::resetAll() {
    theLamC3Ps.clear();
    mvaVals_.clear();
}
