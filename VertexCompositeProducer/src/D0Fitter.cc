////////////////////
// Orgininal author: Abby Wesolek
// last updates: 01 October 2025
// contact:abigail.leigh.wesolek@cern.ch
////////////////////
////////////////////
////////////////////
// this fitter.cc reconstructs D0 candidates from 2 tracks in accordance with the D0->k,pi channnel
// theD0 is a pat::CompositeCandidate
// the daughters are pat::PackedCandidates
// the output from this file can be saved to the edm.root file
// i.e.) vector<pat::CompositeCandidate>      "generalD0CandidatesNew"   "D0"
// this can be configured in (VertexCompositeProducer/test/run_edm_and_ttree_DATA_forD0.py)
// note the output from this file is all reconstructed D0 candidates, cuts are applied next in the selector.cc
////////////////////

// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      D0Fitter
//
/**\class D0Fitter D0Fitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/D0Fitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/D0Fitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TLorentzVector.h"

const float piMassD0 = 0.13957018;
const float piMassD0Squared = piMassD0 * piMassD0;
const float kaonMassD0 = 0.493677;
const float kaonMassD0Squared = kaonMassD0 * kaonMassD0;
const float d0MassD0 = 1.86484;
float piMassD0_sigma = 3.5E-7f;
float kaonMassD0_sigma = 1.6E-5f;

// Constructor and (empty) destructor
D0Fitter::D0Fitter(const edm::ParameterSet &theParameters, edm::ConsumesCollector &&iC) : bField_esToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks_pf = iC.consumes<std::vector<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData>>(edm::InputTag("dedxHarmonic2"));
  chi2Map_ = iC.consumes<edm::ValueMap<float>>(theParameters.getParameter<edm::InputTag>("TrackChi2Label"));

  // Second, initialize post-fit cuts
  mPiKCutMin = theParameters.getParameter<double>(string("mPiKCutMin"));
  mPiKCutMax = theParameters.getParameter<double>(string("mPiKCutMax"));
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  tkPtErrCut = theParameters.getParameter<double>(string("tkPtErrCut"));
  tkEtaCut = theParameters.getParameter<double>(string("tkEtaCut"));
  tkPtSumCut = theParameters.getParameter<double>(string("tkPtSumCut"));
  tkEtaDiffCut = theParameters.getParameter<double>(string("tkEtaDiffCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut2D = theParameters.getParameter<double>(string("collinearityCut2D"));
  collinCut3D = theParameters.getParameter<double>(string("collinearityCut3D"));
  d0MassCut = theParameters.getParameter<double>(string("d0MassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPtCut = theParameters.getParameter<double>(string("dPtCut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  alpha2DCut = theParameters.getParameter<double>(string("alpha2DCut"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));
}

D0Fitter::~D0Fitter()
{
}

// Method containing the algorithm for vertex reconstruction
void D0Fitter::fitAll(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  // if (iEvent.id().event() != 42371649) return;
  using std::cout;
  using std::endl;
  using std::vector;
  using namespace reco;
  using namespace edm;
  using namespace std;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  Handle<pat::PackedCandidateCollection> packedHandle;
  Handle<VertexCollection> theVertexHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;

  ESHandle<MagneticField> bFieldHandle;
  Handle<edm::ValueMap<reco::DeDxData>> dEdxHandle;
  edm::Handle<edm::ValueMap<float>> chi2Handle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks_pf, packedHandle);
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
  iEvent.getByToken(token_dedx, dEdxHandle);
  iEvent.getByToken(chi2Map_, chi2Handle);

  if (!packedHandle.isValid() || packedHandle->empty())
    return;

  // magnetic field needec for other things
  bFieldHandle = iSetup.getHandle(bField_esToken_);
  magField = bFieldHandle.product();

  bool isVtxPV = 0;
  double xVtx = -99999.0;
  double yVtx = -99999.0;
  double zVtx = -99999.0;
  double xVtxError = -999.0;
  double yVtxError = -999.0;
  double zVtxError = -999.0;
  double dzvtx = -999, dxyvtx = -999;
  double dzerror = -999, dxyerror = -999;
  double trk_chi2 = -999;

  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();

  if (vtxCollection.size() > 0 && !vtxPrimary->isFake())
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else
  {
    // cout << "fake vertex" << endl;
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }
  math::XYZPoint bestvtx(xVtx, yVtx, zVtx);

  // Loop over PackedCandidates and apply preselection cuts to identify good D0 daughter tracks

  std::vector<edm::Ptr<pat::PackedCandidate>> input_daughter_tracks; // today change, we are chaging from copy of pat::packedCand to pointers of the real packedCandidates, helpful for memory management and speed

  // --note input_daughter_tracks holds smart pointers edm::Ptr<> to them orginal objects
  for (size_t i = 0; i < packedHandle->size(); ++i) // for all tracks
  {
    edm::Ptr<pat::PackedCandidate> ptr(packedHandle, i);
    const auto &cand = *ptr; // get track

    // Make sure the candidate is charged and has tracking information
    if (cand.charge() == 0 || !cand.hasTrackDetails())
      continue;

    // Basic kinematic cuts
    if (cand.pt() < tkPtCut || fabs(cand.eta()) > tkEtaCut)
      continue;

    if (chi2Handle.isValid() && !chi2Handle.failedToGet())
    {
      trk_chi2 = (*chi2Handle)[ptr];
    }
    else
    {
      trk_chi2 = cand.pseudoTrack().normalizedChi2();
    }

    // Optional quality proxies via pseudoTrack and bestTrack

    // if  (cand.pseudoTrack().normalizedChi2() < tkChi2Cut &&
    if (trk_chi2 < tkChi2Cut &&
        cand.pseudoTrack().numberOfValidHits() >= tkNhitsCut &&
        cand.pseudoTrack().ptError() / cand.pt() < tkPtErrCut &&
        cand.pt() > tkPtCut && fabs(cand.eta()) < tkEtaCut)
    {

      // Impact parameter significance cuts
      // primary cahnges for today!!
      dzvtx = cand.pseudoTrack().dz(bestvtx);
      dxyvtx = cand.pseudoTrack().dxy(bestvtx);
      dzerror = TMath::Sqrt(cand.pseudoTrack().dzError() * cand.pseudoTrack().dzError() + zVtxError * zVtxError);
      dxyerror = TMath::Sqrt(cand.pseudoTrack().dxyError() * cand.pseudoTrack().dxyError() + xVtxError * yVtxError);
      double dauLongImpactSig = dzvtx / dzerror;
      double dauTransImpactSig = dxyvtx / dxyerror;

      if (fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut)
        input_daughter_tracks.push_back(ptr);
    }
  } // prelimanry loop

  float posCandMass[2] = {piMassD0, kaonMassD0};
  float negCandMass[2] = {kaonMassD0, piMassD0};
  float posCandMass_sigma[2] = {piMassD0_sigma, kaonMassD0_sigma};
  float negCandMass_sigma[2] = {kaonMassD0_sigma, piMassD0_sigma};
  int pdg_id[2] = {421, -421};
  int pos_pdg_id[2] = {211, 321};
  int neg_pdg_id[2] = {-321, -211};

  for (unsigned int trdx1 = 0; trdx1 < input_daughter_tracks.size(); trdx1++)
  {
    for (unsigned int trdx2 = trdx1 + 1; trdx2 < input_daughter_tracks.size(); trdx2++)
    {
      edm::Ptr<pat::PackedCandidate> dau1 = input_daughter_tracks[trdx1];
      edm::Ptr<pat::PackedCandidate> dau2 = input_daughter_tracks[trdx2];

      pat::PackedCandidate tk1, tk2;
      if (dau1->charge() > 0 && dau2->charge() < 0)
      {
        tk1 = *dau1;
        tk2 = *dau2;
      }
      else if (dau1->charge() < 0 && dau2->charge() > 0)
      {
        tk1 = *dau2;
        tk2 = *dau1;
      }
      else
        continue;

      if (tk1.pt() + tk2.pt() < tkPtSumCut)
        continue;
      if (fabs(tk1.eta() - tk2.eta()) > tkEtaDiffCut)
        continue;

      int q1 = tk1.charge();
      int q2 = tk2.charge();

      if (!isWrongSign && q1 * q2 != -1)
        continue;
      if (isWrongSign && q1 * q2 != 1)
        continue;

      reco::TransientTrack tt1(tk1.pseudoTrack(), magField);
      reco::TransientTrack tt2(tk2.pseudoTrack(), magField);

      if (!tt1.impactPointTSCP().isValid() || !tt2.impactPointTSCP().isValid())
        continue;

      // DCA calculation
      FreeTrajectoryState state1 = tt1.impactPointTSCP().theState();
      FreeTrajectoryState state2 = tt2.impactPointTSCP().theState();

      ClosestApproachInRPhi cApp;
      cApp.calculate(state1, state2);
      if (!cApp.status())
        continue;

      float dca = fabs(cApp.distance());
      GlobalPoint cxPt = cApp.crossingPoint();
      if (dca < 0. || dca > tkDCACut)
        continue;

      if (sqrt(cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.)
        continue;

      // Assign mass hypotheses for the two combinations: [K-, π+] and [π-, K+]

      TLorentzVector p4_kaon, p4_pion;

      if (!isWrongSign && q1 < 0 && q2 > 0)
      {
        p4_kaon.SetPtEtaPhiM(tk1.pt(), tk1.eta(), tk1.phi(), kaonMassD0); // tk1 = K-
        p4_pion.SetPtEtaPhiM(tk2.pt(), tk2.eta(), tk2.phi(), piMassD0);   // tk2 = π+
      }
      else if (!isWrongSign && q1 > 0 && q2 < 0)
      {
        p4_kaon.SetPtEtaPhiM(tk2.pt(), tk2.eta(), tk2.phi(), kaonMassD0); // tk2 = K+
        p4_pion.SetPtEtaPhiM(tk1.pt(), tk1.eta(), tk1.phi(), piMassD0);   // tk1 = π-
      }

      TLorentzVector sum = p4_kaon + p4_pion;

      auto ts1 = tt1.trajectoryStateClosestToPoint(cxPt);
      auto ts2 = tt2.trajectoryStateClosestToPoint(cxPt);

      double totalE1 = sqrt(ts1.momentum().mag2() + kaonMassD0Squared) +
                       sqrt(ts2.momentum().mag2() + piMassD0Squared);
      double totalE1Sq = totalE1 * totalE1;

      double totalE2 = sqrt(ts1.momentum().mag2() + piMassD0Squared) +
                       sqrt(ts2.momentum().mag2() + kaonMassD0Squared);
      double totalE2Sq = totalE2 * totalE2;

      double totalPSq =
          (ts1.momentum() + ts2.momentum()).mag2();

      auto sumMom = (ts1.momentum() + ts2.momentum());
      double totalPt = sumMom.perp();

      double mass1 = sqrt(totalE1Sq - totalPSq);
      double mass2 = sqrt(totalE2Sq - totalPSq);

      if ((mass1 > mPiKCutMax || mass1 < mPiKCutMin) && (mass2 > mPiKCutMax || mass2 < mPiKCutMin))
        continue;

      if (totalPt < dPtCut)
        continue;
      // the above keeps the pairs that align with k mass and pi mass and assigned pid

      // Create the vertex fitter object and vertex the tracks

      float posCandTotalE[2] = {0.0};
      float negCandTotalE[2] = {0.0};
      float d0TotalE[2] = {0.0};

      // create TransientTracks from the packedcandidates, need them for the kinematic fit
      reco::TransientTrack posTT(tk1.pseudoTrack(), magField);
      reco::TransientTrack negTT(tk2.pseudoTrack(), magField);

      // Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;
      float chi = 0.0;
      float ndf = 0.0;

      // loop over the two hypotheses, k-. pi+ and pi- k+
      for (int i = 0; i < 2; i++)
      {
        std::vector<RefCountedKinematicParticle> d0Particles;

        d0Particles.push_back(pFactory.particle(posTT, posCandMass[i], chi, ndf, posCandMass_sigma[i]));
        d0Particles.push_back(pFactory.particle(negTT, negCandMass[i], chi, ndf, negCandMass_sigma[i]));

        KinematicParticleVertexFitter d0Fitter;
        RefCountedKinematicTree d0Vertex;
        d0Vertex = d0Fitter.fit(d0Particles);

        if (!d0Vertex->isValid())
          continue;

        d0Vertex->movePointerToTheTop();
        RefCountedKinematicParticle d0Cand = d0Vertex->currentParticle();
        if (!d0Cand->currentState().isValid())
          continue;

        RefCountedKinematicVertex d0DecayVertex = d0Vertex->currentDecayVertex();
        if (!d0DecayVertex->vertexIsValid())
          continue;

        float d0C2Prob = TMath::Prob(d0DecayVertex->chiSquared(), d0DecayVertex->degreesOfFreedom());
        if (d0C2Prob < VtxChiProbCut)
          continue;

        d0Vertex->movePointerToTheFirstChild();
        RefCountedKinematicParticle posCand = d0Vertex->currentParticle();
        d0Vertex->movePointerToTheNextChild();
        RefCountedKinematicParticle negCand = d0Vertex->currentParticle();

        if (!posCand->currentState().isValid() || !negCand->currentState().isValid())
          continue;

        KinematicParameters posCandKP = posCand->currentState().kinematicParameters();
        KinematicParameters negCandKP = negCand->currentState().kinematicParameters();
        // Grab basic objects

        GlobalVector d0TotalP = GlobalVector(d0Cand->currentState().globalMomentum().x(),
                                             d0Cand->currentState().globalMomentum().y(),
                                             d0Cand->currentState().globalMomentum().z());

        GlobalVector posCandTotalP = GlobalVector(posCandKP.momentum().x(), posCandKP.momentum().y(), posCandKP.momentum().z());
        GlobalVector negCandTotalP = GlobalVector(negCandKP.momentum().x(), negCandKP.momentum().y(), negCandKP.momentum().z());

        posCandTotalE[i] = sqrt(posCandTotalP.mag2() + posCandMass[i] * posCandMass[i]);
        negCandTotalE[i] = sqrt(negCandTotalP.mag2() + negCandMass[i] * negCandMass[i]);
        d0TotalE[i] = posCandTotalE[i] + negCandTotalE[i];

        const Particle::LorentzVector d0P4(d0TotalP.x(), d0TotalP.y(), d0TotalP.z(), d0TotalE[i]);

        Particle::Point d0Vtx((*d0DecayVertex).position().x(), (*d0DecayVertex).position().y(), (*d0DecayVertex).position().z());

        std::vector<double> d0VtxEVec;
        d0VtxEVec.push_back(d0DecayVertex->error().cxx());
        d0VtxEVec.push_back(d0DecayVertex->error().cyx());
        d0VtxEVec.push_back(d0DecayVertex->error().cyy());
        d0VtxEVec.push_back(d0DecayVertex->error().czx());
        d0VtxEVec.push_back(d0DecayVertex->error().czy());
        d0VtxEVec.push_back(d0DecayVertex->error().czz());
        SMatrixSym3D d0VtxCovMatrix(d0VtxEVec.begin(), d0VtxEVec.end());
        const Vertex::CovarianceMatrix d0VtxCov(d0VtxCovMatrix);
        double d0VtxChi2(d0DecayVertex->chiSquared());
        double d0VtxNdof(d0DecayVertex->degreesOfFreedom());
        double d0NormalizedChi2 = d0VtxChi2 / d0VtxNdof;

        double rVtxMag = 99999.0;
        double lVtxMag = 99999.0;
        double sigmaRvtxMag = 999.0;
        double sigmaLvtxMag = 999.0;
        double d0Angle3D = -100.0;
        double d0Angle2D = -100.0;

        GlobalVector d0LineOfFlight = GlobalVector(d0Vtx.x() - xVtx,
                                                   d0Vtx.y() - yVtx,
                                                   d0Vtx.z() - zVtx);

        SMatrixSym3D d0TotalCov;
        if (isVtxPV)
          d0TotalCov = d0VtxCovMatrix + vtxPrimary->covariance();
        else
          d0TotalCov = d0VtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

        SVector3 distanceVector3D(d0LineOfFlight.x(), d0LineOfFlight.y(), d0LineOfFlight.z());
        SVector3 distanceVector2D(d0LineOfFlight.x(), d0LineOfFlight.y(), 0.0);

        d0Angle3D = angle(d0LineOfFlight.x(), d0LineOfFlight.y(), d0LineOfFlight.z(),
                          d0TotalP.x(), d0TotalP.y(), d0TotalP.z());
        d0Angle2D = angle(d0LineOfFlight.x(), d0LineOfFlight.y(), (float)0.0,
                          d0TotalP.x(), d0TotalP.y(), (float)0.0);

        lVtxMag = d0LineOfFlight.mag();
        rVtxMag = d0LineOfFlight.perp();
        sigmaLvtxMag = sqrt(ROOT::Math::Similarity(d0TotalCov, distanceVector3D)) / lVtxMag;
        sigmaRvtxMag = sqrt(ROOT::Math::Similarity(d0TotalCov, distanceVector2D)) / rVtxMag;

        if (d0NormalizedChi2 > chi2Cut ||
            rVtxMag < rVtxCut ||
            rVtxMag / sigmaRvtxMag < rVtxSigCut ||
            lVtxMag < lVtxCut ||
            lVtxMag / sigmaLvtxMag < lVtxSigCut ||
            cos(d0Angle3D) < collinCut3D || cos(d0Angle2D) < collinCut2D || d0Angle3D > alphaCut || d0Angle2D > alpha2DCut)
          continue;

        // // 3D IP wrt PV
        AnalyticalImpactPointExtrapolator extrap(magField);
        TrajectoryStateOnSurface tsos =
            extrap.extrapolate(d0Cand->currentState().freeTrajectoryState(),
                               RecoVertex::convertPos(vtxPrimary->position()));
        if (!tsos.isValid())
          continue;

        Measurement1D cur3DIP;
        VertexDistance3D a3d;
        GlobalPoint refPoint = tsos.globalPosition();
        GlobalError refPointErr = tsos.cartesianError().position();
        GlobalPoint vertexPosition = RecoVertex::convertPos(vtxPrimary->position());
        GlobalError vertexPositionErr = RecoVertex::convertError(vtxPrimary->error());
        cur3DIP = (a3d.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));

        // Two trkDCA

        FreeTrajectoryState posStateNew = posTT.impactPointTSCP().theState();
        FreeTrajectoryState negStateNew = negTT.impactPointTSCP().theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(posStateNew, negStateNew);
        if (!cApp.status())
          continue;
        float dca = fabs(cApp.distance());
        TwoTrackMinimumDistance minDistCalculator;
        minDistCalculator.calculate(state1, state2);
        dca = minDistCalculator.distance();
        cxPt = minDistCalculator.crossingPoint();
        GlobalError posErr = posStateNew.cartesianError().position();
        GlobalError negErr = negStateNew.cartesianError().position();

        // DCA error propagation
        double sigma_x2 = posErr.cxx() + negErr.cxx();
        double sigma_y2 = posErr.cyy() + negErr.cyy();

        double dcaError = sqrt(sigma_x2 * cxPt.x() * cxPt.x() + sigma_y2 * cxPt.y() * cxPt.y()) / dca;

        // Create CompositeCandidate
        pat::CompositeCandidate theD0;
        theD0.setP4(d0P4);
        theD0.setPdgId(pdg_id[i]);
        theD0.addUserFloat("track3DDCA", dca);
        theD0.addUserFloat("track3DDCAErr", dcaError);
        theD0.addUserFloat("ip3d", cur3DIP.value());
        theD0.addUserFloat("ip3derr", cur3DIP.error());

        // Add the two daughters

        tk1.setPdgId(pos_pdg_id[i]);
        tk2.setPdgId(neg_pdg_id[i]);
        theD0.addDaughter(tk1);
        theD0.addDaughter(tk2);

        // Vertex coordinates & covariance
        theD0.addUserFloat("vtxX", d0Vtx.x());
        theD0.addUserFloat("vtxY", d0Vtx.y());
        theD0.addUserFloat("vtxZ", d0Vtx.z());
        theD0.addUserFloat("bestvtxX", xVtx);
        theD0.addUserFloat("bestvtxY", yVtx);
        theD0.addUserFloat("bestvtxZ", zVtx);
        theD0.addUserFloat("xVtxError", xVtxError);
        theD0.addUserFloat("yVtxError", yVtxError);
        theD0.addUserFloat("zVtxError", zVtxError);
        theD0.addUserFloat("vtxChi2", d0VtxChi2);
        theD0.addUserFloat("vtxNdof", d0VtxNdof);
        theD0.addUserFloat("dzvtx", dzvtx);
        theD0.addUserFloat("dxyvtx", dxyvtx);
        theD0.addUserFloat("dzError", dzerror);
        theD0.addUserFloat("dxyError", dxyerror);

        for (int ii = 0; ii < 3; ++ii)
        {
          for (int jj = 0; jj < 3; ++jj)
          {
            theD0.addUserFloat("vertexCovariance_" + std::to_string(ii) + "_" + std::to_string(jj), d0VtxCov(ii, jj));
          }
        }
        if (fabs(theD0.mass() - d0MassD0) < d0MassCut)
        {
          theD0s.push_back(theD0);
        }
      }
    }
  }
}
// Get methods

const pat::CompositeCandidateCollection &D0Fitter::getD0() const
{
  return theD0s;
}

const std::vector<float> &D0Fitter::getMVAVals() const
{
  return mvaVals_;
}

void D0Fitter::resetAll()
{
  theD0s.clear();
  mvaVals_.clear();
}
