//+++++++++++++++++++++++++++++++++++                                                                                                                    
//Code owner: Nihar Ranjan saha
//contact mail: nihar.ranjan.saha@cern.ch                                                                                                                
//Date: 24 Sept 2025                                                                                                                                     
//+++++++++++++++++++++++++++++++++++
/**\class LamC3PProducer LamC3PProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/LamC3PProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/


// system include files
#include <memory>

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PProducer.h"

// Constructor
LamC3PProducer::LamC3PProducer(const edm::ParameterSet& iConfig) :
theCandidates(iConfig, consumesCollector(), selectedTkhidxSetVec)
{
  useAnyMVA_ = false;
  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< std::vector<pat::CompositeCandidate> >("LamC3P");
  if(useAnyMVA_) produces<MVACollection>("MVAValuesLamC3P");
}



// (Empty) Destructor
LamC3PProducer::~LamC3PProducer() {
}



// Producer Method
void LamC3PProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   theCandidates.fitAll(iEvent, iSetup);

   auto LamC3PCandidates = std::make_unique<pat::CompositeCandidateCollection>();
   LamC3PCandidates->reserve( theCandidates.getLamC3P().size() );

   std::copy( theCandidates.getLamC3P().begin(),
              theCandidates.getLamC3P().end(),
              std::back_inserter(*LamC3PCandidates) );

   // Write the collections to the Event
   iEvent.put( std::move(LamC3PCandidates), std::string("LamC3P") );
    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theCandidates.getMVAVals().begin(),theCandidates.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesLamC3P"));
   }

   theCandidates.resetAll();
}


void LamC3PProducer::beginJob() {
}


void LamC3PProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(LamC3PProducer);
