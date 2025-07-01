// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      D0Producer
//
/**\class D0Producer D0Producer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/D0Producer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//

// system include files
#include <memory>

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/D0Producer.h"

// Constructor
D0Producer::D0Producer(const edm::ParameterSet &iConfig) : theCandidates(iConfig, consumesCollector())
{
   produces<pat::CompositeCandidateCollection>("D0");
}

// (Empty) Destructor
D0Producer::~D0Producer()
{
}

//
// Methods
//

// Producer Method
void D0Producer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   using namespace edm;

   // Create D0Fitter object which reconstructs the vertices and creates
   theCandidates.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
   auto d0Candidates = std::make_unique<pat::CompositeCandidateCollection>();
   d0Candidates->reserve(theCandidates.getD0().size());

   std::copy(theCandidates.getD0().begin(),
             theCandidates.getD0().end(),
             std::back_inserter(*d0Candidates));

   // Write the collections to the Event
   iEvent.put(std::move(d0Candidates), std::string("D0"));

   theCandidates.resetAll();
}

// void D0Producer::beginJob() {
void D0Producer::beginJob()
{
}

void D0Producer::endJob()
{
}

// define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(D0Producer);
