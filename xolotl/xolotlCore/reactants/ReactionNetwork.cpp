#include "ReactionNetwork.h"
#include "Reactant.h"
#include <xolotlPerf.h>
#include <iostream>
#include <cassert>

using namespace xolotlCore;

ReactionNetwork::ReactionNetwork()
  : properties(new std::map<std::string, std::string>())
{
//    concUpdateCounter = xolotlPerf::getHandlerRegistry()->getEventCounter("net_conc_updates");
	return;
}

ReactionNetwork::ReactionNetwork(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry)
  : properties(new std::map<std::string, std::string>()), handlerRegistry(registry)
{
	// Counter for the number of times the network concentration is updated.
    concUpdateCounter = handlerRegistry->getEventCounter("net_conc_updates");

	return;
}

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	// The copy constructor of std::map copies each of the keys and values.
	properties.reset(new std::map<std::string, std::string>(*other.properties));

	handlerRegistry = other.handlerRegistry;

    // TODO - do we copy the source ReactionNetwork's counter also?
    // Or should we have our own counter?  How to distinguish them by name?
//    concUpdateCounter = xolotlPerf::getHandlerRegistry()->getEventCounter("net_conc_updates");

	// Counter for the number of times the network concentration is updated.
    concUpdateCounter = handlerRegistry->getEventCounter("net_conc_updates");

	return;
}

void ReactionNetwork::fillConcentrationsArray(double * concentrations) {
	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();
	int id = 1;

	// Fill the array
	for (int i = 0; i < size; i++) {
		id = reactants->at(i)->getId() - 1;
		concentrations[id] = reactants->at(i)->getConcentration();
	}

	return;
}

void ReactionNetwork::updateConcentrationsFromArray(double * concentrations) {
	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();
	int id = 1;

	// Set the concentrations
    concUpdateCounter->increment();		// increment the update concentration counter
	for (int i = 0; i < size; i++) {
		id = reactants->at(i)->getId() - 1;
		reactants->at(i)->setConcentration(concentrations[id]);
	}

	return;
}

void ReactionNetwork::askReactantsToReleaseNetwork(void) {
	// Get all the reactants
    auto allReactants = this->getAll();    

    // Loop on each reactant to release the network
    for(auto iter = allReactants->begin(); iter != allReactants->end(); ++iter) {
        Reactant* currReactant = *iter;
        assert( currReactant != NULL );

        currReactant->releaseReactionNetwork();
    }
}
