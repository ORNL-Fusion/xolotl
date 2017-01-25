#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <cassert>

using namespace xolotlCore;

ReactionNetwork::ReactionNetwork() :
		temperature(0.0), networkSize(0), reactionsEnabled(true), dissociationsEnabled(
				true), numVClusters(0), numIClusters(0), numSuperClusters(0), maxVClusterSize(0), maxIClusterSize(0) {
//    concUpdateCounter = xolotlPerf::getHandlerRegistry()->getEventCounter("net_conc_updates");
	// Setup the vector to hold all of the reactants
	allReactants = make_shared<std::vector<IReactant *>>();
	return;
}

ReactionNetwork::ReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		handlerRegistry(registry), temperature(0.0), networkSize(0), reactionsEnabled(
				true), dissociationsEnabled(true), numVClusters(0), numIClusters(0), numSuperClusters(0), maxVClusterSize(0), maxIClusterSize(0) {
	// Counter for the number of times the network concentration is updated.
	concUpdateCounter = handlerRegistry->getEventCounter("net_conc_updates");
	// Setup the vector to hold all of the reactants
	allReactants = make_shared<std::vector<IReactant *>>();

	return;
}

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	handlerRegistry = other.handlerRegistry;
	allReactants = other.allReactants;
	temperature = other.temperature;
	networkSize = other.networkSize;
	names = other.names;
	compoundNames = other.compoundNames;
	reactionsEnabled = other.reactionsEnabled;
	dissociationsEnabled = other.dissociationsEnabled;
	numVClusters = other.numVClusters;
	numIClusters = other.numIClusters;
	numSuperClusters = other.numSuperClusters;
	maxVClusterSize = other.maxVClusterSize;
	maxIClusterSize = other.maxIClusterSize;

	// TODO - do we copy the source ReactionNetwork's counter also?
	// Or should we have our own counter?  How to distinguish them by name?

	// Counter for the number of times the network concentration is updated.
	concUpdateCounter = handlerRegistry->getEventCounter("net_conc_updates");

	return;
}

void ReactionNetwork::fillConcentrationsArray(double * concentrations) {
	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();
	int id = -1;

	// Fill the array
	for (int i = 0; i < size; i++) {
		id = reactants->at(i)->getId() - 1;
		concentrations[id] = reactants->at(i)->getConcentration();
	}

	return;
}

void ReactionNetwork::updateConcentrationsFromArray(double * concentrations) {
	// Local Declarations
	auto allReactants = this->getAll();
	int id = -1;

	// Set the concentrations
	concUpdateCounter->increment();	// increment the update concentration counter
	for (auto iter = allReactants->begin(); iter != allReactants->end();
			++iter) {
		id = (*iter)->getId() - 1;
		(*iter)->setConcentration(concentrations[id]);
	}

	return;
}

void ReactionNetwork::askReactantsToReleaseNetwork(void) {
	// Get all the reactants
	auto allReactants = this->getAll();

	// Loop on each reactant to release the network
	for (auto iter = allReactants->begin(); iter != allReactants->end();
			++iter) {
		IReactant* currReactant = *iter;
		assert(currReactant != NULL);

		currReactant->releaseReactionNetwork();
	}
}

void ReactionNetwork::setTemperature(double temp) {
	// Set the temperature
	temperature = temp;

	// Update the temperature for all of the clusters
	for (int i = 0; i < networkSize; i++) {
		// This part will set the temperature in each reactant
		// and recompute the diffusion coefficient
		allReactants->at(i)->setTemperature(temp);
	}

	return;
}

double ReactionNetwork::getTemperature() const {
	return temperature;
}

IReactant * ReactionNetwork::get(const std::string& type,
		const int size) const {
	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	return (IReactant *) retReactant.get();
}

IReactant * ReactionNetwork::getCompound(const std::string& type,
		const std::vector<int>& sizes) const {
	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	return (IReactant *) retReactant.get();
}

const std::shared_ptr<std::vector<IReactant *>> & ReactionNetwork::getAll() const {
	return allReactants;
}

std::vector<IReactant *> ReactionNetwork::getAll(
		const std::string& name) const {
	// Local Declarations
	std::vector<IReactant *> reactants;

	return reactants;
}

const std::vector<std::string> & ReactionNetwork::getNames() const {
	return names;
}

const std::vector<std::string> & ReactionNetwork::getCompoundNames() const {
	return compoundNames;
}
