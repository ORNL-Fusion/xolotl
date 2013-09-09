#include "ReactionNetwork.h"
#include "Reactant.h"

using namespace xolotlCore;

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	// The copy constructor of std::map copies each of the keys and values.
	properties.reset(new std::map<std::string, std::string>(*other.properties));
}

std::map<std::string, int> ReactionNetwork::toClusterMap(int index) const {
	// This base class returns an empty map
	return std::map<std::string, int>();
}

int ReactionNetwork::toClusterIndex(
		std::map<std::string, int> clusterMap) const {
	// This base class returns a zero index value
	return 0;
}

std::map<std::string, int> Reactant::getComposition() const {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Return an empty map
	return clusterMap;
}

/**
 * This operation fills an array of doubles with the concentrations of all
 * of the Reactants in the network.
 * @param concentrations The array that will be filled with the
 * concentrations. This operation does NOT create, destroy or resize the
 * array. If the array is to small to hold the concentrations, SIGSEGV will
 * be thrown.
 */
void ReactionNetwork::fillConcentrationsArray(double * concentrations) {

	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();

	// Fill the array
	for (int i = 0; i < size; i++) {
		concentrations[i] = reactants->at(i)->getConcentration();
	}

	return;
}

/**
 * This operation updates the concentrations for all Reactants in the
 * network from an array.
 * @param concentrations The array of doubles that will be for the
 * concentrations. This operation does NOT create, destroy or resize the
 * array. Properly aligning the array in memory so that this operation
 * does not overrun is up to the caller.
 */
void ReactionNetwork::updateConcentrationsFromArray(double * concentrations) {

	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();

	// Set the concentrations
	for (int i = 0; i < size; i++) {
		reactants->at(i)->setConcentration(concentrations[i]);
	}

	return;
}
