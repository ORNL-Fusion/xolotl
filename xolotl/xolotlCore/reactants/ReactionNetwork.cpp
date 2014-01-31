#include "ReactionNetwork.h"
#include "Reactant.h"
#include <iostream>

using namespace xolotlCore;

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	// The copy constructor of std::map copies each of the keys and values.
	properties.reset(new std::map<std::string, std::string>(*other.properties));
}

/**
 * This operation fills an array of doubles with the concentrations of all
 * of the Reactants in the network. The array is ordered by increasing
 * reactant id, starting at 1.
 * @param concentrations The array that will be filled with the
 * concentrations. This operation does NOT create, destroy or resize the
 * array. If the array is to small to hold the concentrations, SIGSEGV will
 * be thrown.
 */
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
	int id = 1;

	// Set the concentrations
	for (int i = 0; i < size; i++) {
		id = reactants->at(i)->getId() - 1;
		reactants->at(i)->setConcentration(concentrations[id]);
	}

	return;
}
