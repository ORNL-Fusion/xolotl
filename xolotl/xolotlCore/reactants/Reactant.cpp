// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include <iostream>

// Namespaces
using namespace xolotlCore;

Reactant::Reactant() :
		concentration(0.0), name(""), id(0) {
}

Reactant::Reactant(const Reactant &other) :
		concentration(other.concentration), name(other.name), id(other.id) {
}

Reactant::~Reactant() {
}

Reactant::Reactant(double conc) :
		concentration(0.0), name(""), id(0) {
}

double Reactant::getConcentration() {
	return concentration;
}

void Reactant::increaseConcentration(double deltaConc) {
	concentration += deltaConc;
}

void Reactant::decreaseConcentration(double deltaConc) {
	concentration -= deltaConc;
}

void Reactant::setConcentration(double conc) {
	concentration = conc;
}

void Reactant::zero() {
	concentration = 0.0;
}

double Reactant::getTotalFlux(const double temperature) {
	return 0;
}

void Reactant::setReactionNetwork(
		const std::shared_ptr<ReactionNetwork> reactionNetwork) {
	network = reactionNetwork;
	return;
}

std::vector<int> Reactant::getConnectivity() {
	
	std::vector<int> reactionConn = getReactionConnectivity();
	std::vector<int> dissConn = getDissociationConnectivity();
	
	// The vectors must be the same length.
	
	if (reactionConn.size() != dissConn.size()) {
		throw std::string("The reaction and dissociation vectors "
			"must be the same length");
	}
	
	// Merge the two vectors such that the final vector contains
	// a 1 at a positioin if either of the connectivity arrays
	// have a 1
	
	int connLength = reactionConn.size();
	for (int i = 0; i < connLength; i++) {
		
		// We can modify the reaction vector in place since it
		// will not be needed by this method again.
		reactionConn[i] |= dissConn[i];
	}
	
	return reactionConn;
}

const std::string Reactant::toString() {
	throw std::string("Reactant::toString() is unimplemented");
}


std::vector<int> Reactant::getReactionConnectivity()
{
	// By default, return an array with a zero for each reactant
	// in the network
	
	std::vector<int> reactionConn(network->reactants->size(), 0);
	return reactionConn;
}


std::vector<int> Reactant::getDissociationConnectivity()
{
	// By default, return an array with a zero for each reactant
	// in the network
	
	std::vector<int> dissConn(network->reactants->size(), 0);
	return dissConn;
}
