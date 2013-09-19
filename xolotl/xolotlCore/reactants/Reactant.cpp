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

std::shared_ptr<Reactant> Reactant::clone() {
	std::shared_ptr<Reactant> reactant(new Reactant(*this));
	return reactant;
}

Reactant::Reactant(double conc) :
		concentration(0.0), name(""), id(0) {
}

double Reactant::getConcentration() const {
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
	std::shared_ptr<ReactionNetwork> reactionNetwork) {
	
	network = reactionNetwork;
}


std::vector<int> Reactant::getConnectivity() {
	// The connectivity array by default is filled with
	// zeroes.
	
	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	return connectivity;
}

/**
 * This operation returns the name of the reactant.
 * @return the name
 */
const std::string Reactant::getName() {
	return name;
}

const std::string Reactant::toString() {
	throw std::string("Reactant::toString() is unimplemented");
}
