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
	std::vector<int> dummy;
	return dummy;
}

const std::string Reactant::toString() {
	throw std::string("Reactant::toString() is unimplemented");
}
