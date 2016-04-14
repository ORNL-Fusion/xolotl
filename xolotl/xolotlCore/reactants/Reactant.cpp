// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>

// Namespaces
using namespace xolotlCore;

Reactant::Reactant() :
		concentration(0.0),
		name(""), id(0),
		temperature(0.0),
		network(nullptr),
		handlerRegistry(nullptr) {}

Reactant::Reactant(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		concentration(0.0),
		name(""),
		id(0),
		temperature(0.0),
		network(nullptr),
		handlerRegistry(registry) {}

Reactant::Reactant(const Reactant &other) :
		concentration(other.concentration),
		name(other.name),
		typeName(other.typeName),
		id(other.id),
		temperature(other.temperature),
		network(other.network),
		compositionMap(other.compositionMap),
		handlerRegistry(other.handlerRegistry) {}

std::shared_ptr<Reactant> Reactant::clone() {
	std::shared_ptr<Reactant> reactant(new Reactant(*this));
	return reactant;
}

Reactant::Reactant(double conc, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		concentration(conc),
		name(""),
		id(0),
		temperature(0.0),
		network(nullptr),
		handlerRegistry(registry) {}

double Reactant::getConcentration() const {
	return concentration;
}

void Reactant::increaseConcentration(double deltaConc) {
	concentration += deltaConc;

	return;
}

void Reactant::decreaseConcentration(double deltaConc) {
	concentration -= deltaConc;

	return;
}

void Reactant::setConcentration(double conc) {
	concentration = conc;

	return;
}

void Reactant::zero() {
	concentration = 0.0;

	return;
}

double Reactant::getTotalFlux() {
	return 0.0;
}

void Reactant::setReactionNetwork(
	std::shared_ptr<ReactionNetwork> reactionNetwork) {
	network = reactionNetwork;
	
	return;
}


std::vector<int> Reactant::getConnectivity() const {
	// The connectivity array by default is filled with
	// zeros.
	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	return connectivity;
}

std::vector<double> Reactant::getPartialDerivatives() const {
	// The partial derivatives array by default is filled with
	// zeros.
	int length = network->size();
	std::vector<double> partial_derivatives = std::vector<double>(length, 0.0);

	return partial_derivatives;
}

void Reactant::getPartialDerivatives(std::vector<double> &) const {
	// Do nothing.
	return;
}

const std::string Reactant::getName() const {
	return name;
}

std::string Reactant::getType() const {
	return typeName;
}

const std::map<std::string, int> & Reactant::getComposition() const {
	// Return the composition map
	return compositionMap;
}

namespace xolotlCore {
std::ostream& operator<< (std::ostream& out, const Reactant& reactant) {
	// Output the name and Id of the reactant
	out << reactant.name << "  Reactant Id: " << reactant.id-1;

	return out;
}

}
