// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include "xolotlPerf/HandlerRegistryFactory.h"
#include <iostream>

// Namespaces
using namespace xolotlCore;

Reactant::Reactant() :
		concentration(0.0e-16), name(""), id(0) {

	// Temporary definition used to specify the type of performance handler registry
	// that will be used (standard or dummy).  Currently, the handler registry is set via
	// commandline argument and the only way to access which handler registry will be used
	// is through getHandlerRegistry
//	handlerRegistry = xolotlPerf::getHandlerRegistry();
}

Reactant::Reactant(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		concentration(0.0e-16), name(""), id(0), handlerRegistry(registry) {

}

Reactant::Reactant(const Reactant &other) :
		concentration(other.concentration), name(other.name), id(other.id),
		compositionMap(other.compositionMap), handlerRegistry(other.handlerRegistry) {
}

Reactant::~Reactant() {
}

std::shared_ptr<Reactant> Reactant::clone() {
	std::shared_ptr<Reactant> reactant(new Reactant(*this));
	return reactant;
}

Reactant::Reactant(double conc, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		concentration(conc), name(""), id(0), handlerRegistry(registry) {
}

double Reactant::getConcentration() const {
	return concentration;
}

void Reactant::increaseConcentration(double deltaConc) {
	concentration += deltaConc;
	std::cout << "Conc = " << concentration << ", delta = " << deltaConc << std::endl;
	return;
}

void Reactant::decreaseConcentration(double deltaConc) {
	concentration -= deltaConc;
	return;
}

void Reactant::setConcentration(double conc) {
	concentration = std::max(conc, 0.);
	return;
}

void Reactant::zero() {
	concentration = 0.0e-16;
	return;
}

double Reactant::getTotalFlux(const double temperature) {
	return 0.0;
}

void Reactant::setReactionNetwork(
	std::shared_ptr<ReactionNetwork> reactionNetwork) {
	
	network = reactionNetwork;
}


std::vector<int> Reactant::getConnectivity() const {
	// The connectivity array by default is filled with
	// zeroes.
	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	return connectivity;
}


/**
 * This operation returns the list of partial derivatives of this Reactant
 * with respect to all other reactants in the network. The combined lists
 * of partial derivatives from all of the reactants in the network can be
 * used to form, for example, a Jacobian.
 *
 * @param the temperature at which the reactions are occurring
 * @return The partial derivatives for this reactant where index zero
 * corresponds to the first reactant in the list returned by the
 * ReactionNetwork::getAll() operation.
 */
std::vector<double> Reactant::getPartialDerivatives(double temperature) const {
	// The partial derivatives array by default is filled with
	// zeroes.
	int length = network->size();
	std::vector<double> partial_derivatives = std::vector<double>(length, 0.0);

	return partial_derivatives;
}

/**
 * This operation returns the name of the reactant.
 * @return the name
 */
const std::string Reactant::getName() const {
	return name;
}

const std::string Reactant::toString() {
	throw std::string("Reactant::toString() is unimplemented");
}

const std::map<std::string, int> & Reactant::getComposition() const {
	// Return the composition map
	return compositionMap;
}
