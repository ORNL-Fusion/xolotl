// Includes
#include "Reactant.h"
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

double Reactant::getTotalFlux() {
	//TODO Auto-generated method stub
	return 0;
}

void Reactant::setReactionNetwork(const std::shared_ptr<ReactionNetwork> reactionNetwork) {
	network = reactionNetwork;
	std::cout << "Number of Reactants = " << network->reactants->size() << std::endl;
	return;
}

std::vector<int> Reactant::getConnectivity() {
	std::vector<int> dummy;
	return dummy;
}

const std::string Reactant::toString() {
	//TODO Auto-generated method stub
	return "fluffy";
}

ReactionNetwork::ReactionNetwork(const ReactionNetwork &other) {
	// The copy constructor of std::map copies each of the keys and values.
	properties.reset(new std::map<std::string, std::string>(*other.properties));

	// Copy the reactants list by pushing a copy of each element onto the
	// reactants vector
	reactants.reset(new std::vector<std::shared_ptr<Reactant>>);

	for (std::vector<std::shared_ptr<Reactant>>::iterator reactantIt =
			other.reactants->begin(); reactantIt != other.reactants->end();
			reactantIt++) {

		std::shared_ptr<Reactant> reactant(new Reactant(**reactantIt));
		reactants->push_back(reactant);
	}
}
