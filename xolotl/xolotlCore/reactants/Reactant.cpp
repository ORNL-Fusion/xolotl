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


std::map<std::string, int> ReactionNetwork::toClusterMap(int index) {
	// This method defines the ordering and placement of all
	// reactions in a ReactionNetwork
	
	std::map<std::string, int> clusterMap;
	
	// He clusters
	// size of `numHeClusters`
	
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	
	if (index < numHeClusters) {
		clusterMap["He"] = index + 1;
		return clusterMap;
	}
	
	index -= numHeClusters;
	
	// V (vacancy) clusters
	// size of `numVClusters`
	
	int numVClusters = std::stoi(properties->at("numVClusters"));
	
	if (index < numVClusters) {
		clusterMap["V"] = index + 1;
		return clusterMap;
	}
	
	index -= numVClusters;
	
	// I (interstitial) clusters
	// size of `numIClusters`
	
	int numIClusters = std::stoi(properties->at("numIClusters"));
	
	if (index < numIClusters) {
		clusterMap["I"] = index + 1;
		return clusterMap;
	}
	
	index -= numIClusters;
	
	// HeV clusters
	// size of `numHeVClusters`
	
	int numHeVClusters = std::stoi(properties->at("numHeVClusters"));
	int maxMixedClusterSize = std::stoi(properties->at("maxMixedClusterSize"));
	
	// Iterate through all posibilities of quantities for He and V
	// under the condition that
	// numHe + numV <= maxMixedClusterSize.
	
	// This doesn't need to be a loop, but it's less error prone than
	// a closed form.
	
	for (int numV = 1; numV <= maxMixedClusterSize; numV++) {
		int maxHe = maxMixedClusterSize - numV;
		int numHe = index + 1;
		
		if (numHe <= maxHe) {
			// The index must be referring to this quantity of He and V.
			// Finalize and return the cluster map
			
			clusterMap["He"] = numHe;
			clusterMap["V"] = numV;
			return clusterMap;
		}
		
		// Decrease the index by the number of vacancies
		// possible with the current quantity of helium
		index -= maxHe;
	}
	
	// HeI clusters
	// TODO
	
	throw std::string("HeI clusters cannot yet be indexed");
}


int ReactionNetwork::toClusterIndex(std::map<std::string, int> clusterMap) {
	throw std::string("Unimplemented: ReactionNetwork::toClusterIndex()");
}
