// Includes
#include "HeVCluster.h"

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV) :
	PSICluster(1), numHe(numHe), numV(numV) {
	
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Set the reactant name appropriately
	name = "HeV Cluster";
}

HeVCluster::~HeVCluster() {}

double HeVCluster::getGenByEm() {
	return 0;
}

double HeVCluster::getAnnByEm() {
	return 0;
}

int HeVCluster::getSpeciesSize(const std::string speciesName) {
	if (speciesName == "He") {
		return numHe;
	}
	else if (speciesName == "V") {
		return numV;
	}
	else {
		return 0;
	}
}

std::vector<int> HeVCluster::getConnectivity() {
	
	// Extract some of the properties from the network
	
	std::shared_ptr<std::map<std::string, std::string>> properties =
		network->properties;
	
	int maxHeClusterSize = std::stoi(properties->at("maxHeClusterSize"));
	int maxVClusterSize = std::stoi(properties->at("maxVClusterSize"));
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	int numVClusters = std::stoi(properties->at("numVClusters"));
	int numIClusters = std::stoi(properties->at("numIClusters"));
	
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;
	
	// Initialize the return array with zeroes
	std::vector<int> connectivityArray(network->reactants->size(), 0);
	
	// This cluster's index in the reactants array
	int clusterIndex = numSingleSpeciesClusters + (numHe - 1) + (numV - 1) * numHeClusters;
	
	
	// This cluster is involved in the following interactions:
	
	// HeV[x, y] + He[z] --> HeV[x + z, y]
	for (int z = 1; numHe + z <= maxHeClusterSize; z++) {
		// Select He[z]
		connectivityArray.at(z - 1) = 1;
	}
	
	// HeV[x, y] + V[1]  --> HeV[x, y + 1]
	if (numV + 1 <= maxVClusterSize) {
		// Select V[1]
		connectivityArray.at(numHeClusters) = 1;
	}
	
	// HeV[x, y] + I[z]  --> HeV[x, y - z]
	for (int z = 1; numV - z >= 1; z++) {
		// Select I[z]
		connectivityArray.at(numHeClusters + numVClusters + z - 1) = 1;
	}
	
	
	return connectivityArray;
}

double HeVCluster::getDissociationFlux(const double temperature) {
	return 0.0;
}

double HeVCluster::getProductionFlux(const double temperature) {
	// Local declarations
	double fluxOne = 0.0, fluxTwo = 0.0, kPlus = 0.0;
	int thisClusterIndex = 0;
	int numHeClusters = std::stoi(network->properties->at("numHeClusters"));
	int numVClusters = std::stoi(network->properties->at("numVClusters"));
	int numIClusters = std::stoi(network->properties->at("numIClusters"));
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;

	// This cluster's index in the reactants array
	thisClusterIndex = numSingleSpeciesClusters + (numHe - 1) + (numV - 1) * numHeClusters;

	// Loop over all possible clusters
	for (int j = 0; j < network->reactants->size(); j++) {
		for (int k = 0; k < network->reactants->size(); k++) {
			// If the jth and kth reactants react to produce this reactant...
			if (network->isConnected(j, k)
					&& (network->getReactionProduct(j, k).get() == this)) {

				// This fluxOne term considers all reactions that
				// produce C_i
				fluxOne = fluxOne + calculateReactionRateConstant(j, k, temperature)
								* network->reactants->at(j)->getConcentration()
								* network->reactants->at(k)->getConcentration();
			}
		}

		// Calculate Second term of production flux
		// this acts to take away from the current reactant
		// as it is reacting with others, thus decreasing itself.
		// This considers all populations that are produced by C_i
		if (network->isConnected(j, thisClusterIndex)) {
			fluxTwo = fluxTwo
					+ calculateReactionRateConstant(thisClusterIndex, j, temperature)
							* network->reactants->at(j)->getConcentration();
		}
	}

	// Return the production flux
	return fluxOne - (fluxTwo * getConcentration());

}
