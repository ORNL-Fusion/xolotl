// Includes
#include "HeVCluster.h"
#include <iostream>
#include <Constants.h>

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

std::vector<int> HeVCluster::getReactionConnectivity() {
	
	// Extract some of the properties from the network
	
	std::shared_ptr<std::map<std::string, std::string>> properties =
		network->properties;
	
	int maxMixedClusterSize = std::stoi(properties->at("maxMixedClusterSize"));
	
	// Initialize the return array with zeroes
	std::vector<int> connectivityArray(network->reactants->size(), 0);
	
	// This cluster is involved in the following interactions:
	
	// xHe*yV + zHe --> (x + z)He*yV
	for (int z = 1; numHe + numV + z <= maxMixedClusterSize; z++) {
		// Select the zHe index
		std::map<std::string, int> speciesMap;
		speciesMap["He"] = z;
		int i = network->toClusterIndex(speciesMap);
		
		connectivityArray.at(i) = 1;
	}
	
	// xHe*yV + V   --> xHe*(y + 1)V
	if (numHe + numV + 1 <= maxMixedClusterSize) {
		// Select the single V index
		std::map<std::string, int> speciesMap;
		speciesMap["V"] = 1;
		int i = network->toClusterIndex(speciesMap);
		
		connectivityArray.at(i) = 1;
	}
	
	// xHe*yV + zI  --> xHe*(y - z)V
	
	// Here I am assuming that the HeV and Interstitial can only interact if
	// they would produce a positive number of vacancy species
	
	for (int numIOther = 1; numV - numIOther >= 1; numIOther++) {
		// Select the zI index
		std::map<std::string, int> speciesMap;
		speciesMap["I"] = numIOther;
		int i = network->toClusterIndex(speciesMap);
		
		connectivityArray.at(i) = 1;
	}
	
	// Everything else is 0 (not connected)
	
	return connectivityArray;
}

std::vector<int> HeVCluster::getDissociationConnectivity() {
	// Local Declarations
	int nReactants = network->reactants->size();
	std::vector<int> dissConnections(nReactants, 0);
	std::map<std::string, int> clusterMap;

	// Vacancy Dissociation
	clusterMap["He"] = numHe-1; clusterMap["V"] = numV; clusterMap["I"] = 0;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["V"] = 0; clusterMap["He"] = 1;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;

	// Trap Mutation
	clusterMap["V"] = numV + 1; clusterMap["He"] = numHe;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["I"] = 1; clusterMap["V"] = 0; clusterMap["He"] = 0;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;

	// Vacancy Dissociation
	clusterMap["He"] = numHe; clusterMap["V"] = numV - 1; clusterMap["I"] = 0;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["He"] = 0; clusterMap["V"] = 1; clusterMap["I"] = 0;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;

	// Return the connections
	return dissConnections;
}

bool HeVCluster::isProductReactant(int reactantI, int reactantJ) {
	// Local Declarations, integers for species number for I, J reactants
	int rI_I = 0, rJ_I = 0, rI_He = 0, rJ_He = 0, rI_V = 0, rJ_V = 0;

	// Get the ClusterMap corresponding to
	// the given reactants
	std::map<std::string, int> reactantIMap = network->toClusterMap(reactantI);
	std::map<std::string, int> reactantJMap = network->toClusterMap(reactantJ);

	// Grab the numbers for each species
	// from each Reactant
	rI_I = reactantIMap["I"]; rJ_I = reactantJMap["I"];
	rI_He = reactantIMap["He"]; rJ_He = reactantJMap["He"];
	rI_V = reactantIMap["V"]; rJ_V = reactantJMap["V"];

	// We should have no interstitials, a
	// total of numHe Helium, and a total of
	// numV Vacancies
	return ((rI_I + rJ_I) == 0)
			&& ((rI_He + rJ_He) == numHe)
			&& ((rI_V + rJ_V) == numV);
}


std::map<std::string, int> HeVCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = numHe;
	clusterMap["V"] = numV;
	clusterMap["I"] = 0;

	// Return it
	return clusterMap;
}

double HeVCluster::getReactionRadius() {
	return (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
				+ pow((3.0 * pow(xolotlCore::latticeConstant, 3.0) * numV)
								/ (8.0 * xolotlCore::pi),
								(1.0 / 3.0))
				- pow((3.0 * pow(xolotlCore::latticeConstant, 3.0))
								/ (8.0 * xolotlCore::pi),
								(1.0 / 3.0));
}
