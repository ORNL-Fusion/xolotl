// Includes
#include "HeVCluster.h"
#include <iostream>

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
	
	// xHe * yV + zHe --> (x + z)He + yV
	for (int z = 1; numHe + z <= maxHeClusterSize; z++) {
		// Select zHe
		connectivityArray.at(z - 1) = 1;
	}
	
	// xHe * yV + V   --> xHe + (y + 1)V
	if (numV + 1 <= maxVClusterSize) {
		// Select V
		connectivityArray.at(numHeClusters) = 1;
	}
	
	// xHe * yV + zI  --> xHe + (y - z)V
	for (int z = 1; numV - z >= 1; z++) {
		// Select zI
		connectivityArray.at(numHeClusters + numVClusters + z - 1) = 1;
	}
	
	
	return connectivityArray;
}

double HeVCluster::getDissociationFlux(const double temperature) {

	// Local Declarations
	int vIndex = -1, iIndex = -1, heIndex = -1;
	int thisIndex = network->toClusterIndex(getClusterMap());
	std::map<std::string, int> oneHe;
	std::map<std::string, int> oneV;
	std::map<std::string, int> oneI;

	// Set the cluster map data for 1 of each species
	oneHe["He"] = 1; oneHe["V"] = 0; oneHe["I"] = 0;
	oneV["He"] = 0; oneV["V"] = 1; oneV["I"] = 0;
	oneI["He"] = 0; oneI["V"] = 0; oneI["I"] = 1;

	// Get their indices in the array
	int oneHeIndex = network->toClusterIndex(oneHe);
	int oneVIndex = network->toClusterIndex(oneV);
	int oneIIndex = network->toClusterIndex(oneI);

	// Find the indices such that they are the index of the
	// concentration in C_bar that contain one additional species
	// than this reactant
	for (int k = 0; k < network->reactants->size(); k++) {
		if ((network->toClusterMap(k)["He"] - numHe) == 1 && heIndex == -1) {
			heIndex = k;
		}
		if ((network->toClusterMap(k)["V"] - numV) == 1 && vIndex == -1) {
			vIndex = k;
		}
		if (network->toClusterMap(k)["I"] == 1 && iIndex == -1) {
			iIndex = k;
		}
	}

	// Calculate and return the dissociation constant
	return calculateDissociationConstant(vIndex, oneVIndex, temperature)
			* network->reactants->at(vIndex)->getConcentration()
			+ calculateDissociationConstant(iIndex, oneIIndex, temperature)
					* network->reactants->at(iIndex)->getConcentration()
			+ calculateDissociationConstant(heIndex, oneHeIndex, temperature)
					* network->reactants->at(heIndex)->getConcentration()
			- (calculateDissociationConstant(thisIndex, oneVIndex, temperature)
					+ calculateDissociationConstant(thisIndex, oneIIndex, temperature) +
					calculateDissociationConstant(thisIndex, oneHeIndex, temperature)) * getConcentration();
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
