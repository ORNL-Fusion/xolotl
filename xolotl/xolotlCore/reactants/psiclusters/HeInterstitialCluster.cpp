// Includes
#include "HeInterstitialCluster.h"

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(std::map<std::string, int> speciesMap) :
		PSICluster(1) {
	numHe = speciesMap["He"];
	numI = speciesMap["I"];
	size = numHe + numI;
	// Set the reactant name appropriately
	name = "He-Interstitial Cluster";
}
HeInterstitialCluster::~HeInterstitialCluster() {
	//TODO Auto-generated method stub
}
double HeInterstitialCluster::getGenByEm() {
	//TODO Auto-generated method stub
	return 0;
}
double HeInterstitialCluster::getAnnByEm() {
	//TODO Auto-generated method stub
	return 0;
}
int HeInterstitialCluster::getSpeciesSize(const std::string speciesName) {
	//TODO Auto-generated method stub
	return 0;
}


double HeInterstitialCluster::getDissociationFlux(const double temperature) {

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
		if (network->toClusterMap(k)["V"] == 1 && vIndex == -1) {
			vIndex = k;
		}
		if ((network->toClusterMap(k)["I"] - numI) == 1 && iIndex == -1) {
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
					+ calculateDissociationConstant(thisIndex, oneIIndex,
							temperature)
					+ calculateDissociationConstant(thisIndex, oneHeIndex,
							temperature)) * getConcentration();
}

std::vector<int> HeInterstitialCluster::getConnectivity() {
	
	// Local Declarations
	std::shared_ptr<std::map<std::string, std::string>> properties =
		network->properties;
	
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	int numVClusters = std::stoi(properties->at("numVClusters"));
	int numIClusters = std::stoi(properties->at("numIClusters"));
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;
	
	int numHeVClusters = std::stoi(properties->at("numHeVClusters"));
	
	// Initialize the return array with zeroes
	std::vector<int> connectivityArray(network->reactants->size(), 0);
	
	// This cluster's index in the reactants array
	int clusterIndex = numSingleSpeciesClusters + numHeVClusters +
		(numHe - 1) + (numI - 1) * numHeVClusters;
	
	
	// This cluster is involved in the following interactions:
	
	// Identity reaction
	connectivityArray.at(clusterIndex) = 1;
	
	// xHe * yI + zHe --> [x+z]He * yI
	for (int z = 1; numHe + z <= numHeClusters; z++) {
		connectivityArray.at(z - 1) = 1;
	}
	
	// xHe * yI + V --> xHe * [y - 1]I
	if (numI - 1 <= numVClusters) {
		connectivityArray.at(numHeClusters) = 1;
	}
	
	// xHe * yI + zI  --> xHe * [y + z]I
	for (int z = 1; numI + z <= numVClusters; z++) {
		connectivityArray.at(numHeClusters + numVClusters + z - 1) = 1;
	}
	
	
	return connectivityArray;
}
