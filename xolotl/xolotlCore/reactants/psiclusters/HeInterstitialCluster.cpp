// Includes
#include "HeInterstitialCluster.h"

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(int numHelium, int numInterstitial) :
		PSICluster(1) {
	numHe = numHelium;
	numI = numInterstitial;
	size = numHe + numI;
	name = "He-Interstitial Cluster";
}

HeInterstitialCluster::HeInterstitialCluster(
		std::map<std::string, int> speciesMap) :
		PSICluster(1) {
	numHe = speciesMap["He"];
	numI = speciesMap["I"];
	size = numHe + numI;
	// Set the reactant name appropriately
	name = "He-Interstitial Cluster";
}

HeInterstitialCluster::~HeInterstitialCluster() {
}

double HeInterstitialCluster::getGenByEm() {
	return 0;
}

double HeInterstitialCluster::getAnnByEm() {
	return 0;
}

int HeInterstitialCluster::getSpeciesSize(const std::string speciesName) {
	return 0;
}

std::vector<int> HeInterstitialCluster::getReactionConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> &properties = *network->properties;

	int numHeClusters = std::stoi(properties["numHeClusters"]);
	int numVClusters = std::stoi(properties["numVClusters"]);
	int numIClusters = std::stoi(properties["numIClusters"]);
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;

	int maxMixedClusterSize = std::stoi(properties["maxMixedClusterSize"]);

	// Initialize the return array with zeroes
	std::vector<int> connectivityArray(network->reactants->size(), 0);

	// This cluster is involved in the following interactions:

	// Growth through helium absorption
	// xHe * yI + zHe --> (x+z)He * yI
	// under the condition that x + y + z <= maxSize
	for (int numHeOther = 1; numHe + numI + numHeOther <= numHeClusters;
			numHeOther++) {

		std::map<std::string, int> speciesMap;
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}

	// Interstitial absorption (single)
	// xHe * yI + I --> xHe * (y + 1)I
	// if x + y + 1 <= maxSize
	if (numHe + numI + 1 <= maxMixedClusterSize) {

		std::map<std::string, int> speciesMap;
		speciesMap["I"] = 1;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}

	// Reduction through vacancy absorption
	// xHe * yI + zV --> xHe * (y - z)I
	for (int numVOther = 1; numI - numVOther >= 1; numVOther++) {

		std::map<std::string, int> speciesMap;
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}

	return connectivityArray;
}

std::vector<int> HeInterstitialCluster::getDissociationConnectivity() {
	// Local Declarations
	int nReactants = network->reactants->size();
	std::vector<int> dissConnections(nReactants, 0);
	std::map<std::string, int> clusterMap;

	return dissConnections;
}

std::map<std::string, int> HeInterstitialCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = numHe;
	clusterMap["V"] = 0;
	clusterMap["I"] = numI;

	// Return it
	return clusterMap;
}

double HeInterstitialCluster::getReactionRadius() {
	return 0.0;
}
