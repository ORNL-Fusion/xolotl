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

void HeInterstitialCluster::createReactionConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> properties = *network->properties;

	int numHeClusters = std::stoi(properties["numHeClusters"]);
	int numVClusters = std::stoi(properties["numVClusters"]);
	int numIClusters = std::stoi(properties["numIClusters"]);
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;
	int maxMixedClusterSize = std::stoi(properties["maxMixedClusterSize"]);
	std::map<std::string, int> speciesMap;

	// This cluster is involved in the following interactions:

	// Growth through helium absorption
	// xHe * yI + zHe --> (x+z)He * yI
	// under the condition that x + y + z <= maxSize
	for (int numHeOther = 1; numHe + numI + numHeOther <= numHeClusters;
			numHeOther++) {
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}

	// Interstitial absorption (single)
	// xHe * yI + I --> xHe * (y + 1)I
	// if x + y + 1 <= maxSize
	if (numHe + numI + 1 <= maxMixedClusterSize) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["I"] = 1;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}

	// Reduction through vacancy absorption
	// xHe * yI + zV --> xHe * (y - z)I
	for (int numVOther = 1; numI - numVOther >= 1; numVOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}
}

void HeInterstitialCluster::createDissociationConnectivity() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Vacancy Dissociation
	clusterMap["He"] = numHe - 1;
	clusterMap["V"] = 0;
	clusterMap["I"] = numI;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["I"] = 0;
	clusterMap["He"] = 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;

	// Trap Mutation
	clusterMap["I"] = numI + 1;
	clusterMap["He"] = numHe;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["I"] = 1;
	clusterMap["V"] = 0;
	clusterMap["He"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;

	// Vacancy Dissociation
	clusterMap["He"] = numHe;
	clusterMap["V"] = 0;
	clusterMap["I"] = numI - 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["He"] = 0;
	clusterMap["V"] = 1;
	clusterMap["I"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
}

double HeInterstitialCluster::getDissociationFlux(double temperature) {
	// Local Declarations
	std::map<std::string, int> oneHe, oneV, oneI, dissMap;
	double f4 = 0.0, f3 = 0.0;

	// Set the cluster map data for 1 of each species
	oneHe["He"] = 1;
	oneHe["V"] = 0;
	oneHe["I"] = 0;
	oneV["He"] = 0;
	oneV["V"] = 1;
	oneV["I"] = 0;
	oneI["He"] = 0;
	oneI["V"] = 0;
	oneI["I"] = 1;

	// Get this PSICluster or subclasses' cluster map
	std::map<std::string, int> thisMap = getClusterMap();

	// Get the various indices
	int thisIndex = network->toClusterIndex(thisMap);
	int oneIIndex = network->toClusterIndex(oneI);
	int oneVIndex = network->toClusterIndex(oneV);
	int oneHeIndex = network->toClusterIndex(oneHe);

	// Calculate the much easier f4 term...
	f4 = calculateDissociationConstant(thisIndex, oneIIndex, temperature)
			+ calculateDissociationConstant(thisIndex, oneVIndex, temperature)
			+ calculateDissociationConstant(thisIndex, oneHeIndex, temperature);

	// Loop over all the elements of the dissociation
	// connectivity to find where this mixed species dissociates...
	for (int i = 0; i < dissociationConnectivity.size(); i++) {
		if (dissociationConnectivity[i] == 1) {
			// Get the cluster map of this connection
			dissMap = network->toClusterMap(i);

			// We need to find if this is a Helium dissociation,
			// Vacancy dissociation, or a trap mutation.
			if (numHe - dissMap["He"] == 1 && numI == dissMap["I"]
					&& dissMap["V"] == 0) {
				f3 = f3
						+ calculateDissociationConstant(i, oneHeIndex,
								temperature)
								* network->reactants->at(i)->getConcentration();
			} else if (numHe == dissMap["He"] && numI - dissMap["V"] == 1
					&& dissMap["V"] == 0) {
				f3 = f3
						+ calculateDissociationConstant(i, oneVIndex,
								temperature)
								* network->reactants->at(i)->getConcentration();
			} else if (numHe == dissMap["He"] && dissMap["I"] - numI == 1
					&& dissMap["V"] == 0) {
				f3 = f3
						+ calculateDissociationConstant(i, oneIIndex,
								temperature)
								* network->reactants->at(i)->getConcentration();
			}

		}
	}

	return f3 - f4 * getConcentration();
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
