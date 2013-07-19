#include "InterstitialCluster.h"
#include <Constants.h>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI) :
		PSICluster(nI) {
	// Set the reactant name appropriately
	name = "Interstitial";
}
InterstitialCluster::~InterstitialCluster() {
}

void InterstitialCluster::createReactionConnectivity() {

	std::map<std::string, std::string> &props = *network->properties;
	
	int numI = size;
	int numVClusters = std::stoi(props["numVClusters"]);
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	std::map<std::string, int> speciesMap;
	
	// Initialize the connectivity row with zeroes
	int reactantsLength = network->reactants->size();
	for (int i = 0; i < reactantsLength; i++) {
		reactionConnectivity.push_back(0);
	}
	//reactionConnectivity.resize(reactantsLength, 0);
	
	// Interstitials can interact with other interstitials, vacancies,
	// helium, and mixed-species clusters. They cannot cluster with other
	// interstitials that are so large that the combination of the two would
	// produce an interstitial above the maximum size.
	
	// xHe + yI --> xHe*yI
	for (int numHeOther = 1; numHeOther + numI <= maxMixedClusterSize; numHeOther++) {
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}
	
	//----- A*I + B*V
	// → (A-B)*I, if A > B
	// → (B-I)*V, if A < B
	// → 0, if A = B
	// Annihilation
	for (int numVOther = 1; numVOther <= numVClusters; numVOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}
	
	// A*I + B*I → (A+B)*I -----
	for (int numIOther = 1; numI + numIOther <= maxIClusterSize; numIOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["I"] = numIOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}
	
	// ----- (A*He)(B*V) + (C*I) --> (A*He)[(B-C)V] -----
	// Interstitials interact with all mixed-species clusters by
	// annihilating vacancies.
	for (int numVOther = 1; numVOther <= maxMixedClusterSize; numVOther++) {
		for (int numHeOther = 1; numVOther + numHeOther <= maxMixedClusterSize;
			numHeOther++) {
			// Clear the map since we are reusing it
			speciesMap.clear();
			bool connected = numVOther - numI >= 1;
			speciesMap["He"] = numHeOther;
			speciesMap["V"] = numVOther;
			int indexOther = network->toClusterIndex(speciesMap);
			reactionConnectivity[indexOther] = (int) connected;
		}
	}
	
	// Interstitial absorption
	// xHe*yI + I --> xHe*(y + 1)I
	// Under the condition that (x + y + 1) <= maxSize
	if (numI == 1) {
		for (int numIOther = 1; numIOther <= maxMixedClusterSize; numIOther++) {
			for (int numHeOther = 1; numIOther + numHeOther + 1 <=
				maxMixedClusterSize; numHeOther++) {
				// Clear the map since we are reusing it
				speciesMap.clear();
				speciesMap["He"] = numHeOther;
				speciesMap["I"] = numIOther;
				int indexOther = network->toClusterIndex(speciesMap);
				reactionConnectivity[indexOther] = 1;
			}
		}
	}
}


void InterstitialCluster::createDissociationConnectivity() {
	
	// Resize the connectivity row with zeroes
	int reactantsLength = network->reactants->size();
	dissociationConnectivity.resize(reactantsLength, 0);
}


bool InterstitialCluster::isProductReactant(int reactantI, int reactantJ) {

	// Local Declarations, integers for species number for I, J reactants
	int rI_I = 0, rJ_I = 0, rI_He = 0, rJ_He = 0, rI_V = 0, rJ_V = 0;

	// Get the ClusterMap corresponding to
	// the given reactants
	std::map<std::string, int> reactantIMap = network->toClusterMap(reactantI);
	std::map<std::string, int> reactantJMap = network->toClusterMap(reactantJ);

	// Grab the numbers for each species
	// from each Reactant
	rI_I = reactantIMap["I"];
	rJ_I = reactantJMap["I"];
	rI_He = reactantIMap["He"];
	rJ_He = reactantJMap["He"];
	rI_V = reactantIMap["V"];
	rJ_V = reactantJMap["V"];

	// We should have this->size interstitials, a
	// total of 0 Helium, and a total of
	// 0 Vacancies
	return ((rI_I + rJ_I) == size) && ((rI_He + rJ_He) == 0)
			&& ((rI_V + rJ_V) == 0);
}

std::map<std::string, int> InterstitialCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = 0;
	clusterMap["I"] = size;

	// Return it
	return clusterMap;
}

double InterstitialCluster::getReactionRadius() {

	double EightPi = 8.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3.0);
	double termOne = 1.15*(sqrt(3.0)/4.0) * xolotlCore::latticeConstant;
	double termTwo = pow( (3.0/EightPi) * aCubed * size, (1.0/3.0));
	double termThree = pow( (3.0/EightPi) * aCubed, (1.0/3.0));

	return termOne + termTwo - termThree;
}
