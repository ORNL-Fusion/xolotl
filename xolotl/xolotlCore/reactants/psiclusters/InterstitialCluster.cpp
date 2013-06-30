#include "InterstitialCluster.h"

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI) :
		PSICluster(nI) {
	// Set the reactant name appropriately
	name = "Interstitial";
}
InterstitialCluster::~InterstitialCluster() {
}

std::vector<int> InterstitialCluster::getConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> props = *(network->properties);
	int numHe = strtol(props["numHeClusters"].c_str(), NULL, 10);
	int numV = strtol(props["numVClusters"].c_str(), NULL, 10);
	int numI = strtol(props["numIClusters"].c_str(), NULL, 10);
	std::vector<int> connectivityArray(network->reactants->size(), 0);

	// Interstitials can interact with other interstitials, vacancies and
	// mixed-species clusters, but not helium. They cannot cluster with other
	// interstitials that are so large that the combination of the two would
	// produce an interstitial above the maximum size.

	//----- A*I + B*V
	// → (A-B)*I, if A > B
	// → (B-I)*V, if A < B
	// → 0, if A = B
	//
	// AND A*I + B*I → (A+B)*I -----
	//
	// Vacancies and interstitials always annihilate each other. Since
	// vacancies and interstitials are next to each other in the array,
	// the connectivity can be combined to go up to the largest
	// interstitial that can combine with this one without breaking the
	// size limit.
	for (int i = numHe; i < numHe + numV + numI - size; i++) {
		connectivityArray.at(i) = 1;
	}

	// ----- (A*He)(B*V) + (C*I) → (A*He)[(B-C)V] -----
	// Interstitials interact with all mixed-species clusters by
	// annihilating vacancies.
	for (int i = numHe + numV + numI; i < connectivityArray.size(); i++) {
		connectivityArray.at(i) = 1;
	}

	return connectivityArray;
}
