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

double InterstitialCluster::getDissociationFlux(const double temperature) {
	// Local Declarations
	double diss = 0.0;
	int numI = 0, deltaIndex = -1;

	for (int j = 0; j < network->reactants->size(); j++) {
		numI = network->toClusterMap(j)["I"];
		// If the Jth reactant contains Intersitials, then we calculate
		if (numI > 0) {
			// Search for the index of the cluster that contains exactly
			// one less Interstitial than reactant->at(j)
			for (int k = 0; k < network->reactants->size(); k++) {
				if ((network->toClusterMap(k)["I"] - numI) == 1) {
					deltaIndex = k;
				}
			}

			// There may not have been an index that had one less
			// Interstitial, if so, we won't add to the dissociation flux
			if (deltaIndex != -1) {
				// Calculate the dissociation, with K^- evaluated
				// at deltaIndex and this Intersitial Cluster's index.
				diss = diss + calculateDissociationConstant(j,
								network->toClusterIndex(getClusterMap()),
								temperature) * network->reactants->at(j)->getConcentration();
			}
		}
	}

	// Return the dissociation
	return diss;
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
