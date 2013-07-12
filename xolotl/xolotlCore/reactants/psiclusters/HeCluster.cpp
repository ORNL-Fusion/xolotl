// Includes
#include "HeCluster.h"

using namespace xolotlCore;

HeCluster::HeCluster(int nHe) :
		PSICluster(nHe) {
	// Set the reactant name appropriately
	name = "Helium";
}

HeCluster::~HeCluster() {
}

std::vector<int> HeCluster::getConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> props = *(network->properties);
	int maxMixedSize = strtol(props["maxMixedClusterSize"].c_str(), NULL, 10);
	int numHe = strtol(props["numHeClusters"].c_str(), NULL, 10);
	int numV = strtol(props["numVClusters"].c_str(), NULL, 10);
	int numI = strtol(props["numIClusters"].c_str(), NULL, 10);
	int numMixed = strtol(props["numMixedClusters"].c_str(), NULL, 10);
	int clusterIndex = 0;
	int vSize = 1, mixedSize = 1;
	std::vector<int> connectivityArray(network->reactants->size(), 0);

	// ----- A*He + B*He --> (A+B)*He -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	for (int i = 1; i + numHe <= size; i++)
		connectivityArray.at(i - 1) = 1;

	// -----  A*He + B*V → (A*He)(B*V) -----
	// Helium clusters can interact with any vacancy cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	clusterIndex = numHe;
	// Get the size of the first vacancy cluster
	vSize = (std::dynamic_pointer_cast < PSICluster
			> (network->reactants->at(clusterIndex)))->getSize();
	// Loop over all of the valid vacancy clusters
	while (size + vSize <= maxMixedSize && clusterIndex < numHe + numV) {
		// Set the connectivity to 1
		connectivityArray.at(clusterIndex) = 1;
		// Increment the counter
		clusterIndex++;
		// Get the size of the next vacancy cluster
		vSize = (std::dynamic_pointer_cast < PSICluster
				> (network->reactants->at(clusterIndex)))->getSize();
	}

	// ----- (A*He)(B*V) + C*He → [(A+C)He](B*V) -----
	// Helium can interact with a mixed-species cluster so long as the sum of
	// the number of helium atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	clusterIndex = numHe + numV + numI;
	// Get the size of the first vacancy cluster
	mixedSize = (std::dynamic_pointer_cast < PSICluster
			> (network->reactants->at(clusterIndex)))->getSize();
	// Loop over all of the mixed clusters -- FIX BOUNDS!
	while (size + mixedSize <= maxMixedSize
			&& clusterIndex < numHe + numV + numI + numMixed) {
		// Set the connectivity to 1
		connectivityArray.at(clusterIndex) = 1;
		// Increment the counter
		clusterIndex++;
		// Get the size of the next mixed cluster
		mixedSize = (std::dynamic_pointer_cast < PSICluster
				> (network->reactants->at(clusterIndex)))->getSize();
	}

	return connectivityArray;
}

double HeCluster::getDissociationFlux(const double temperature) {

	// Local Declarations
	double diss = 0.0;
	int numHelium = 0, deltaIndex = -1;

	// Loop over all reactants
	for (int j = 0; j < network->reactants->size(); j++) {

		// Get the number of helium species in the jth reactant
		numHelium = network->toClusterMap(j)["He"];

		// If the Jth reactant contains Helium, then we calculate
		if (numHelium > 0) {
			// Search for the index of the cluster that contains exactly
			// one less helium than reactant->at(j)
			for (int k = 0; k < network->reactants->size(); k++) {
				if ((network->toClusterMap(k)["He"] - numHelium) == 1) {
					// Once found, get the current index
					deltaIndex = k;
				}
			}

			// There may not have been an index that had one less
			// helium, if so, we won't add to the dissociation flux
			if (deltaIndex != -1) {
				// Calculate the dissociation, with K^- evaluated
				// at deltaIndex and this Helium Cluster's index.
				diss = diss + calculateDissociationConstant(j,
								network->toClusterIndex(getClusterMap()),
								temperature) * network->reactants->at(j)->getConcentration();
			}
		}
	}

	// Return the dissociation
	return diss;
}

bool HeCluster::isProductReactant(int reactantI, int reactantJ) {

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

	// We should have no interstitials, a
	// total of size Helium, and a total of
	// 0 Vacancies
	return ((rI_I + rJ_I) == 0) && ((rI_He + rJ_He) == size)
			&& ((rI_V + rJ_V) == 0);
}

std::map<std::string, int> HeCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = size;
	clusterMap["V"] = 0;
	clusterMap["I"] = 0;

	// Return it
	return clusterMap;
}
