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
	
	// Note the reference to the properties map
	std::map<std::string, std::string> &props = *(network->properties);
	
	int numHe = size;
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	
	// Initialize the connectivity row with zeroes
	int reactantsLength = network->reactants->size();
	std::vector<int> connectivityArray(reactantsLength, 0);
	
	// ----- A*He + B*He --> (A+B)*He -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	for (int numHeOther = 1; numHe + numHeOther <= maxHeClusterSize; numHeOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	// -----  A*He + B*V --> (A*He)(B*V) -----
	// Helium clusters can interact with any vacancy cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	
	for (int numVOther = 1; numHe + numVOther <= maxMixedClusterSize; numVOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	// ----- (A*He)(B*V) + C*He --> [(A+C)He](B*V) -----
	// Helium can interact with a mixed-species cluster so long as the sum of
	// the number of helium atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	
	// Get the index of the first HeV cluster in the reactants list
	
	std::map<std::string, int> speciesMap;
	speciesMap["He"] = 1;
	speciesMap["V"] = 1;
	int firstHeVIndex = network->toClusterIndex(speciesMap);
	
	for (int indexOther = firstHeVIndex; indexOther < reactantsLength; indexOther++) {
		
		std::map<std::string, int> speciesMap = network->toClusterMap(indexOther);
		int numHeOther = speciesMap["He"];
		int numVOther = speciesMap["V"];
		
		// Check if the sum of this He and HeV is no larger than the
		// maximum allowed 
		if (numHe + numHeOther + numVOther <= maxMixedClusterSize) {
			connectivityArray[indexOther] = 1;
		}
	}
	
	return connectivityArray;
}

double HeCluster::getDissociationFlux(const double temperature) {

	// Local Declarations
	double diss = 0.0;
	int numHelium = 0, deltaIndex = -1;
	std::map<std::string, int> oneHe;

	// Set the cluster map data for 1 of each species
	oneHe["He"] = 1; oneHe["V"] = 0; oneHe["I"] = 0;

	// Get their indices in the array
	int oneHeIndex = network->toClusterIndex(oneHe);

	// Loop over all reactants
	for (int j = 0; j < network->reactants->size(); j++) {

		// Get the number of helium species in C_j
		numHelium = network->toClusterMap(j)["He"];

		// If the C_j contains Helium, then we calculate
		if (numHelium > 0) {
			// Search for the index of the cluster that contains exactly
			// one less helium than C_j, then break from the loop
			for (int k = 0; k < network->reactants->size(); k++) {
				if ((numHelium - network->toClusterMap(k)["He"]) == 1) {
					// Once found, get the current index
					deltaIndex = k;
					break;
				}
			}

			// There may not have been an index that had one less
			// helium, if so, we won't add to the dissociation flux
			if (deltaIndex != -1) {
				// Calculate the dissociation, with K^- evaluated
				// at deltaIndex and this Helium Cluster's index.
				diss = diss + calculateDissociationConstant(deltaIndex, oneHeIndex,
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
