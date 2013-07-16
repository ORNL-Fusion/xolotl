// Includes
#include "VCluster.h"
#include <Constants.h>

using namespace xolotlCore;

VCluster::VCluster(int nV) :
		PSICluster(nV) {
	// Set the reactant name appropriately
	name = "Vacancy";
}

VCluster::~VCluster() {
}

std::vector<int> VCluster::getConnectivity() {
	
	std::map<std::string, std::string> &props = *(network->properties);
	
	int numV = size;
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int numVClusters = std::stoi(props["numVClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	
	// Initialize the connectivity row with zeroes
	int reactantsLength = network->reactants->size();
	std::vector<int> connectivityArray(reactantsLength, 0);
	
	
	// Vacancies interact with everything except for vacancies bigger than they
	// would combine with to form vacancies larger than the size limit.
	
	// -----  A*He + B*V → (A*He)(B*V) -----
	// Vacancy clusters can interact with any helium cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	for (int numHeOther = 1; numV + numHeOther <= maxMixedClusterSize; numHeOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	//----- A*V + B*V --> (A+B)*V -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	for (int numVOther = 1; numV + numVOther <= maxVClusterSize; numVOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	//----- A*I + B*V
	// → (A-B)*I, if A > B
	// → (B-I)*V, if A < B
	// → 0, if A = B -----
	// Vacancies always annihilate interstitials.
	for (int numIOther = 1; numIOther <= numIClusters; numIOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["I"] = numIOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	// ----- (A*He)(B*V) + C*V → (A*He)[(B+C)*V] -----
	// Vacancies can interact with a mixed-species cluster so long as the sum of
	// the number of vacancy atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	if (numV == 1) {
		for (int numVOther = 1; numVOther <= maxMixedClusterSize; numVOther++) {
			for (int numHeOther = 1; (numHeOther + numVOther + numV) <=
				maxMixedClusterSize; numHeOther++) {
				
				std::map<std::string, int> speciesMap;
				speciesMap["He"] = numHeOther;
				speciesMap["V"] = numVOther;
				int indexOther = network->toClusterIndex(speciesMap);
				connectivityArray[indexOther] = 1;
			}
		}
	}
	
	return connectivityArray;
}

double VCluster::getDissociationFlux(const double temperature) {

	// Local Declarations
	double diss = 0.0;
	int numV = 0, deltaIndex = -1;
	std::map<std::string, int> oneV;

	// Set the cluster map data for 1 of each species
	oneV["He"] = 0; oneV["V"] = 1; oneV["I"] = 0;

	// Get their indices in the array
	int oneVIndex = network->toClusterIndex(oneV);

	// Loop over all Reactants
	for (int j = 0; j < network->reactants->size(); j++) {

		// Get the number of vacancies in C_j
		numV = network->toClusterMap(j)["V"];

		// If the C_j contains a Vacancy, then we calculate
		if (numV > 0) {
			// Search for the index of the cluster that contains exactly
			// one less Vacancy than C_j, once found break from the loop
			for (int k = 0; k < network->reactants->size(); k++) {
				if ((numV - network->toClusterMap(k)["V"]) == 1) {
					deltaIndex = k;
					break;
				}
			}

			// There may not have been an index that had one less
			// Vacancy, if so, we won't add to the dissociation flux
			if (deltaIndex != -1) {
				// Calculate the dissociation, with K^- evaluated
				// at deltaIndex and this Vacancy Cluster's index.
				diss = diss	+ calculateDissociationConstant(deltaIndex, oneVIndex,
								temperature) * network->reactants->at(j)->getConcentration();
			}
		}
	}

	// Return the dissociation
	return diss;
}

bool VCluster::isProductReactant(int reactantI, int reactantJ) {

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
	// total of 0 Helium, and a total of
	// size Vacancies
	return ((rI_I + rJ_I) == 0) && ((rI_He + rJ_He) == 0)
			&& ((rI_V + rJ_V) == size);
}

std::map<std::string, int> VCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = size;
	clusterMap["I"] = 0;

	// Return it
	return clusterMap;
}

double VCluster::getReactionRadius() {
	// FIXME Not right...
	return (sqrt(3)/4) * xolotlCore::latticeConstant;
}
