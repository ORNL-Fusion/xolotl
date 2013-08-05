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

void VCluster::createReactionConnectivity() {

	std::map<std::string, std::string> props = *(network->properties);

	int numV = size;
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int numVClusters = std::stoi(props["numVClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	std::map<std::string, int> speciesMap;

	// Vacancies interact with everything except for vacancies bigger than they
	// would combine with to form vacancies larger than the size limit.

	// -----  A*He + B*V → (A*He)(B*V) -----
	// Vacancy clusters can interact with any helium cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	for (int numHeOther = 1; numV + numHeOther <= maxMixedClusterSize;
			numHeOther++) {
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}

	//----- A*V + B*V --> (A+B)*V -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	for (int numVOther = 1; numV + numVOther <= maxVClusterSize; numVOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}

	//----- A*I + B*V
	// → (A-B)*I, if A > B
	// → (B-I)*V, if A < B
	// → 0, if A = B -----
	// Vacancies always annihilate interstitials.
	for (int numIOther = 1; numIOther <= numIClusters; numIOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["I"] = numIOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
	}

	// ----- (A*He)(B*V) + C*V → (A*He)[(B+C)*V] -----
	// Vacancies can interact with a mixed-species cluster so long as the sum of
	// the number of vacancy atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	if (numV == 1) {
		for (int numVOther = 1; numVOther <= maxMixedClusterSize; numVOther++) {
			for (int numHeOther = 1;
					(numHeOther + numVOther + numV) <= maxMixedClusterSize;
					numHeOther++) {
				// Clear the map since we are reusing it
				speciesMap.clear();
				speciesMap["He"] = numHeOther;
				speciesMap["V"] = numVOther;
				int indexOther = network->toClusterIndex(speciesMap);
				reactionConnectivity[indexOther] = 1;
			}
		}
	}

	// Vacancy absorption by HeI:
	// xHe*yI + zV --> xHe*(y - z)V
	// under the condition that y - z >= 1
	if (numHeIClusters > 0) {
		for (int numIOther = 1; numIOther <= maxMixedClusterSize; numIOther++) {
			for (int numHeOther = 1;
					numIOther + numHeOther <= maxMixedClusterSize;
					numHeOther++) {
				// Clear the map since we are reusing it
				speciesMap.clear();
				bool connects = numIOther - numV >= 1;
				speciesMap["He"] = numHeOther;
				speciesMap["I"] = numIOther;
				int indexOther = network->toClusterIndex(speciesMap);
				reactionConnectivity[indexOther] = (int) connects;
			}
		}
	}

}

void VCluster::createDissociationConnectivity() {
	// Local Declarations
	int nReactants = network->reactants->size();
	std::map<std::string, int> clusterMap;

	// Vacancy Dissociation
	clusterMap["He"] = 0;
	clusterMap["V"] = size - 1;
	clusterMap["I"] = 0;
	if (size != 1) {
		dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
		clusterMap["V"] = 1;
		dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	}

	// Trap Mutation
	clusterMap["V"] = size + 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["I"] = 1;
	clusterMap["V"] = 0;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
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
	return (sqrt(3) / 4) * xolotlCore::latticeConstant;
}
