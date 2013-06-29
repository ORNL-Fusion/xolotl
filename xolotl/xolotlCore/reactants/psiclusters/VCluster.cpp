// Includes
#include "VCluster.h"
#include <iostream>

using namespace xolotlCore;

VCluster::VCluster(int nV) :
		PSICluster(nV) {
	// Set the reactant name appropriately
	name = "Vacancy";
}

VCluster::~VCluster() {
}

std::vector<int> VCluster::getConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> props = *(network->properties);
	int maxHeSize = strtol(props["maxHeClusterSize"].c_str(), NULL, 10);
	int maxMixedSize = strtol(props["maxMixedClusterSize"].c_str(), NULL, 10);
	int numHe = strtol(props["numHeClusters"].c_str(), NULL, 10);
	int numV = strtol(props["numVClusters"].c_str(), NULL, 10);
	int numI = strtol(props["numIClusters"].c_str(), NULL, 10);
	int numMixed = strtol(props["numMixedClusters"].c_str(), NULL, 10);
	int clusterIndex = 0;
	int heSize = 1, mixedSize = 1;
	std::vector<int> connectivityArray(network->reactants->size(), 0);

	// Vacancies interact with everything except for vacancies bigger than they
	// would combine with to form vacancies larger than the size limit.

	// -----  A*He + B*V → (A*He)(B*V) -----
	// Vacancy clusters can interact with any helium cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	clusterIndex = 0;
	// Get the size of the first helium cluster
	heSize = (std::dynamic_pointer_cast < PSICluster
			> (network->reactants->at(clusterIndex)))->getSize();
	// Loop over all of the valid helium clusters
	while (size + heSize <= maxMixedSize && clusterIndex < numHe) {
		// Set the connectivity to 1
		connectivityArray.at(clusterIndex) = 1;
		// Increment the counter
		clusterIndex++;
		// Get the size of the next vacancy cluster
		heSize = (std::dynamic_pointer_cast < PSICluster
				> (network->reactants->at(clusterIndex)))->getSize();
	}

	//----- A*V + B*V --> (A+B)*V -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	std::cout << "starting index = " << numHe << std::endl;
	for (int i = numHe; i < numHe + numV - size; i++) {
		std::cout << "vacancy index = " << i << std::endl;
		connectivityArray.at(i) = 1;
	}

	//----- A*I + B*V
	// → (A-B)*I, if A > B
	// → (B-I)*V, if A < B
	// → 0, if A = B -----
	// Vacancies always annihilate interstitials.
	for (int i = numHe + numV; i < numHe + numV + numI; i++) {
		connectivityArray.at(i) = 1;
	}

	// ----- (A*He)(B*V) + C*V → (A*He)[(B+C)*V] -----
	// Vacancies can interact with a mixed-species cluster so long as the sum of
	// the number of vacancy atoms and the size of the mixed-species cluster
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
