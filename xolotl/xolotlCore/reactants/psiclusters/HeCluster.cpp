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
	for (int i = 0; i < numHe - size; i++)
		connectivityArray.at(i) = 1;

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

