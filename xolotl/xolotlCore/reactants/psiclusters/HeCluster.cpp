// Includes
#include "HeCluster.h"
#include <iostream>

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
	int maxHeSize = strtol(props["maxHeClusterSize"].c_str(),NULL,10);
	int maxMixedSize = strtol(props["maxMixedClusterSize"].c_str(),NULL,10);
	int numHe = strtol(props["numHeClusters"].c_str(),NULL,10);
	int numV = strtol(props["numVClusters"].c_str(),NULL,10);
	int vSize = 1;
	std::vector<int> connectivityArray(network->reactants->size(),0);

	// ----- A*He + B*He --> (A+B)*He -----
	// This cluster should interact with all other clusters of the same type up
	// the max size minus one to produce larger clusters.
	for (int i = 0; i < numHe - 1; i++)
		connectivityArray.at(i) = 1;

	// -----  A*He + B*V → (A*He)(B*V) -----
	// Helium clusters can interact with any vacancy cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	int clusterIndex = numHe;
	// Get the size of the first vacancy cluster
	vSize = (std::dynamic_pointer_cast<PSICluster>(network->reactants->at(clusterIndex)))->getSize();
	while(size + vSize < maxMixedSize && clusterIndex < numHe + numV) {
		std::cout << clusterIndex << std::endl;
		// Set the connectivity to 1
		connectivityArray.at(clusterIndex) = 1;
		// Increment the counter
		clusterIndex++;
		// Get the size of the next vacancy cluster
		vSize = (std::dynamic_pointer_cast<PSICluster>(network->reactants->at(clusterIndex)))->getSize();
	}

	return connectivityArray;
}

