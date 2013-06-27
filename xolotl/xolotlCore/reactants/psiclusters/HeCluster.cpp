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
	std::vector<int> connectivityArray(maxHeSize,0);

	std::cout << "Max He Cluster Size = " << props["maxHeClusterSize"] << std::endl;
	std::cout << maxHeSize << std::endl;
	std::cout << "Number of Reactants = " << network->reactants->size() << std::endl;

	// ----- A*He + B*He -> (A+B)*He -----
	// This cluster should interact with all other clusters of the same type up
	// the max size minus one to produce larger clusters.
	for (int i = 0; i < maxHeSize - 1; i++)
		connectivityArray.at(i) = 1;

	return connectivityArray;
}

