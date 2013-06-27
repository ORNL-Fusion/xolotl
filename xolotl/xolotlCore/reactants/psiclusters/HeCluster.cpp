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
	std::vector<int> connectivityArray;
	std::map<std::string, std::string> props = *(network->properties);
	int maxHeSize = strtol(props["maxHeClusterSize"].c_str(),NULL,10);

	// This cluster should interact with all other clusters of the same type up
	// the max size minus one.
	//for (int i = 0; i < )

	return connectivityArray;
}

