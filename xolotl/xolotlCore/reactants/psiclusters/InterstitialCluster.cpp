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
	int maxHeSize = strtol(props["maxHeClusterSize"].c_str(),NULL,10);
	int maxVSize = strtol(props["maxVClusterSize"].c_str(),NULL,10);
	std::vector<int> connectivityArray(maxVSize,0);

	// ----- A*He + B*He -> (A+B)*He -----
	// This cluster should interact with all other clusters of the same type up
	// the max size minus one to produce larger clusters.
	for (int i = 0; i < maxHeSize - 1; i++)
		connectivityArray.at(i) = 1;

	return connectivityArray;
}
