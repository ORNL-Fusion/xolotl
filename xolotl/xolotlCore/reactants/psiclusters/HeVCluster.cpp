// Includes
#include "HeVCluster.h"

using namespace xolotlCore;

HeVCluster::HeVCluster(std::map<std::string, int> speciesMap) :
		PSICluster(1) {
	numHe = speciesMap["He"];
	numV = speciesMap["V"];
	size = numHe + numV;
	// Set the reactant name appropriately
	name = "HeV Cluster";
}
HeVCluster::~HeVCluster() {
	//TODO Auto-generated method stub
}
double HeVCluster::getGenByEm() {
	//TODO Auto-generated method stub
	return 0;
}
double HeVCluster::getAnnByEm() {
	//TODO Auto-generated method stub
	return 0;
}
int HeVCluster::getSpeciesSize(const std::string speciesName) {
	//TODO Auto-generated method stub
	return 0;
}

std::vector<int> HeVCluster::getConnectivity() {

	// Local Declarations
	std::map<std::string, std::string> props = *(network->properties);
	int totalNumHeClusters = strtol(props["numHeClusters"].c_str(), NULL, 10);
	int totalNumVClusters = strtol(props["numVClusters"].c_str(), NULL, 10);
	int totalNumIClusters = strtol(props["numIClusters"].c_str(), NULL, 10);
	int maxMixedSize = strtol(props["maxMixedClusterSize"].c_str(), NULL, 10);
	int numSingleSpeciesClusters = totalNumHeClusters + totalNumVClusters;
	std::vector<int> connectivityArray(network->reactants->size(), 0);
	
	// This cluster's index in reactants. It accounts for the interlaced array
	// and the zero-indexing.
	int clusterIndex = numSingleSpeciesClusters + (numHe + numV - 1) - 1;
	int oneLessHeIndex = 0, oneLessVIndex = 0, oneMoreVIndex = 0;

	// This cluster should only interact with a few specific clusters:
	// >single He - helium dissociation
	// >single V - vacancy dissociation
	// >single I - interstitial absorption
	// >[(A-1)*He](B*V) - helium dissociation
	// >(A*He)*[(B-1)*V] - vacancy dissociation
	// >(A*He)*[(B+1)*V] - interstitial absorption
	// via the following reactions:
	//(A*He)(B*V) →  [(A-1)*He](B*V) + He
	//(A*He)(B*V) →  (A*He)[(B-1)*V] + V
	//(A*He)(B*V) →  (A*He)[(B+1)*V] + I

	// Handle all of the single-species interactions first.
	connectivityArray.at(0) = 1;
	connectivityArray.at(totalNumHeClusters - 1) = 1;
	connectivityArray.at(totalNumHeClusters + totalNumVClusters - 1) = 1;

	// Figure out which mixed-clusters we need to update. There can never
	// be a case where V > He because of the ordering.
	if (numHe == numV) { // Equal He and V
		oneLessHeIndex = clusterIndex - 2;
		oneLessVIndex = clusterIndex - 1;
		oneMoreVIndex = clusterIndex + 2;
	} else { // More He
		oneLessHeIndex = clusterIndex - 1;
		oneLessVIndex = clusterIndex - 2;
		oneMoreVIndex = clusterIndex + 1;
	}
	// Set the connectivity.
	connectivityArray.at(oneLessHeIndex) = 1;
	connectivityArray.at(oneLessVIndex) = 1;
	// An additional vacancy can be added only if the cluster size limit is not
	// broken.
	if (size + 1 <= maxMixedSize)
		connectivityArray.at(oneMoreVIndex) = 1;

	return connectivityArray;
}
