// Includes
#include "HeVCluster.h"

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV) :
	PSICluster(1), numHe(numHe), numV(numV) {
	
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
	if (speciesName == "He") {
		return numHe;
	}
	else if (speciesName == "V") {
		return numV;
	}
	else {
		return 0;
	}
}


std::vector<int> HeVCluster::getConnectivity() {
	
	// Extract some of the properties from the network
	
	std::shared_ptr<std::map<std::string, std::string>> properties =
		network->properties;
	
	int maxHeClusterSize = std::stoi(properties->at("maxHeClusterSize"));
	int maxVClusterSize = std::stoi(properties->at("maxVClusterSize"));
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	int numVClusters = std::stoi(properties->at("numVClusters"));
	int numIClusters = std::stoi(properties->at("numIClusters"));
	
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;
	
	// Initialize the return array with zeroes
	std::vector<int> connectivityArray(network->reactants->size(), 0);
	
	// This cluster's index in the reactants array
	int clusterIndex = numSingleSpeciesClusters + (numHe - 1) + (numV - 1) * numHeClusters;
	
	
	// This cluster is involved in the following interactions:
	
	// HeV[x, y] + He[z] --> HeV[x + z, y]
	for (int z = 1; numHe + z <= maxHeClusterSize; z++) {
		// Select He[z]
		connectivityArray.at(z - 1) = 1;
	}
	
	// HeV[x, y] + V[1]  --> HeV[x, y + 1]
	if (numV + 1 <= maxVClusterSize) {
		// Select V[1]
		connectivityArray.at(numHeClusters) = 1;
	}
	
	// HeV[x, y] + I[z]  --> HeV[x, y - z]
	for (int z = 1; numV - z >= 1; z++) {
		// Select I[z]
		connectivityArray.at(numHeClusters + numVClusters + z - 1) = 1;
	}
	
	
	return connectivityArray;
}
