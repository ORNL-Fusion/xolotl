// Includes
#include "HeInterstitialCluster.h"

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(std::map<std::string, int> speciesMap) :
		PSICluster(1) {
	numHe = speciesMap["He"];
	numI = speciesMap["I"];
	size = numHe + numI;
	// Set the reactant name appropriately
	name = "He-Interstitial Cluster";
}
HeInterstitialCluster::~HeInterstitialCluster() {
	//TODO Auto-generated method stub
}
double HeInterstitialCluster::getGenByEm() {
	//TODO Auto-generated method stub
	return 0;
}
double HeInterstitialCluster::getAnnByEm() {
	//TODO Auto-generated method stub
	return 0;
}
int HeInterstitialCluster::getSpeciesSize(const std::string speciesName) {
	//TODO Auto-generated method stub
	return 0;
}

std::vector<int> HeInterstitialCluster::getConnectivity() {
	
	// Local Declarations
	std::shared_ptr<std::map<std::string, std::string>> properties =
		network->properties;
	
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	int numVClusters = std::stoi(properties->at("numVClusters"));
	int numIClusters = std::stoi(properties->at("numIClusters"));
	int numSingleSpeciesClusters = numHeClusters + numVClusters + numIClusters;
	
	int numHeVClusters = std::stoi(properties->at("numHeVClusters"));
	
	// Initialize the return array with zeroes
	std::vector<int> connectivityArray(network->reactants->size(), 0);
	
	// This cluster's index in the reactants array
	int clusterIndex = numSingleSpeciesClusters + numHeVClusters +
		(numHe - 1) + (numI - 1) * numHeVClusters;
	
	
	// This cluster is involved in the following interactions:
	
	// Identity reaction
	connectivityArray.at(clusterIndex) = 1;
	
	// xHe * yI + zHe --> [x+z]He * yI
	for (int z = 1; numHe + z <= numHeClusters; z++) {
		connectivityArray.at(z - 1) = 1;
	}
	
	// xHe * yI + V --> xHe * [y - 1]I
	if (numI - 1 <= numVClusters) {
		connectivityArray.at(numHeClusters) = 1;
	}
	
	// xHe * yI + zI  --> xHe * [y + z]I
	for (int z = 1; numI + z <= numVClusters; z++) {
		connectivityArray.at(numHeClusters + numVClusters + z - 1) = 1;
	}
	
	
	return connectivityArray;
}
