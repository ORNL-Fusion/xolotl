// Includes
#include "HeCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

HeCluster::HeCluster(int nHe) :
		PSICluster(nHe) {
	// Set the reactant name appropriately
	name = "Helium";
}

HeCluster::~HeCluster() {
}

std::vector<int> HeCluster::getReactionConnectivity() {
	
	// Note the reference to the properties map
	std::map<std::string, std::string> &props = *(network->properties);
	
	int numHe = size;
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	
	// Initialize the connectivity row with zeroes
	int reactantsLength = network->reactants->size();
	std::vector<int> connectivityArray(reactantsLength, 0);
	
	// ----- A*He + B*He --> (A+B)*He -----
	// This cluster should interact with all other clusters of the same type up
	// to the max size minus the size of this one to produce larger clusters.
	for (int numHeOther = 1; numHe + numHeOther <= maxHeClusterSize; numHeOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	// -----  A*He + B*V --> (A*He)(B*V) -----
	// Helium clusters can interact with any vacancy cluster so long as the sum
	// of the number of helium atoms and vacancies does not produce a cluster
	// with a size greater than the maximum mixed-species cluster size.
	
	for (int numVOther = 1; numHe + numVOther <= maxMixedClusterSize; numVOther++) {
		
		std::map<std::string, int> speciesMap;
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		connectivityArray[indexOther] = 1;
	}
	
	// ----- (A*He)(B*V) + C*He --> [(A+C)He](B*V) -----
	// Helium can interact with a mixed-species cluster so long as the sum of
	// the number of helium atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	
	// Get the index of the first HeV cluster in the reactants list
	
	for (int numVOther = 1; numVOther <= maxMixedClusterSize; numVOther++) {
		for (int numHeOther = 1; numVOther + numHeOther + numHe <=
			maxMixedClusterSize; numHeOther++) {
			std::map<std::string, int> speciesMap;
			speciesMap["He"] = numHeOther;
			speciesMap["V"] = numVOther;
			int indexOther = network->toClusterIndex(speciesMap);
			connectivityArray[indexOther] = 1;
		}
	}
	
	// (A*He)(B*I) + C*He --> ([A + C]*He)(B*I)
	
	for (int numIOther = 1; numIOther <= maxMixedClusterSize; numIOther++) {
		for (int numHeOther = 1; numIOther + numHeOther + numHe <=
			maxMixedClusterSize; numHeOther++) {
			
			std::map<std::string, int> speciesMap;
			speciesMap["He"] = numHeOther;
			speciesMap["I"] = numIOther;
			int indexOther = network->toClusterIndex(speciesMap);
			connectivityArray[indexOther] = 1;
		}
	}
	
	return connectivityArray;
}

std::vector<int> HeCluster::getDissociationConnectivity() {

	// Local Declarations
	int nReactants = network->reactants->size();
	std::vector<int> dissConnections(nReactants, 0);
	std::map<std::string, int> clusterMap;

	// He_x -> He_(x-1) + He, so a connection
	// to the Helium cluster with one helium,
	// and the Helium cluster with (x-1) helium
	clusterMap["He"] = size - 1;
	clusterMap["V"] = 0;
	clusterMap["I"] = 0;

	// For Helium dissociation make sure we have
	// more than one helium
	if (size != 1) {
		// He_x -> He_(x-1) + He
		dissConnections[network->toClusterIndex(clusterMap)] = 1;
		clusterMap["He"] = 1;
		dissConnections[network->toClusterIndex(clusterMap)] = 1;
	}

	// Trap Mutation...
	clusterMap["He"] = size; clusterMap["V"] = 1;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["He"] = 0; clusterMap["V"] = 0; clusterMap["I"] = 1;
	dissConnections[network->toClusterIndex(clusterMap)] = 1;

	// Return the connections
	return dissConnections;
}

bool HeCluster::isProductReactant(int reactantI, int reactantJ) {

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
	// total of size Helium, and a total of
	// 0 Vacancies
	return ((rI_I + rJ_I) == 0) && ((rI_He + rJ_He) == size)
			&& ((rI_V + rJ_V) == 0);
}

std::map<std::string, int> HeCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = size;
	clusterMap["V"] = 0;
	clusterMap["I"] = 0;

	// Return it
	return clusterMap;
}

double HeCluster::getReactionRadius() {
	double FourPi = 4.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3);
	double termOne = pow((3.0/FourPi)*(1.0/10.0)*aCubed*size,(1.0/3.0));
	double termTwo = pow((3.0/FourPi)*(1.0/10.0)*aCubed,(1.0/3.0));
	return .3 + termOne - termTwo;
}
