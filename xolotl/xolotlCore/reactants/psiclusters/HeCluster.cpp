// Includes
#include "HeCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

HeCluster::HeCluster(int nHe) :
		PSICluster(nHe) {
	// Set the reactant name appropriately
	name = "He";
}

HeCluster::~HeCluster() {
}

std::shared_ptr<Reactant> HeCluster::clone() {
	std::shared_ptr<Reactant> reactant(new HeCluster(*this));
	return reactant;
}

void HeCluster::createReactionConnectivity() {

	// Local Declarations - Note the reference to the properties map
	std::map<std::string, std::string> props = *(network->properties);
	int indexOther, otherNumHe, otherNumV, otherNumI;
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	std::map<std::string, int> speciesMap, otherSpeciesMap;
	int totalSize = 1, firstSize = 0, secondSize = 0;
	int firstIndex = -1, secondIndex = -1;
	std::map<std::string, int> firstSpeciesMap, secondSpeciesMap;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>>reactants =
			network->reactants;

	/*
	 * This section fills the array of reacting pairs that combine to produce
	 * this cluster. The only reactions that produce He clusters are those He
	 * clusters that are smaller than this.size. Each cluster i combines with
	 * a second cluster of size this.size - i.size.
	 *
	 * Total size starts with a value of one so that clusters of size one are
	 * not considered in this loop.
	 */
	while (totalSize < size) {
		// Increment the base sizes
		++firstSize;
		secondSize = size - firstSize;
		// Update the maps
		firstSpeciesMap["He"] = firstSize;
		secondSpeciesMap["He"] = secondSize;
		// Get the first and second reactants for the reaction
		// first + second = this.
		firstIndex = network->toClusterIndex(firstSpeciesMap);
		firstReactant = reactants->at(firstIndex);
		secondIndex = network->toClusterIndex(secondSpeciesMap);
		secondReactant = reactants->at(secondIndex);
		// Create a ReactingPair with the two reactants
		ReactingPair pair;
		pair.first = std::dynamic_pointer_cast < PSICluster > (firstReactant);
		pair.second = std::dynamic_pointer_cast < PSICluster > (secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Update the total size. Do not delete this or you'll have an infinite
		// loop!
		totalSize = firstSize + secondSize;
	}

	/* ----- A*He + B*He --> (A+B)*He -----
	 * This cluster should interact with all other clusters of the same type up
	 * to the max size minus the size of this one to produce larger clusters.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numHeOther = 1; size + numHeOther <= maxHeClusterSize;
			numHeOther++) {
		speciesMap["He"] = numHeOther;
		indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* -----  A*He + B*V --> (A*He)(B*V) -----
	 * Helium clusters can interact with any vacancy cluster so long as the sum
	 * of the number of helium atoms and vacancies does not produce a cluster
	 * with a size greater than the maximum mixed-species cluster size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numVOther = 1; size + numVOther <= maxMixedClusterSize;
			numVOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["V"] = numVOther;
		indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- A*He + B*I --> (A*He)(B*I)
	 * Helium clusters can interact with any interstitial cluster so long as
	 * the sum of the number of helium atoms and interstitials does not produce
	 * a cluster with a size greater than the maximum mixed-species cluster
	 * size.

	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numIOther = 1; numIOther + size <= maxMixedClusterSize;
			numIOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["I"] = numIOther;
		indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- (A*He)(B*V) + C*He --> [(A+C)He](B*V) -----
	 * Helium can interact with a mixed-species cluster so long as the sum of
	 * the number of helium atoms and the size of the mixed-species cluster
	 * does not exceed the maximum mixed-species cluster size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 *
	 * Find the clusters by looping over all size combinations of HeV clusters.
	 */
	if (numHeVClusters > 0) {
		// Get the index of the first HeV cluster. Clear the map since we are
		// reusing it.
		speciesMap.clear();
		speciesMap["He"] = 1;
		speciesMap["V"] = 1;
		int heVIndex = network->toClusterIndex(speciesMap);
		// Connect to the HeV clusters if possible
		connectWithMixedClusters(heVIndex,heVIndex+numHeVClusters,"V");
	}

	/* ----- (A*He)(B*I) + C*He --> ([A + C]*He)(B*I) -----
	 * Helium-interstitial clusters can absorb single-species helium clusters
	 * so long as the maximum cluster size limit is not violated.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 *
	 * Find the clusters by looping over all size combinations of HeI clusters.
	 */
	if (numHeIClusters > 0) {
		// Get the index of the first HeI cluster. Clear the map since we are
		// reusing it.
		speciesMap.clear();
		speciesMap["He"] = 1;
		speciesMap["I"] = 1;
		int heIIndex = network->toClusterIndex(speciesMap);
		// Connect to the HeV clusters if possible
		connectWithMixedClusters(heIIndex,heIIndex+numHeIClusters,"V");
	}

	return;
}

void HeCluster::connectWithMixedClusters(int startIndex, int stopIndex,
		std::string mixedSpecies) {

	// Local Declarations
	int otherNumHe, otherNumMixed;
	std::shared_ptr < std::vector<std::shared_ptr<Reactant>>>reactants =
			network->reactants;
	int maxMixedClusterSize = std::stoi(network->properties->at("maxMixedClusterSize"));
	std::map<std::string, int> otherSpeciesMap;
	std::shared_ptr<PSICluster> mixedCluster;

	// Loop over the mixed clusters
	for (int i = startIndex; i < stopIndex; i++) {
		mixedCluster = std::dynamic_pointer_cast < PSICluster
				> (reactants->at(i));
		// Get the cluster sizes for the mixed cluster
		otherSpeciesMap = (std::dynamic_pointer_cast < PSICluster
				> (reactants->at(i)))->getClusterMap();
		otherNumHe = otherSpeciesMap["He"];
		otherNumMixed = otherSpeciesMap[mixedSpecies];
		// React with it if the sizes are compatible.
		if (otherNumHe + otherNumMixed + size <= maxMixedClusterSize) {
			reactionConnectivity[i] = 1;
			combiningReactants.push_back(reactants->at(i));
		}
	}

	return;
}

void HeCluster::createDissociationConnectivity() {

	// Local Declarations
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
		dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
		clusterMap["He"] = 1;
		dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	}

	// Trap Mutation...
	clusterMap["He"] = size;
	clusterMap["V"] = 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
	clusterMap["He"] = 0;
	clusterMap["V"] = 0;
	clusterMap["I"] = 1;
	dissociationConnectivity[network->toClusterIndex(clusterMap)] = 1;
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

const std::map<std::string,int> HeCluster::getComposition() {
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
	double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
			(1.0 / 3.0));
	double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
	return 0.3 + termOne - termTwo;
}
