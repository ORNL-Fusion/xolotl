#include "InterstitialCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI) :
		PSICluster(nI) {
	// Set the reactant name appropriately
	name = "Interstitial";
}
InterstitialCluster::~InterstitialCluster() {
}

void InterstitialCluster::createReactionConnectivity() {

	// Local Declarations - Note the reference to the properties map
	std::map<std::string, std::string> props = *(network->properties);
	int numI = size, indexOther;
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int numVClusters = std::stoi(props["numVClusters"]);
	int maxIClusterSize = std::stoi(props["maxIClusterSize"]);
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	std::map<std::string, int> speciesMap;
	int totalSize = 1, firstSize = 0, secondSize = 0;
	int firstIndex = -1, secondIndex = -1;
	std::map<std::string, int> firstSpeciesMap, secondSpeciesMap;
	std::shared_ptr<Reactant> firstReactant, secondReactant;
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>reactants =
	network->reactants;

	/*
	 * This section fills the array of reacting pairs that combine to produce
	 * this cluster. The only reactions that produce I clusters are those I
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
		firstSpeciesMap["I"] = firstSize;
		secondSpeciesMap["I"] = secondSize;
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

	/* ----- (A*He) + (B*I) --> (A*He)*(B*I)
	 * Interstitials can interact with other interstitials, vacancies,
	 * helium, and mixed-species clusters. They cannot cluster with other
	 * interstitials that are so large that the combination of the two would
	 * produce an interstitial above the maximum size.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numHeOther = 1; numHeOther + numI <= maxMixedClusterSize;
			numHeOther++) {
		speciesMap["He"] = numHeOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- A*I + B*V -----
	 * → (A-B)*I, if A > B
	 * → (B-I)*V, if A < B
	 * → 0, if A = B
	 * Interstitial-Vacancy Annihilation
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numVOther = 1; numVOther <= numVClusters; numVOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["V"] = numVOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- A*I + B*I → (A+B)*I -----
	 *	Interstitial absorption
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numIOther = 1; numI + numIOther <= maxIClusterSize; numIOther++) {
		// Clear the map since we are reusing it
		speciesMap.clear();
		speciesMap["I"] = numIOther;
		int indexOther = network->toClusterIndex(speciesMap);
		reactionConnectivity[indexOther] = 1;
		combiningReactants.push_back(reactants->at(indexOther));
	}

	/* ----- (A*He)(B*V) + (C*I) --> (A*He)[(B-C)V] -----
	 * Interstitials interact with all mixed-species clusters by
	 * annihilating vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	for (int numVOther = 1; numVOther <= maxMixedClusterSize; numVOther++) {
		for (int numHeOther = 1; numVOther + numHeOther <= maxMixedClusterSize;
				numHeOther++) {
			// Clear the map since we are reusing it
			speciesMap.clear();
			bool connected = numVOther - numI >= 1;
			speciesMap["He"] = numHeOther;
			speciesMap["V"] = numVOther;
			int indexOther = network->toClusterIndex(speciesMap);
			if (indexOther >= reactants->size()) {
				break;
			}
			reactionConnectivity[indexOther] = (int) connected;
			combiningReactants.push_back(reactants->at(indexOther));
		}
	}

	/* ----- (A*He)*(B*I) + I --> (A*He)*(B + 1)*I -----
	 * Interstitial absorption by a He under the condition that (x + y + 1)
	 * <= maxSize
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	if (numI == 1 && numHeIClusters > 0) {
		for (int numIOther = 1; numIOther <= maxMixedClusterSize; numIOther++) {
			for (int numHeOther = 1;
					numIOther + numHeOther + 1 <= maxMixedClusterSize;
					numHeOther++) {
				// Clear the map since we are reusing it
				speciesMap.clear();
				speciesMap["He"] = numHeOther;
				speciesMap["I"] = numIOther;
				int indexOther = network->toClusterIndex(speciesMap);
				if (indexOther >= reactants->size()) {
					break;
				}
				reactionConnectivity[indexOther] = 1;
				combiningReactants.push_back(reactants->at(indexOther));
			}
		}
	}

	return;
}

void InterstitialCluster::createDissociationConnectivity() {

	// Commented out the below because it is wrong! FIXME!

	// Resize the connectivity row with zeroes
//	int reactantsLength = network->reactants->size();
//	dissociationConnectivity.resize(reactantsLength, 0);
}

bool InterstitialCluster::isProductReactant(int reactantI, int reactantJ) {

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

	// We should have this->size interstitials, a
	// total of 0 Helium, and a total of
	// 0 Vacancies
	return ((rI_I + rJ_I) == size) && ((rI_He + rJ_He) == 0)
			&& ((rI_V + rJ_V) == 0);
}

std::map<std::string, int> InterstitialCluster::getClusterMap() {
	// Local Declarations
	std::map<std::string, int> clusterMap;

	// Set the number of each species
	clusterMap["He"] = 0;
	clusterMap["V"] = 0;
	clusterMap["I"] = size;

	// Return it
	return clusterMap;
}

double InterstitialCluster::getReactionRadius() {

	double EightPi = 8.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::latticeConstant, 3.0);
	double termOne = 1.15 * (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant;
	double termTwo = pow((3.0 / EightPi) * aCubed * size, (1.0 / 3.0));
	double termThree = pow((3.0 / EightPi) * aCubed, (1.0 / 3.0));

	return termOne + termTwo - termThree;
}
