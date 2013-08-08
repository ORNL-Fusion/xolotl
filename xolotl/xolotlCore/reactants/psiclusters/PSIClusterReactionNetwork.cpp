#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include <iostream>
#include <algorithm>

using namespace xolotlCore;
using std::shared_ptr;

PSIClusterReactionNetwork::PSIClusterReactionNetwork() :
		ReactionNetwork() {

	// Initialize default properties
	(*properties)["reactionsEnabled"] = "true";
	(*properties)["dissociationsEnabled"] = "true";
}

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		const PSIClusterReactionNetwork &other) :
		ReactionNetwork(other) {
}

std::map<std::string, int> PSIClusterReactionNetwork::toClusterMap(
		int index) const {

	// This functionality is implemented in PSICluster.
	// Instead of redefining the map for all indices, we simply obtain
	// the cluster map from the reactant in the reactants vector at the
	// position of the index.
	shared_ptr<PSICluster> cluster = std::dynamic_pointer_cast < PSICluster
			> (reactants->at(index));
	return cluster->getClusterMap();
}

int PSIClusterReactionNetwork::toClusterIndex(
		std::map<std::string, int> clusterMap) const {

	// If the clusterMap doesn't have one of the following keys,
	// the value will be set to zero.
	int numHe = 0;
	int numV = 0;
	int numI = 0;
	int numReactants = reactants->size();
	int finalIndex = 0;
	// A handy typedef
	typedef const std::map<std::string, int> clusterMap_t;
	// Convert the property strings so we can use them
	int numHeClusters = std::stoi((*properties)["numHeClusters"]);
	int numVClusters = std::stoi((*properties)["numVClusters"]);
	int numIClusters = std::stoi((*properties)["numIClusters"]);
	int numHeVClusters = std::stoi((*properties)["numHeVClusters"]);
	int numHeIClusters = std::stoi((*properties)["numHeIClusters"]);
	int maxMixedClusterSize = std::stoi((*properties)["maxMixedClusterSize"]);
	int numSpecies = 0;

	// Update the indices
	numHe = std::max(clusterMap["He"],0);
	numV = std::max(clusterMap["V"],0);
	numI = std::max(clusterMap["I"],0);
	numSpecies = (numHe > 0) + (numV > 0) + (numI > 0);

	// Cluster should either contain one or two types of species.
	if (numSpecies == 0) {
		throw std::string("Cluster map contains no species");
	} else if (numSpecies == 1) {
		// Single species
		if (numHe) {
			finalIndex = numHe - 1;
		} else if (numV) {
			finalIndex = numV + numHeClusters - 1;
		} else if (numI) {
			finalIndex = numI + numHeClusters + numVClusters - 1;
		}
		return finalIndex;
	} else if (numSpecies == 2) {
		// HeVCluster
		int indexOffset = numHeClusters + numVClusters + numIClusters;
		if (numHe && numV && numHeVClusters > 0) {
			// Closed form for converting a top-left triangle grid
			// to an index
			int index = (numV - 1) * maxMixedClusterSize - numV * (numV - 1) / 2
					+ numHe - 1;
			std::cout << indexOffset << " " << index << " " << maxMixedClusterSize << " " << numHe << " " << numV << std::endl;
			finalIndex = indexOffset + index;
		}
		// ----- HeICluster -----
		// Increment the offset by the number of HeVClusters
		// indexOffset += maxMixedClusterSize * (maxMixedClusterSize - 1) / 2;
		else if (numHe && numI && numHeIClusters > 0) {
			indexOffset += numHeVClusters;
			// Closed form for converting a top-left triangle grid
			// to an index
			int index = (numI - 1) * maxMixedClusterSize - numI * (numI - 1) / 2
					+ numHe - 1;
			finalIndex = indexOffset + index;
		}
		return finalIndex;
	}

	throw std::string("Reaction index could not be found");
}
