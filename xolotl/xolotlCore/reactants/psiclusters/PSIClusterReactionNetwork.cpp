
#include "PSIClusterReactionNetwork.h"

using namespace xolotlCore;


PSIClusterReactionNetwork::PSIClusterReactionNetwork()
: ReactionNetwork() {
}


PSIClusterReactionNetwork::PSIClusterReactionNetwork(
	const PSIClusterReactionNetwork &other)
: ReactionNetwork(other) {
}


std::map<std::string, int> PSIClusterReactionNetwork::toClusterMap(int index) {
	// This method defines the ordering and placement of all
	// reactions in a PSIClusterReactionNetwork
	
	std::map<std::string, int> clusterMap;
	
	// He clusters
	// size of `numHeClusters`
	
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	
	if (index < numHeClusters) {
		clusterMap["He"] = index + 1;
		return clusterMap;
	}
	
	index -= numHeClusters;
	
	// V (vacancy) clusters
	// size of `numVClusters`
	
	int numVClusters = std::stoi(properties->at("numVClusters"));
	
	if (index < numVClusters) {
		clusterMap["V"] = index + 1;
		return clusterMap;
	}
	
	index -= numVClusters;
	
	// I (interstitial) clusters
	// size of `numIClusters`
	
	int numIClusters = std::stoi(properties->at("numIClusters"));
	
	if (index < numIClusters) {
		clusterMap["I"] = index + 1;
		return clusterMap;
	}
	
	index -= numIClusters;
	
	// HeV clusters
	// size of `numHeVClusters`
	
	int numHeVClusters = std::stoi(properties->at("numHeVClusters"));
	int maxMixedClusterSize = std::stoi(properties->at("maxMixedClusterSize"));
	
	// Iterate through all posibilities of quantities for He and V
	// under the condition that
	// numHe + numV <= maxMixedClusterSize.
	
	// This doesn't need to be a loop, but it's less error prone than
	// a closed form.
	
	for (int numV = 1; numV <= maxMixedClusterSize; numV++) {
		int maxHe = maxMixedClusterSize - numV;
		int numHe = index + 1;
		
		if (numHe <= maxHe) {
			// The index must be referring to this quantity of He and V.
			// Finalize and return the cluster map
			
			clusterMap["He"] = numHe;
			clusterMap["V"] = numV;
			return clusterMap;
		}
		
		// Decrease the index by the number of vacancies
		// possible with the current quantity of helium
		index -= maxHe;
	}
	
	// HeI clusters
	// TODO
	
	throw std::string("HeI clusters cannot yet be indexed");
}


int PSIClusterReactionNetwork::toClusterIndex(std::map<std::string, int> clusterMap) {
	
	// If the clusterMap doesn't have one of these keys, the operator[] method
	// returns zero.
	
	int numHe = clusterMap["He"];
	int numV = clusterMap["V"];
	int numI = clusterMap["I"];
	
	// Convert the property strings so C++ can use them
	
	int numHeClusters = std::stoi(properties->at("numHeClusters"));
	int numVClusters = std::stoi(properties->at("numVClusters"));
	int numIClusters = std::stoi(properties->at("numIClusters"));
	int numHeVClusters = std::stoi(properties->at("numHeVClusters"));
	int maxMixedClusterSize = std::stoi(properties->at("maxMixedClusterSize"));
	
	int numSpecies = (numHe > 0) + (numV > 0) + (numI > 0);
	
	// Cluster should either contain one or two types of species.
	
	if (numSpecies == 0) {
		throw std::string("Cluster map contains no species");
	}
	else if (numSpecies == 1) {
		// Single species
		
		if (numHe) {
			return numHe - 1;
		}
		else if (numV) {
			return numV - numHeClusters - 1;
		}
		else if (numI) {
			return numI - numHeClusters - numVClusters - 1;
		}
	}
	else if (numSpecies == 2) {
		int indexOffset = numHeClusters + numVClusters + numIClusters;
		
		if (numHe && numV) {
			// HeVCluster
			
			// Closed form for converting a top-left triangle grid
			// to an index
			
			int index = (numV - 1) * maxMixedClusterSize -
				numV * (numV - 1) / 2 + numHe;
			
			return indexOffset + index;
		}
	}
	
	throw std::string("Reaction index could not be found");
}
