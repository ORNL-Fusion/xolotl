
#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"

using namespace xolotlCore;
using std::shared_ptr;


PSIClusterReactionNetwork::PSIClusterReactionNetwork()
: ReactionNetwork() {
}


PSIClusterReactionNetwork::PSIClusterReactionNetwork(
	const PSIClusterReactionNetwork &other)
: ReactionNetwork(other) {
}


std::map<std::string, int> PSIClusterReactionNetwork::toClusterMap(int index) {
	
	// This functionality is implemented in PSICluster.
	// Instead of redefining the map for all indices, we simply obtain
	// the cluster map from the reactant in the reactants vector at the
	// position of the index.
	
	shared_ptr<PSICluster> cluster =
		std::dynamic_pointer_cast<PSICluster>(reactants->at(index));
	return cluster->getClusterMap();
}


int PSIClusterReactionNetwork::toClusterIndex(std::map<std::string, int> clusterMap) {
	
	// If the clusterMap doesn't have one of these keys, the operator[] method
	// returns zero.
	
	int numHe = clusterMap["He"];
	int numV = clusterMap["V"];
	int numI = clusterMap["I"];
	
	// Convert the property strings so C++ can use them
	
	std::map<std::string, std::string> &props = *properties;
	int numHeClusters = std::stoi(props["numHeClusters"]);
	int numVClusters = std::stoi(props["numVClusters"]);
	int numIClusters = std::stoi(props["numIClusters"]);
	int numHeVClusters = std::stoi(props["numHeVClusters"]);
	int numHeIClusters = std::stoi(props["numHeIClusters"]);
	int maxMixedClusterSize = std::stoi(props["maxMixedClusterSize"]);
	
	int numSpecies = (numHe > 0) + (numV > 0) + (numI > 0);
	
	// Cluster should either contain one or two types of species.
	
	if (numSpecies == 0) {
		throw std::string("Cluster map contains no species");
	} else if (numSpecies == 1) {
		// Single species
		
		if (numHe) {
			return numHe - 1;
		}
		else if (numV) {
			return numV + numHeClusters - 1;
		}
		else if (numI) {
			return numI + numHeClusters + numVClusters - 1;
		}
	} else if (numSpecies == 2) {
		int indexOffset = numHeClusters + numVClusters + numIClusters;
		
		if (numHe && numV) {
			// HeVCluster
			
			// Closed form for converting a top-left triangle grid
			// to an index
			
			int index = (numV - 1) * maxMixedClusterSize -
				numV * (numV - 1) / 2 + numHe - 1;
			
			return indexOffset + index;
		} else if (numHe && numI) {
			// FIXME... Andrew
			return 0;
		}
	}
	
	throw std::string("Reaction index could not be found");
}
