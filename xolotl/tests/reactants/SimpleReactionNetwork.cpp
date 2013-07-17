/*
 * SimpleReactionNetwork.cpp
 *
 *  Created on: Jun 29, 2013
 *      Author: bkj
 */

#include "SimpleReactionNetwork.h"
#include <PSICluster.h>
#include <HeCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <HeVCluster.h>
#include <HeInterstitialCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using std::shared_ptr;
using namespace xolotlCore;
using namespace testUtils;

SimpleReactionNetwork::SimpleReactionNetwork() {
	
	// Hard code the size of the largest cluster
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;

	// Add He clusters
	for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
		// Create a He cluster with cluster size numHe
		std::shared_ptr<HeCluster> cluster(new HeCluster(numHe));
		
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		// Create a He cluster with cluster size numV
		std::shared_ptr<VCluster> cluster(new VCluster(numV));
		
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		// Create a He cluster with cluster size numI
		std::shared_ptr<InterstitialCluster> cluster(
			new InterstitialCluster(numI));
		
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add HeV clusters, assuming that
	// numHe + numV <= maxMixedClusterSize
	
	int numHeVClusters = 0;
	
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			// Create a HeVCluster with the current amount of He and V
			std::shared_ptr<HeVCluster> cluster(new HeVCluster(numHe, numV));
			
			reactants->push_back(cluster);
			numHeVClusters++;
		}
	}
	
	// Create the HeI clusters in this simple reaction network
	int numHeIClusters = 0;

	// Create all possible combinations of numHe and numI
	// clusters with numHe, numI < maxClusterSize
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			// Create the HeI cluster
			std::shared_ptr<HeInterstitialCluster> cluster(new HeInterstitialCluster(numHe, numI));

			// Add it to the reactants vector
			reactants->push_back(cluster);

			// Increment the number of
			// HeIClusters for the properties map
			numHeIClusters++;
		}
	}

	// Setup the properties map
	(*properties)["maxHeClusterSize"] = std::to_string((long long) maxClusterSize);
	(*properties)["maxVClusterSize"] = std::to_string((long long) maxClusterSize);
	(*properties)["maxIClusterSize"] = std::to_string((long long) maxClusterSize);
	(*properties)["maxMixedClusterSize"] = std::to_string((long long) maxClusterSize);
	
	(*properties)["numHeClusters"] = std::to_string((long long) numClusters);
	(*properties)["numVClusters"] = std::to_string((long long) numClusters);
	(*properties)["numIClusters"] = std::to_string((long long) numClusters);
	(*properties)["numHeVClusters"] = std::to_string((long long) numHeVClusters);
	(*properties)["numHeIClusters"] = std::to_string((long long) numHeIClusters);
}

SimpleReactionNetwork::~SimpleReactionNetwork() {
	// Nothing to do
}

/**
 * This operation creates a SimpleReactionNetwork and makes sure that it is
 * properly registered with the clusters it contains. This operation should
 * always be called instead of constructing a SimpleReactionNetwork manually.
 * @return The reaction network.
 */
std::shared_ptr<xolotlCore::ReactionNetwork> testUtils::getSimpleReactionNetwork() {
	// Create the network
	std::shared_ptr<xolotlCore::ReactionNetwork> network(new SimpleReactionNetwork());
	// Register the reaction network with its clusters
	for (int i = 0; i < network->reactants->size(); i++) {
		network->reactants->at(i)->setReactionNetwork(network);
	}
	
	
	// TEMPORARY
	// Print the reaction connectivity matrix
	
	for (auto reactantIt = network->reactants->begin();
		reactantIt != network->reactants->end(); reactantIt++)
	{
		std::shared_ptr<PSICluster> cluster =
			std::dynamic_pointer_cast<PSICluster>(*reactantIt);
		std::vector<int> conn = cluster->getReactionConnectivity();
		
		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++)
		{
			printf("%s", *connIt ? "* " : "  ");
		}
		printf("\n");
	}
	
	return network;
}
