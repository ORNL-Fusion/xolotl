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
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace xolotlCore;
using namespace testUtils;

SimpleReactionNetwork::SimpleReactionNetwork() :
		ReactionNetwork() {

	// Local Declarations
	int numHe = 1, numV = 1;
	int maxClusterSize = 10, numClusters = maxClusterSize;

	// Fill the ReactionNetwork with 10 He clusters
	for (int i = 0; i < numClusters; i++) {
		// Create a He cluster with cluster size i
		std::shared_ptr<HeCluster> cluster(new HeCluster(i + 1));
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add vacancy clusters
	for (int i = numClusters; i < 2 * numClusters; i++) {
		// Create a He cluster with cluster size i
		std::shared_ptr<VCluster> cluster(new VCluster(i + 1 - numClusters));
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add interstitial clusters
	for (int i = 2 * numClusters; i < 3 * numClusters; i++) {
		// Create a He cluster with cluster size i
		std::shared_ptr<InterstitialCluster> cluster(
				new InterstitialCluster(i + 1 - 2 * numClusters));
		// Add it to the network
		reactants->push_back(cluster);
	}

	// Add mixed-species clusters.
	for (int i = 3 * numClusters; i < 3 * numClusters + 2 * maxClusterSize;
			i++) {
		// Create and configure the species map
		std::map<std::string, int> speciesMap;
		speciesMap["He"] = numHe;
		speciesMap["V"] = numV;
		// Create a He cluster with cluster size i
		std::shared_ptr<HeVCluster> cluster(
				new HeVCluster(speciesMap));
		// Add it to the network
		reactants->push_back(cluster);
		// Figure out which species to increment for the next iteration
		if (numHe <= numV && numHe < maxClusterSize)
			numHe++;
		else if (numV < maxClusterSize)
			numV++;
		else
			break;
	}

	// Setup the properties map
	(*properties)["maxHeClusterSize"] = std::to_string((long long) maxClusterSize);
	(*properties)["maxVClusterSize"] = std::to_string((long long) maxClusterSize);
	(*properties)["maxIClusterSize"] = std::to_string((long long) maxClusterSize);
	(*properties)["maxHeVClusterSize"] = std::to_string(
			(long long) 2 * maxClusterSize);
	(*properties)["numHeClusters"] = std::to_string((long long) numClusters);
	(*properties)["numVClusters"] = std::to_string((long long) numClusters);
	(*properties)["numIClusters"] = std::to_string((long long) numClusters);
	(*properties)["numHeVClusters"] = std::to_string(
			(long long) (numClusters * numClusters));

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
	return network;
}
