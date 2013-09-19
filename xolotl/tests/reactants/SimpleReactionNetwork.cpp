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
#include <iostream>

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
		add(cluster);
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		// Create a He cluster with cluster size numV
		std::shared_ptr<VCluster> cluster(new VCluster(numV));
		// Add it to the network
		add(cluster);
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		// Create a He cluster with cluster size numI
		std::shared_ptr<InterstitialCluster> cluster(
				new InterstitialCluster(numI));
		// Add it to the network
		add(cluster);
	}

	// Add HeV clusters, assuming that
	// numHe + numV <= maxMixedClusterSize
	int numHeVClusters = 0;
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			// Create a HeVCluster with the current amount of He and V
			std::shared_ptr<HeVCluster> cluster(new HeVCluster(numHe, numV));
			add(cluster);
		}
	}

	// Add HeI clusters
	int numHeIClusters = 0;
	// Create all possible combinations of numHe and numI
	// clusters with numHe, numI < maxClusterSize
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			// Create the HeI cluster
			std::shared_ptr<HeInterstitialCluster> cluster(
					new HeInterstitialCluster(numHe, numI));
			// Add it to the reactants vector
			add(cluster);
		}
	}

	return;
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
	std::shared_ptr<xolotlCore::ReactionNetwork> network(
			new SimpleReactionNetwork());
	std::cout << "SimpleReactionNetwork Message: "
			<< "Created network with size " << network->size() << std::endl;
	// Register the reaction network with its clusters
	auto reactants = network->getAll();
	for (int i = 0; i < reactants->size(); i++) {
		reactants->at(i)->setReactionNetwork(network);
	}

	// ----- TEMPORARY DEBUG OUTPUT!!!!! -----
	// Print the reaction connectivity matrix
	for (auto reactantIt = reactants->begin();
			reactantIt != reactants->end(); reactantIt++) {
		std::shared_ptr<PSICluster> cluster = std::dynamic_pointer_cast<
				PSICluster>(*reactantIt);
		std::vector<int> conn = cluster->getConnectivity();

		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++) {
			printf("%s", *connIt ? "* " : "  ");
		}
		printf("\n");
	}

	return network;
}
