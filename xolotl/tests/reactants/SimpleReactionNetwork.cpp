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
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;
using namespace xolotlPerf;

SimpleReactionNetwork::SimpleReactionNetwork(const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
	PSIClusterReactionNetwork(registry){
	// Add He clusters
	for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
		// Create a He cluster with cluster size numHe
		shared_ptr<HeCluster> cluster = make_shared<HeCluster>(numHe, registry);
		// Add it to the network
		add(cluster);
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		// Create a V cluster with cluster size numV
		shared_ptr<VCluster> cluster = make_shared<VCluster>(numV, registry);
		// Add it to the network
		add(cluster);
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		// Create a He cluster with cluster size numI
		shared_ptr<InterstitialCluster> cluster = make_shared<InterstitialCluster>(numI, registry);
		// Add it to the network
		add(cluster);
	}

	// Add HeV clusters, assuming that
	// numHe + numV <= maxMixedClusterSize
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			// Create a HeVCluster with the current amount of He and V
			shared_ptr<HeVCluster> cluster = make_shared<HeVCluster>(numHe, numV, registry);
			add(cluster);
		}
	}

	// Add HeI clusters
	// Create all possible combinations of numHe and numI
	// clusters with numHe, numI < maxClusterSize
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			// Create the HeI cluster
			shared_ptr<HeInterstitialCluster> cluster = make_shared<HeInterstitialCluster>(numHe, numI, registry);
			// Add it to the reactants vector
			add(cluster);
		}
	}

	return;
}

shared_ptr<xolotlCore::ReactionNetwork> testUtils::getSimpleReactionNetwork(const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	// Create the network
	shared_ptr<xolotlCore::ReactionNetwork> network(
			new SimpleReactionNetwork(maxClusterSize, registry));
	cout << "SimpleReactionNetwork Message: "
			<< "Created network with size " << network->size() << endl;
	// Register the reaction network with its clusters
	auto reactants = network->getAll();
	for (int i = 0; i < reactants->size(); i++) {
		reactants->at(i)->setReactionNetwork(network);
	}

	// ----- TEMPORARY DEBUG OUTPUT!!!!! -----
	// Print the reaction connectivity matrix
	for (auto reactantIt = reactants->begin();
			reactantIt != reactants->end(); reactantIt++) {
		auto cluster = (PSICluster *) (*reactantIt);
		vector<int> conn = cluster->getConnectivity();

		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++) {
			printf("%s", *connIt ? "* " : "' ");
		}
		printf("\n");
	}

	return network;
}
