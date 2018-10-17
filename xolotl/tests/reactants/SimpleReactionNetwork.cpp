/*
 * SimpleReactionNetwork.cpp
 *
 *  Created on: Jun 29, 2013
 *      Author: bkj
 */

#include "SimpleReactionNetwork.h"
#include <PSICluster.h>
#include <PSIHeCluster.h>
#include <PSIDCluster.h>
#include <PSITCluster.h>
#include <PSIVCluster.h>
#include <PSIInterstitialCluster.h>
#include <PSIMixedCluster.h>
#include <PSIHeInterstitialCluster.h>
#include <NECluster.h>
#include <NEXeCluster.h>
#include <FeCluster.h>
#include <FeHeCluster.h>
#include <FeVCluster.h>
#include <FeInterstitialCluster.h>
#include <FeHeVCluster.h>
#include <AlloyCluster.h>
#include <AlloyVacCluster.h>
#include <AlloyIntCluster.h>
#include <AlloyVoidCluster.h>
#include <AlloyFaultedCluster.h>
#include <AlloyFrankCluster.h>
#include <AlloyPerfectCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;
using namespace xolotlPerf;

SimplePSIReactionNetwork::SimplePSIReactionNetwork(const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSIClusterReactionNetwork(registry) {
	// Add He clusters
	for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
		// Create a He cluster with cluster size numHe
		auto cluster = new PSIHeCluster(numHe, *this, registry);
		// Set the diffusion factor for some of them to 1.0 so that they can react
		if (numHe < 8)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<PSICluster>(cluster));
	}

	// Add deuterium clusters
	for (int numD = 1; numD <= maxClusterSize; numD++) {
		// Create a D cluster with cluster size numD
		auto cluster = new PSIDCluster(numD, *this, registry);
		// Set the diffusion factor for the first one so that it can react
		if (numD == 1)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<PSICluster>(cluster));
	}

	// Add tritium clusters
	for (int numT = 1; numT <= maxClusterSize; numT++) {
		// Create a T cluster with cluster size numT
		auto cluster = new PSITCluster(numT, *this, registry);
		// Set the diffusion factor for the first one so that it can react
		if (numT == 1)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<PSICluster>(cluster));
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		// Create a V cluster with cluster size numV
		auto cluster = new PSIVCluster(numV, *this, registry);
		// Set the diffusion factor for the first one so that it can react
		if (numV == 1)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<PSICluster>(cluster));
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		// Create a He cluster with cluster size numI
		auto cluster = new PSIInterstitialCluster(numI, *this, registry);
		// Set the diffusion factor for all of them to 1.0 so that they can react
		cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<PSICluster>(cluster));
	}

	// Add HeV clusters, assuming that
	// numHe + numV <= maxMixedClusterSize
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe + numV <= maxClusterSize; numHe++) {
			for (int numD = 0; numD <= 1; numD++) {
				for (int numT = 0; numT <= 1; numT++) {
					// Create a MixedCluster with the current amount of He and V
					auto cluster = new PSIMixedCluster(numHe, numD, numT, numV,
							*this, registry);
					add(std::unique_ptr<PSICluster>(cluster));
				}
			}
		}
	}

	// Add HeI clusters
	// Create all possible combinations of numHe and numI
	// clusters with numHe, numI < maxClusterSize
	for (int numI = 1; numI <= maxClusterSize; numI++) {
		for (int numHe = 1; numHe + numI <= maxClusterSize; numHe++) {
			// Create the HeI cluster
			auto cluster = new PSIHeInterstitialCluster(numHe, numI, *this,
					registry);
			// Add it to the reactants vector
			add(std::unique_ptr<PSICluster>(cluster));
		}
	}

	return;
}

SimpleNEReactionNetwork::SimpleNEReactionNetwork(const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		NEClusterReactionNetwork(registry) {
	// Add Xe clusters
	for (int numXe = 1; numXe <= maxClusterSize; numXe++) {
		// Create a He cluster with cluster size numHe
		auto cluster = new NEXeCluster(numXe, *this, registry);
		// Set the diffusion factor for the first one so that it can react
		if (numXe == 1)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<NECluster>(cluster));
	}

	return;
}

SimpleFeReactionNetwork::SimpleFeReactionNetwork(const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		FeClusterReactionNetwork(registry) {
	// Add He clusters
	for (int numHe = 1; numHe <= 8; numHe++) {
		// Create a He cluster with cluster size numHe
		auto cluster = new FeHeCluster(numHe, *this, registry);
		// Set the diffusion factor for some of them to 1.0 so that they can react
		if (numHe < 8)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<FeCluster>(cluster));
	}

	// Add vacancy clusters
	for (int numV = 1; numV <= 9; numV++) {
		// Create a V cluster with cluster size numV
		auto cluster = new FeVCluster(numV, *this, registry);
		// Set the diffusion factor for the first one so that it can react
		if (numV == 1)
			cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<FeCluster>(cluster));
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= 1; numI++) {
		// Create a He cluster with cluster size numI
		auto cluster = new FeInterstitialCluster(numI, *this, registry);
		// Set the diffusion factor for all of them to 1.0 so that they can react
		cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<FeCluster>(cluster));
	}

	// Add HeV clusters, assuming that
	// numHe + numV <= maxMixedClusterSize
	for (int numV = 1; numV <= maxClusterSize; numV++) {
		for (int numHe = 1; numHe <= maxClusterSize; numHe++) {
			// Create a HeVCluster with the current amount of He and V
			auto cluster = new FeHeVCluster(numHe, numV, *this, registry);
			add(std::unique_ptr<FeCluster>(cluster));
		}
	}

	return;
}

SimpleAlloyReactionNetwork::SimpleAlloyReactionNetwork(const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		AlloyClusterReactionNetwork(registry) {
	// Add vacancy clusters
	for (int numV = 1; numV <= 5; numV++) {
		// Create a vacancy cluster with cluster size numV
		auto cluster = new AlloyVacCluster(numV, *this, registry);
		// Set the diffusion factor for some of them to 1.0 so that they can react
		cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<AlloyCluster>(cluster));
	}

	// Add interstitial clusters
	for (int numI = 1; numI <= 4; numI++) {
		// Create an interstitial cluster with cluster size numI
		auto cluster = new AlloyIntCluster(numI, *this, registry);
		// Set the diffusion factor for all of them to 1.0 so that they can react
		cluster->setDiffusionFactor(1.0);
		// Add it to the network
		add(std::unique_ptr<AlloyCluster>(cluster));
	}

	// Add void clusters
	for (int numV = 6; numV <= maxClusterSize; numV++) {
		// Create a void cluster with cluster size numV
		auto cluster = new AlloyVoidCluster(numV, *this, registry);
		// Add it to the network
		add(std::unique_ptr<AlloyCluster>(cluster));
	}

	// Add faulted clusters
	for (int numV = 6; numV <= maxClusterSize; numV++) {
		// Create a faulted cluster with cluster size numV
		auto cluster = new AlloyFaultedCluster(numV, *this, registry);
		// Add it to the network
		add(std::unique_ptr<AlloyCluster>(cluster));
	}

	// Add perfect clusters
	for (int numI = 5; numI <= maxClusterSize; numI++) {
		// Create a perfect cluster with cluster size numI
		auto cluster = new AlloyPerfectCluster(numI, *this, registry);
		// Add it to the network
		add(std::unique_ptr<AlloyCluster>(cluster));
	}

	// Add frank clusters
	for (int numI = 5; numI <= maxClusterSize; numI++) {
		// Create a frank cluster with cluster size numI
		auto cluster = new AlloyFrankCluster(numI, *this, registry);
		// Add it to the network
		add(std::unique_ptr<AlloyCluster>(cluster));
	}

	return;
}

shared_ptr<xolotlCore::PSIClusterReactionNetwork> testUtils::getSimplePSIReactionNetwork(
		const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	// Create the network
	shared_ptr<xolotlCore::PSIClusterReactionNetwork> network(
			new SimplePSIReactionNetwork(maxClusterSize, registry));
	cout << "SimpleReactionNetwork Message: " << "Created network with size "
			<< network->size() << endl;

	// Update reactants now that they are in network.
	auto reactants = network->getAll();
	for (IReactant& currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Define the phase space for the network
	int nDim = 1;
	Array<int, 5> list;
	list[0] = 0;
	// Now that all the clusters are created
	// Give the information on the phase space to the network
	network->setPhaseSpace(nDim, list);

	// Create the reactions
	network->createReactionConnectivity();
	// Recompute Ids and network size
	network->reinitializeNetwork();
	// Redefine the connectivities
	network->reinitializeConnectivities();

	return network;
}

shared_ptr<xolotlCore::NEClusterReactionNetwork> testUtils::getSimpleNEReactionNetwork(
		const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	// Create the network
	shared_ptr<xolotlCore::NEClusterReactionNetwork> network(
			new SimpleNEReactionNetwork(maxClusterSize, registry));
	cout << "SimpleReactionNetwork Message: " << "Created network with size "
			<< network->size() << endl;

	// Update reactants now that they are in network.
	auto reactants = network->getAll();
	for (IReactant& currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Create the reactions
	network->createReactionConnectivity();
	// Recompute Ids and network size
	network->reinitializeNetwork();
	// Redefine the connectivities
	network->reinitializeConnectivities();

	// ----- TEMPORARY DEBUG OUTPUT!!!!! -----
	// Print the reaction connectivity matrix
	for (IReactant& currCluster : reactants) {
		vector<int> conn = currCluster.getConnectivity();

		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++) {
			printf("%s", *connIt ? "* " : "' ");
		}
		printf("\n");
	}

	return network;
}

shared_ptr<xolotlCore::FeClusterReactionNetwork> testUtils::getSimpleFeReactionNetwork(
		const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	// Create the network
	shared_ptr<xolotlCore::FeClusterReactionNetwork> network(
			new SimpleFeReactionNetwork(maxClusterSize, registry));
	cout << "SimpleReactionNetwork Message: " << "Created network with size "
			<< network->size() << endl;

	// Update reactants now that they are in network.
	auto reactants = network->getAll();
	for (IReactant& currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Create the reactions
	network->createReactionConnectivity();
	// Recompute Ids and network size
	network->reinitializeNetwork();
	// Redefine the connectivities
	network->reinitializeConnectivities();

	// ----- TEMPORARY DEBUG OUTPUT!!!!! -----
	// Print the reaction connectivity matrix
	for (IReactant& currCluster : reactants) {
		vector<int> conn = currCluster.getConnectivity();

		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++) {
			printf("%s", *connIt ? "* " : "' ");
		}
		printf("\n");
	}

	return network;
}

shared_ptr<xolotlCore::AlloyClusterReactionNetwork> testUtils::getSimpleAlloyReactionNetwork(
		const int maxClusterSize,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	// Create the network
	shared_ptr<xolotlCore::AlloyClusterReactionNetwork> network(
			new SimpleAlloyReactionNetwork(maxClusterSize, registry));
	cout << "SimpleReactionNetwork Message: " << "Created network with size "
			<< network->size() << endl;

	// Update reactants now that they are in network.
	auto reactants = network->getAll();
	for (IReactant& currCluster : reactants) {
		currCluster.updateFromNetwork();
	}

	// Create the reactions
	network->createReactionConnectivity();
	// Recompute Ids and network size
	network->reinitializeNetwork();
	// Redefine the connectivities
	network->reinitializeConnectivities();

	// ----- TEMPORARY DEBUG OUTPUT!!!!! -----
	// Print the reaction connectivity matrix
	for (IReactant& currCluster : reactants) {
		vector<int> conn = currCluster.getConnectivity();

		for (auto connIt = conn.begin(); connIt != conn.end(); connIt++) {
			printf("%s", *connIt ? "* " : "' ");
		}
		printf("\n");
	}

	return network;
}
