/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include <HeCluster.h>
#include <VCluster.h>
#include <InterstitialCluster.h>
#include <MixedSpeciesCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>

using namespace std;
using namespace xolotlCore;

/**
 * This operation writes information about the cluster to stdout
 * @param cluster The cluster to dump
 */
void writeCluster(shared_ptr<Reactant> cluster) {
	shared_ptr<PSICluster> psiCluster = static_pointer_cast < PSICluster
			> (cluster);
	BOOST_TEST_MESSAGE(psiCluster->getSize());
	return;
}

/**
 * This suite is responsible for testing the HeCluster.
 */BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/** This operation checks the ability of the HeCluster to describe its connectivity to other clusters. */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	// Local Declarations
	int maxClusterSize = 10, numClusters = maxClusterSize;
	vector<int> connectivityArray;
	shared_ptr<ReactionNetwork> network(new ReactionNetwork());
	shared_ptr < vector<shared_ptr<Reactant> > > reactants = network->reactants;
	shared_ptr<std::map<std::string, std::string>> props = network->properties;

	// Fill the ReactionNetwork with 10 He clusters
	for (int i = 0; i < numClusters; i++) {
		// Create a He cluster with cluster size i
		shared_ptr<HeCluster> cluster(new HeCluster(i+1));
		// Add it to the network
		reactants->push_back(cluster);
		// Register the network with the cluster
		cluster->setReactionNetwork(network);
	}

	// Add vacancy clusters
	for (int i = numClusters; i < 2*numClusters - 1; i++) {
		// Create a He cluster with cluster size i
		shared_ptr<VCluster> cluster(new VCluster(i+1-numClusters));
		// Add it to the network
		reactants->push_back(cluster);
		// Register the network with the cluster
		cluster->setReactionNetwork(network);
	}

	// Add interstitial clusters
	for (int i = 2*numClusters; i < 3*numClusters - 1; i++) {
		// Create a He cluster with cluster size i
		shared_ptr<InterstitialCluster> cluster(new InterstitialCluster(i+1-2*numClusters));
		// Add it to the network
		reactants->push_back(cluster);
		// Register the network with the cluster
		cluster->setReactionNetwork(network);
	}

	// Add mixed-species clusters -- FIX SIZES!
//	for (int i = 3*numClusters; i < 3*numClusters; i++) {
//		// Create a He cluster with cluster size i
//		shared_ptr<MixedSpeciesCluster> cluster(new MixedSpeciesCluster(i+1));
//		// Add it to the network
//		reactants->push_back(cluster);
//		// Register the network with the cluster
//		cluster->setReactionNetwork(network);
//	}

	// Setup the properties map
	(*props)["maxHeClusterSize"] = std::to_string((long long) maxClusterSize);
	(*props)["maxVClusterSize"] = std::to_string((long long) maxClusterSize);
	(*props)["maxIClusterSize"] = std::to_string((long long) maxClusterSize);
	(*props)["maxMixedClusterSize"] = std::to_string((long long) 2*maxClusterSize);
	(*props)["numHeClusters"] = std::to_string((long long) numClusters);
	(*props)["numVClusters"] = std::to_string((long long) numClusters);
	(*props)["numIClusters"] = std::to_string((long long) numClusters);
	(*props)["numMixedClusters"] = std::to_string(
			(long long) (numClusters * numClusters));

	// Write the cluster information to stdout
	BOOST_TEST_MESSAGE("Sizes of clusters in network:");
	for_each(reactants->begin(), reactants->end(), writeCluster);
	BOOST_TEST_MESSAGE(
			"Maximum He Cluster Size = " << (*props)["maxHeClusterSize"]);
	BOOST_TEST_MESSAGE(
			"Maximum V Cluster Size = " << (*props)["maxVClusterSize"]);
	BOOST_TEST_MESSAGE(
			"Maximum Interstitial Cluster Size = " << (*props)["maxIClusterSize"]);
	BOOST_TEST_MESSAGE("Number of He clusters = " << (*props)["numHeClusters"]);
	BOOST_TEST_MESSAGE("Number of V clusters = " << (*props)["numVClusters"]);
	BOOST_TEST_MESSAGE("Number of I clusters = " << (*props)["numIClusters"]);
	BOOST_TEST_MESSAGE(
			"Number of mixed clusters = " << (*props)["numMixedClusters"]);

	// Get the connectivity of the first reactant
	connectivityArray = reactants->at(0)->getConnectivity();
	// The connectivity array should be the same size as the reactants array
	BOOST_TEST_MESSAGE(
			"Connectivity Array Size = " << connectivityArray.size());
	BOOST_REQUIRE(connectivityArray.size() == reactants->size());

	// Since this is a helium cluster of size 1, it should interact with
	// every other helium up to size 9.
	for (int i = 0; i < numClusters - 1; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}
	// And not the last one
	BOOST_REQUIRE(connectivityArray.at(numClusters - 1) == 0);

	// Helium clusters can interact with vacancy clusters so long as the
	// maximum size of the helium atoms and vacancies is not greater than
	// the maximum mixed-species cluster size. This reactant should interact
	// with all of the vacancies.
	for (int i = numClusters; i < 2 * numClusters; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}

	// Helium can interact with a mixed-species cluster so long as the sum of
	// the number of helium atoms and the size of the mixed-species cluster
	// does not exceed the maximum mixed-species cluster size.
	for (int i = 3 * numClusters; i < reactants->size() - 1; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}
	// Single-species Helium can interact with all but the last.
	BOOST_REQUIRE(connectivityArray.at(reactants->size() - 1) == 0);

	return;
}
BOOST_AUTO_TEST_SUITE_END()

