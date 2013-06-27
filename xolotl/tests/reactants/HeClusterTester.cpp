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
 */
BOOST_AUTO_TEST_SUITE (HeCluster_testSuite)

/** This operation checks the ability of the HeCluster to describe its connectivity to other clusters. */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	// Local Declarations
	vector<int> connectivityArray;
	shared_ptr<ReactionNetwork> network = shared_ptr < ReactionNetwork
			> (new ReactionNetwork());
	vector<shared_ptr<Reactant>> reactants = *(network->reactants);
	std::map<std::string, std::string> props = *(network->properties);

	// Load the ReactionNetwork with some He clusters
	for (int i = 1; i < 11; i++) {
		// Create a He cluster with cluster size i
		shared_ptr<HeCluster> cluster = shared_ptr<HeCluster>(new HeCluster(i));
		// Register the network with the cluster
		cluster->setReactionNetwork(network);
		// Add it to the network
		reactants.push_back(cluster);
	}
	// Setup the properties map
	props["maxHeClusterSize"] = "10";
	props["maxVClusterSize"] = "1";
	props["maxIClusterSize"] = "1";
	props["numHeClusters"] = "10";
	props["numVClusters"] = "0";
	props["numIClusters"] = "0";
	props["numMixedClusters"] = "0";

	// Write the cluster information to stdout
	BOOST_TEST_MESSAGE("Sizes of clusters in network:");
	for_each(reactants.begin(),reactants.end(),writeCluster);
	BOOST_TEST_MESSAGE("Maximum He Cluster Size = " << props["maxHeClusterSize"]);
	BOOST_TEST_MESSAGE("Maximum V Cluster Size = " << props["maxVClusterSize"]);
	BOOST_TEST_MESSAGE("Maximum Interstitial Cluster Size = " << props["maxIClusterSize"]);
	BOOST_TEST_MESSAGE("Number of He clusters = " << props["numHeClusters"]);
	BOOST_TEST_MESSAGE("Number of V clusters = " << props["numVClusters"]);
	BOOST_TEST_MESSAGE("Number of I clusters = " << props["numIClusters"]);
	BOOST_TEST_MESSAGE("Number of mixed clusters = " << props["numMixedClusters"]);

	// Get the connectivity of the last reactant
	connectivityArray = reactants.at(9)->getConnectivity();
	// Since this is a helium cluster of size 10, it should interact with
	// everything up to size 9.
	BOOST_REQUIRE(connectivityArray.size() == 9);
	for (int i = 0; i < 9; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}
	// And not the last one
	BOOST_REQUIRE(0 == connectivityArray.at(9));

	return;
}

BOOST_AUTO_TEST_SUITE_END()

