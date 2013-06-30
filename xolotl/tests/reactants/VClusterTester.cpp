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
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

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
 * This suite is responsible for testing the VCluster.
 */BOOST_AUTO_TEST_SUITE (VCluster_testSuite)

/**
 * This operation checks the ability of the VCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int maxClusterSize = 10, numClusters = maxClusterSize;
	vector<int> connectivityArray;
	shared_ptr < vector<shared_ptr<Reactant> > > reactants = network->reactants;
	shared_ptr<std::map<std::string, std::string>> props = network->properties;

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

	// Get the connectivity of the fifth vacancy (index 14)
	connectivityArray = reactants->at(numClusters + 4)->getConnectivity();
	// The connectivity array should be the same size as the reactants array
	BOOST_TEST_MESSAGE(
			"Connectivity Array Size = " << connectivityArray.size());
	BOOST_REQUIRE(connectivityArray.size() == reactants->size());

	// Since this is a vacancy cluster of size 5, it should not interact with
	// vacancies bigger than maxClusterSize - 5 (which is conveniently 5 in
	// this case). So check the small vacancies first...
	for (int i = numClusters; i < numClusters + maxClusterSize - 5; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}
	// ...and the big vacancies second.
	for (int i = numClusters + maxClusterSize - 5;
			i < numClusters + maxClusterSize; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 0);
	}

	// Vacancies can interact with everything else, within size limits.

	// Check single-species He. Since this is a vacancy of size 5, it can
	// interact with all of the single-species helium (all size 10 or less)
	// in this network.
	for (int i = 0; i < numClusters; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}

	// Check single-species interstitials
	for (int i = 2 * numClusters; i < 3 * numClusters; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}

	// Check mixed species
	for (int i = 3 * numClusters; i < reactants->size() - 5; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}

	return;
}
BOOST_AUTO_TEST_SUITE_END()

