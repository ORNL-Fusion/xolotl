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
#include <InterstitialCluster.h>
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
 * This suite is responsible for testing the InterstitialCluster.
 */BOOST_AUTO_TEST_SUITE(InterstitialCluster_testSuite)

/**
 * This operation checks the ability of the InterstitialCluster to describe
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

	// Get the connectivity of the fifth interstitial (index 24)
	connectivityArray = reactants->at(2 * numClusters + 4)->getConnectivity();
	// The connectivity array should be the same size as the reactants array
	BOOST_TEST_MESSAGE(
			"Connectivity Array Size = " << connectivityArray.size());
	BOOST_REQUIRE(connectivityArray.size() == reactants->size());

	// Since this is an interstitial cluster of size 5, it should not interact with
	// interstitials bigger than maxClusterSize - 5 (which is conveniently 5 in
	// this case). So check the small interstitials first...
	for (int i = 2 * numClusters; i < 2 * numClusters + maxClusterSize - 5;
			i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}
	// ...and the big interstitials second.
	for (int i = 2 * numClusters + maxClusterSize - 5;
			i < 2 * numClusters + maxClusterSize; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 0);
	}

	// Interstitials can interact with other interstitials, vacancies and
	// mixed-species clusters, but not helium. They cannot cluster with other
	// interstitials that are so large that the combination of the two would
	// produce an interstitial above the maximum size.

	// Check single-species He. Interstititals should never interact with He.
	for (int i = 0; i < numClusters; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 0);
	}

	// Check single-species V. Interstitials should always interact with V.
	for (int i = numClusters; i < 2 * numClusters; i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}

	// Check mixed species. Interstitials should always interact with
	// mixed-species clusters because they only reduce the size, never
	// increase it.
	for (int i = 3 * numClusters; i < reactants->size(); i++) {
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}

	return;
}

/**
 * This operation checks the reaction radius for InterstitialCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	std::vector<std::shared_ptr<InterstitialCluster>> clusters;
	std::shared_ptr<InterstitialCluster> cluster;
	double expectedRadii[] = { 0.4979646072, 0.6259425872, 0.7157161386,
			0.7871847381, 0.8475372467, 0.9002923252, 0.9474668259,
			0.9903371181, 1.0297681911, 1.0663765142 };

	for (int i = 1; i <= 10; i++) {
		cluster = std::shared_ptr<InterstitialCluster>(
				new InterstitialCluster(i));
		BOOST_CHECK_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				.000001);
	}
}

BOOST_AUTO_TEST_SUITE_END()

