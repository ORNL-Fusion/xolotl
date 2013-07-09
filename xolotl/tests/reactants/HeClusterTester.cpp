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
 * This suite is responsible for testing the HeCluster.
 */BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/**
 * This operation checks the ability of the HeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int numHe = 1, numV = 1;
	int maxClusterSize = 10, numClusters = maxClusterSize;
	vector<int> connectivityArray;
	shared_ptr<vector<shared_ptr<Reactant> > > reactants = network->reactants;
	shared_ptr<std::map<std::string, std::string>> props = network->properties;

	// Write the cluster information to stdout
	BOOST_TEST_MESSAGE("Sizes of clusters in network:");
	for_each(reactants->begin(), reactants->end(), writeCluster);
	BOOST_TEST_MESSAGE("Maximum He Cluster Size = " << (*props)["maxHeClusterSize"]);
	BOOST_TEST_MESSAGE("Maximum V Cluster Size = " << (*props)["maxVClusterSize"]);
	BOOST_TEST_MESSAGE("Maximum Interstitial Cluster Size = " << (*props)["maxIClusterSize"]);
	BOOST_TEST_MESSAGE("Number of He clusters = " << (*props)["numHeClusters"]);
	BOOST_TEST_MESSAGE("Number of V clusters = " << (*props)["numVClusters"]);
	BOOST_TEST_MESSAGE("Number of I clusters = " << (*props)["numIClusters"]);
	BOOST_TEST_MESSAGE("Number of mixed clusters = " << (*props)["numMixedClusters"]);

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
	// does not exceed the maximum mixed-species cluster size. Since
	// MixedSpeciesClusters start with a minimum size of two in this case, He
	// does not interact with the last two entries in the set.
	for (int i = 3 * numClusters; i < reactants->size() - 1; i++) {
		BOOST_TEST_MESSAGE("mixed index = " << i);
		BOOST_REQUIRE(connectivityArray.at(i) == 1);
	}
	// Single-species Helium can interact with all but the last.
	BOOST_REQUIRE(connectivityArray.at(reactants->size() - 1) == 0);

	return;
}
BOOST_AUTO_TEST_SUITE_END()

