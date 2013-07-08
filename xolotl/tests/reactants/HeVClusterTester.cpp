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
 * This suite is responsible for testing the HeVCluster.
 */BOOST_AUTO_TEST_SUITE(HeVCluster_testSuite)

/**
 * This operation checks the ability of the HeVCluster to describe
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
	BOOST_TEST_MESSAGE("Number of He clusters = " << (*props)["numHeClusters"]);
	BOOST_TEST_MESSAGE("Number of V clusters = " << (*props)["numVClusters"]);
	BOOST_TEST_MESSAGE(
			"Number of mixed clusters = " << (*props)["numMixedClusters"]);

	// Get the connectivity of the fifth mixed-species cluster (index 34). It
	// has three helium atoms and three vacancies.
	connectivityArray = reactants->at(3 * numClusters + 4)->getConnectivity();
	// The connectivity array should be the same size as the reactants array
	BOOST_TEST_MESSAGE(
			"Connectivity Array Size = " << connectivityArray.size());
	BOOST_REQUIRE(connectivityArray.size() == reactants->size());

	// This cluster should only interact with a few specific clusters:
	// >single He - helium dissociation
	// >single V - vacancy dissociation
	// >single I - interstitial absorption
	// >[(A-1)*He](B*V) - helium dissociation
	// >(A*He)*[(B-1)*V] - vacancy dissociation
	// >(A*He)*[(B+1)*V] - interstitial absorption
	
	for (int a : connectivityArray)
		printf("%d\n", a);
	
	BOOST_REQUIRE(connectivityArray.at(0) == 1);
	BOOST_REQUIRE(connectivityArray.at(numClusters-1) == 1);
	BOOST_REQUIRE(connectivityArray.at(2*numClusters-1) == 1);
	BOOST_REQUIRE(connectivityArray.at(3*numClusters+2) == 1);
	BOOST_REQUIRE(connectivityArray.at(3*numClusters+3) == 1);
	BOOST_REQUIRE(connectivityArray.at(3*numClusters+6) == 1);

	return;
}
BOOST_AUTO_TEST_SUITE_END()

