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
#include <math.h>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

/**
 * This suite is responsible for testing the Reactant.
 */
BOOST_AUTO_TEST_SUITE(Reactant_testSuite)

/**
 * This operation tests the copy constructor.
 */
BOOST_AUTO_TEST_CASE(checkCopying) {

	// Create a reference Reactant
	shared_ptr<Reactant> reactant(new Reactant);
	reactant->setConcentration(10.0);

	// Copy the Reactant
	shared_ptr<Reactant> reactantCopy(new Reactant(*reactant));

	// Check that the pointers are different
	BOOST_REQUIRE_NE(reactant.get(), reactantCopy.get());

	reactantCopy->increaseConcentration(5.0);

	// The values should now be different,
	// so check them against the known values
	BOOST_REQUIRE_CLOSE(reactant->getConcentration(), 10.0, 1e-7);
	BOOST_REQUIRE_CLOSE(reactantCopy->getConcentration(), 15.0, 1e-7);

	// Try cloning the Reactant
	auto reactantClone = reactant->clone();
	BOOST_REQUIRE_CLOSE(10.0, reactantClone->getConcentration(), 1e-7);
}

BOOST_AUTO_TEST_CASE(checkManipulateConcentration) {

	// Create a Reactant
	shared_ptr<Reactant> reactant(new Reactant);
	reactant->setConcentration(1.0);

	// Make sure it was set correctly
	BOOST_REQUIRE_EQUAL(1.0, reactant->getConcentration());

	// Increase it
	reactant->increaseConcentration(3.3);

	// Make sure its correct
	BOOST_REQUIRE_EQUAL(4.3, reactant->getConcentration());

	// Decrease it
	reactant->decreaseConcentration(1.3);

	// Make sure its correct
	BOOST_REQUIRE_EQUAL(3.0, reactant->getConcentration());

	// Zero it
	reactant->zero();

	// Check it was zeroed
	BOOST_REQUIRE_EQUAL(0.0, reactant->getConcentration());

	// Make sure the base class getTotalFlux returns 0 for now
	BOOST_REQUIRE_EQUAL(0.0, reactant->getTotalFlux(0.0));

}

BOOST_AUTO_TEST_CASE(toClusterMap) {

 	BOOST_TEST_MESSAGE("ReactantTester Message: toClusterMap returns an empty map");

//	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
//
//	map<string, int> cluster;
//
//	// Test a couple of the HeClusters
//
//	cluster = network->toClusterMap(0);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 1);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
//
//	cluster = network->toClusterMap(9);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 10);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
//
//	// Test VClusters
//
//	cluster = network->toClusterMap(10);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 1);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
//
//	cluster = network->toClusterMap(19);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 10);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
//
//	// Test IClusters
//
//	cluster = network->toClusterMap(20);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 1);
//
//	cluster = network->toClusterMap(29);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 10);
//
//	// Test HeVClusters
//
//	cluster = network->toClusterMap(40);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 2);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 2);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
//
//	cluster = network->toClusterMap(60);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 1);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 5);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 0);
//
//	// Test HeInterstitialClusters
//
//	cluster = network->toClusterMap(80);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 6);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 1);
//
//	cluster = network->toClusterMap(100);
//	BOOST_REQUIRE_EQUAL(cluster["He"], 2);
//	BOOST_REQUIRE_EQUAL(cluster["V"], 0);
//	BOOST_REQUIRE_EQUAL(cluster["I"], 4);
}

BOOST_AUTO_TEST_CASE(checkIsConnected) {
	// Create a reaction network containing only clusters with maximum size 2
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork(2);

	// Check the connectivity matrix (8 * 8)
	int connectivityExpected[8][8] = {
			1, 1, 1, 0, 1, 0, 1, 1,
			1, 1, 0, 0, 0, 0, 0, 0,
			1, 0, 1, 1, 1, 1, 1, 0,
			0, 0, 1, 1, 1, 1, 0, 0,
			1, 0, 1, 1, 1, 1, 0, 1,
			0, 0, 1, 1, 1, 1, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 1,
			0, 0, 0, 0, 0, 0, 1, 0
	};

	// Initialyze i and j to access the connectivityExpected matrix
	int i = 0;
	int j = 0;

	// Get the connectivity matrix from the network
	auto reactants = network->getAll();
	for (auto reactantIt = reactants->begin();
			reactantIt != reactants->end(); reactantIt++) {
		shared_ptr<PSICluster> cluster = dynamic_pointer_cast<
				PSICluster>(*reactantIt);
		vector<int> reactionConnectivity = cluster->getConnectivity();

		for (auto connIt = reactionConnectivity.begin(); connIt != reactionConnectivity.end(); connIt++) {
			// Compare values between what is in the network and expected
			BOOST_REQUIRE_EQUAL(*connIt, connectivityExpected[i][j]);
			j++;
		}
		i++;
		// j has to be reset to 0 when i is incremented
		j = 0;
	}

}


BOOST_AUTO_TEST_SUITE_END()

