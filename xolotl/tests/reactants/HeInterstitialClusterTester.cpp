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
#include <HeVCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;


/**
 * This suite is responsible for testing the HeVCluster.
 */
BOOST_AUTO_TEST_SUITE(HeInterstitialCluster_testSuite)


BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	HeVCluster cluster(4, 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("He"), 4);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("V"), 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("I"), 0);
}

/**
 * This operation checks the ability of the HeInterstitialCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	std::vector<shared_ptr<Reactant>> reactants = *network->reactants;
	std::map<std::string, std::string> props = *network->properties;
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the HeI cluster
	// with 5He and 3I
	
	{
		// Get the index of the 5He*3I reactant
		std::map<std::string, int> species;
		species["He"] = 5;
		species["I"] = 3;
		int index = network->toClusterIndex(species);
		
		// Get the connectivity array from the reactant
		
		shared_ptr<PSICluster> reactant =
			std::dynamic_pointer_cast<PSICluster>(reactants.at(index));
		shared_ptr<std::vector<int>> reactionConnectivity =
			reactant->getConnectivity();
		
		BOOST_REQUIRE_EQUAL(reactant->getClusterMap()["He"], 5);
		BOOST_REQUIRE_EQUAL(reactant->getClusterMap()["I"], 3);
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			
			// V
			1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
			
			// I
			// Only single-I clusters react with HeI
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			
			// HeV
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0,
			0, 0,
			0,
			
			// HeI
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0,
			0, 0,
			0
		};
		
		for (int i = 0; i < reactionConnectivity->size(); i++) {
			BOOST_REQUIRE_EQUAL(reactionConnectivity->at(i), connectivityExpected[i]);
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()
