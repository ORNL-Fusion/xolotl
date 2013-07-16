/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include "SimpleReactionNetwork.h"
#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include <map>
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the ReactionNetwork
 */
BOOST_AUTO_TEST_SUITE(ReactionNetwork_testSuite)

/**
 * This operation tests the copy constructor.
 */
BOOST_AUTO_TEST_CASE(checkCopying) {
	
	ReactionNetwork network;
	
	// Set some properties
	(*network.properties)["numMixedClusters"] = "4";
	
	// Add a reactant
	network.reactants->push_back(std::shared_ptr<Reactant>(new Reactant));
	network.reactants->at(0)->setConcentration(50.0);
	
	
	// Copy the network
	ReactionNetwork network2 = network;
	
	// Check that the ReactionNetwork fields are copied
	BOOST_REQUIRE_NE(network.properties.get(), network2.properties.get());
	BOOST_REQUIRE_EQUAL(network.properties->size(), network2.properties->size());
	
	BOOST_REQUIRE_NE(network.reactants.get(), network2.reactants.get());
	BOOST_REQUIRE_EQUAL(network.reactants->size(), network2.reactants->size());
	
	// Change the values of the copied ReactionNetwork
	// This shouldn't modify the Reactants contained inside the
	// first ReactionNetwork.
	
	(*network2.properties)["numMixedClusters"] = "5";
	network2.reactants->at(0)->setConcentration(52.0);
	
	BOOST_REQUIRE_NE(network.properties->at("numMixedClusters"),
		network2.properties->at("numMixedClusters"));
	BOOST_REQUIRE_CLOSE(network.reactants->at(0)->getConcentration(), 50.0, 1e-5);
	BOOST_REQUIRE_CLOSE(network2.reactants->at(0)->getConcentration(), 52.0, 1e-5);
}


BOOST_AUTO_TEST_CASE(clusterIndexConversions) {
	
	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	
	// Convert each index to a cluster map and back
	
	printf("i\tHe\tV\tI\n");
	printf("===\t===\t===\t===\n");
	
	int reactionLength = network->reactants->size();
	
	for (int index = 0; index < reactionLength; index++) {
		std::map<std::string, int> speciesMap = network->toClusterMap(index);
		printf("%d\t%d\t%d\t%d\n", index,
			speciesMap["He"], speciesMap["V"], speciesMap["I"]);
		int convertedIndex = network->toClusterIndex(speciesMap);
		
		BOOST_REQUIRE_EQUAL(convertedIndex, index);
	}
}

BOOST_AUTO_TEST_SUITE_END()
