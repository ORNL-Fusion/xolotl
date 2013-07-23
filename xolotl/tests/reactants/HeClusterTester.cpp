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
#include <HeCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using std::shared_ptr;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HeCluster.
 */
BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/**
 * This operation checks the ability of the HeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	
	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	std::vector<shared_ptr<Reactant>> reactants = *network->reactants;
	std::map<std::string, std::string> props = *network->properties;
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the 6th He reactant (numHe=6)
	
	{
		// Get the index of the 6He reactant
		std::map<std::string, int> species;
		species["He"] = 6;
		int index = network->toClusterIndex(species);
		
		// Get the connectivity array from the reactant
		
		shared_ptr<PSICluster> reactant =
			std::dynamic_pointer_cast<PSICluster>(reactants.at(index));
		shared_ptr<std::vector<int>> reactionConnectivity =
			reactant->getConnectivity();
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// V
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// I
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// HeV
			1, 1, 1, 0, 0, 0, 0, 0, 0,
			1, 1, 0, 0, 0, 0, 0, 0,
			1, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0,
			0, 0,
			0,
			
			// HeI
			1, 1, 1, 0, 0, 0, 0, 0, 0,
			1, 1, 0, 0, 0, 0, 0, 0,
			1, 0, 0, 0, 0, 0, 0,
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

/**
 * This operation checks the HeCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
}

/**
 * This operation checks the reaction radius for HeCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	std::vector<std::shared_ptr<HeCluster>> clusters;
	std::shared_ptr<HeCluster> cluster;

	double expectedRadii[] = { 0.3, 0.3748419767, 0.4273418681, 0.4691369586,
			0.5044313198, 0.5352826768, 0.5628704922, 0.5879411911,
			0.6110006225, 0.6324092998 };

	for (int i = 1; i <= 10; i++) {
		cluster = std::shared_ptr<HeCluster>(new HeCluster(i));
		BOOST_CHECK_CLOSE(expectedRadii[i-1], cluster->getReactionRadius(), .000001);
	}
}


BOOST_AUTO_TEST_SUITE_END()

