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
BOOST_AUTO_TEST_SUITE(HeVCluster_testSuite)


BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	HeVCluster cluster(4, 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("He"), 4);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("V"), 5);
	BOOST_REQUIRE_EQUAL(cluster.getSpeciesSize("I"), 0);
}

/**
 * This operation checks the ability of the HeVCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	shared_ptr<ReactionNetwork> network = testUtils::getSimpleReactionNetwork();
	auto reactants = network->getAll();
	auto props = network->getProperties();
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the HeV cluster
	// with 3He and 2V
	
	{
		// Get the index of the 3He*2V reactant
		std::map<std::string, int> species;
		species["He"] = 3;
		species["V"] = 2;
		int index = network->toClusterIndex(species);
		
		// Get the connectivity array from the reactant
		
		shared_ptr<PSICluster> reactant =
			std::dynamic_pointer_cast<PSICluster>(reactants->at(index));
		shared_ptr<std::vector<int>> reactionConnectivity =
			reactant->getConnectivity();
		
		BOOST_REQUIRE_EQUAL(reactant->getClusterMap()["He"], 3);
		BOOST_REQUIRE_EQUAL(reactant->getClusterMap()["V"], 2);
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			
			// V
			// Only single-V clusters react with HeV
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			
			// I
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

BOOST_AUTO_TEST_CASE(checkGetFlux) {
	/*shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	int maxClusterSize = 10;
	int numClusters = maxClusterSize;

	shared_ptr<vector<shared_ptr<Reactant>>> reactants = network->reactants;
	shared_ptr<map<string, string>> props = network->properties;

	for (int i = 0; i < reactants->size(); i++) {
		reactants->at(i)->setConcentration(2.0*i);
		(dynamic_pointer_cast<PSICluster>(reactants->at(i)))->setMigrationEnergy(1.0);
		(dynamic_pointer_cast<PSICluster>(reactants->at(i)))->setDiffusionFactor(1.0);
	}
	// Get the connectivity of the 20th HeV cluster (index 59).
	shared_ptr<Reactant> reactant = reactants->at(3 * numClusters + 20 - 1);
	shared_ptr<HeVCluster> cluster = dynamic_pointer_cast<HeVCluster>(reactant);

	double flux = cluster->getProductionFlux(1.0);

	std::cout << "Flux is " << flux << "\n";*/
}

/**
 * This operation checks the reaction radius for HeVCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	std::vector<std::shared_ptr<HeVCluster>> clusters;
	std::shared_ptr<HeVCluster> cluster;
	double expectedRadii[] = { 0.4330127019, 0.5609906819, 0.6507642333,
			0.7222328328, 0.7825853415 };

	for (int i = 1; i <= 5; i++) {
		cluster = std::shared_ptr<HeVCluster>(new HeVCluster(1, i));
		BOOST_CHECK_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				.000001);
	}
}

BOOST_AUTO_TEST_SUITE_END()
