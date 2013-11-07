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
#include <math.h>
#include <limits>

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
		std::vector<int> reactionConnectivity =
			reactant->getConnectivity();
		
		BOOST_REQUIRE_EQUAL(reactant->getComposition()["He"], 3);
		BOOST_REQUIRE_EQUAL(reactant->getComposition()["V"], 2);
		
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
		
		for (int i = 0; i < reactionConnectivity.size(); i++) {
			BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
		}
	}
}

BOOST_AUTO_TEST_CASE(checkTotalFlux) {

	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

	// Get an HeV cluster with sizes 1,1,0.
	std::vector<int> composition = {1,1,0};
	auto cluster = std::dynamic_pointer_cast<PSICluster>(network->getCompound("HeV",composition));
	// Get one that it combines with
	composition.at(1) = 2;
	auto secondCluster = std::dynamic_pointer_cast<PSICluster>(network->getCompound("HeV",composition));

	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(0.0);
	cluster->setMigrationEnergy(std::numeric_limits<double>::infinity());
	std::vector<double> energies = {5.09,5.09,std::numeric_limits<double>::infinity()};
	cluster->setBindingEnergies(energies);

	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(0.0);
	secondCluster->setMigrationEnergy(std::numeric_limits<double>::infinity());
	energies = {5.39,0.725,std::numeric_limits<double>::infinity()};
	secondCluster->setBindingEnergies(energies);

	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(1000.0);
	std::cout.precision(15);
	std::cout << "HeVClusterTester Message: " << " Flux is " << flux << "\n";
	BOOST_REQUIRE(!std::isnan(flux));
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
