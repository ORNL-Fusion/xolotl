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

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

/**
 * This suite is responsible for testing the HeCluster.
 */
BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/**
 * This operation checks the ability of the HeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	auto reactants = network->getAll();
	auto props = network->getProperties();
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the 6th He reactant (numHe=6)
	{
		// Get the connectivity array from the reactant
		auto reactant = dynamic_pointer_cast < PSICluster
				> (network->get("He", 6));
		auto reactionConnectivity = reactant->getConnectivity();
		
		// Check the connectivity for He, V, and I
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			
			// V
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// I
			1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			
			// HeV
			1, 1, 1, 0, 0, 1, 1, 1, 1,
			1, 1, 0, 0, 0, 1, 1, 1,
			1, 0, 0, 0, 0, 1, 1,
			0, 0, 0, 0, 0, 1,
			0, 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0,
			0, 0,
			0,
			
			// HeI
			1, 1, 1, 0, 0, 1, 1, 1, 1,
			1, 1, 0, 0, 0, 1, 1, 1,
			1, 0, 0, 0, 0, 1, 1,
			0, 0, 0, 0, 0, 1,
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

/**
 * This operation checks the HeCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

	// Get an He cluster with compostion 1,0,0.
	auto cluster = dynamic_pointer_cast<PSICluster>(network->get("He", 1));
	// Get one that it combines with (He2)
	auto secondCluster = dynamic_pointer_cast<PSICluster>(network->get("He", 2));
	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.950E+10);
	cluster->setMigrationEnergy(0.13);
	vector<double> energies = {numeric_limits<double>::infinity(),
			numeric_limits<double>::infinity(), numeric_limits<double>::infinity(), 8.27};
	cluster->setBindingEnergies(energies);
	cluster->setConcentration(0.5);

	// Set the diffusion factor, migration and binding energies based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(3.240E+010);
	secondCluster->setMigrationEnergy(0.2);
	energies = {0.864, numeric_limits<double>::infinity(),
			numeric_limits<double>::infinity(), 6.12};
	secondCluster->setBindingEnergies(energies);
	secondCluster->setConcentration(0.5);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(1000.0);
	BOOST_TEST_MESSAGE("HeClusterTester Message: \n" << "Total Flux is " << flux << "\n"
			  << "   -Production Flux: " << cluster->getProductionFlux(1000.0) << "\n"
			  << "   -Combination Flux: " << cluster->getCombinationFlux(1000.0) << "\n"
			  << "   -Dissociation Flux: " << cluster->getDissociationFlux(1000.0) << "\n");
	BOOST_REQUIRE_CLOSE(-43623893263., flux, 10);
}

/**
 * This operation checks the reaction radius for HeCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	vector<shared_ptr<HeCluster>> clusters;
	shared_ptr<HeCluster> cluster;

	double expectedRadii[] = { 0.3, 0.3237249066, 0.3403673722, 0.3536164159,
				0.3648047284, 0.3745846085, 0.3833299460, 0.3912773576,
				0.3985871973, 0.4053737480 };

	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<HeCluster>(new HeCluster(i));
		BOOST_REQUIRE_CLOSE(expectedRadii[i-1], cluster->getReactionRadius(), .000001);
	}
}


BOOST_AUTO_TEST_SUITE_END()

