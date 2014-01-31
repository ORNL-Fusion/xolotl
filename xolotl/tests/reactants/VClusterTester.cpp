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
#include <VCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

/**
 * This suite is responsible for testing the VCluster.
 */BOOST_AUTO_TEST_SUITE(VCluster_testSuite)

/**
 * This operation checks the ability of the VCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {

	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	auto props = network->getProperties();

	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";

	// Get the connectivity array from the reactant for a vacancy cluster of size 2.
	auto reactant = dynamic_pointer_cast < PSICluster
			> (network->get("V", 2));
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
			// V
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

			// HeV
			// The VCluster type only reacts with HeV for
			// single-V clusters.
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0,
			0, 0,
			0,

			// HeI
			1, 1, 1, 1, 1, 1, 1, 0, 0,
			1, 1, 1, 1, 1, 1, 0, 0,
			1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1,
			1, 1, 1, 1,
			1, 1, 1,
			1, 1,
			1
	};

	for (int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_TEST_MESSAGE("Connectivity [" << i << "] = " << reactionConnectivity[i]);
	}
	for (int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i],
				connectivityExpected[i]);
	}
}

 /**
  * This operation checks the VCluster get*Flux methods.
  */
 BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
 	// Local Declarations
 	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

 	// Get an V cluster with compostion 0,1,0.
 	auto cluster = dynamic_pointer_cast<PSICluster>(network->get("V", 1));
 	// Get one that it combines with (V2)
 	auto secondCluster = dynamic_pointer_cast<PSICluster>(network->get("V", 2));
 	// Set the diffusion factor, migration and binding energies based on the
 	// values from the tungsten benchmark for this problem.
 	cluster->setDiffusionFactor(2.41E+11);
 	cluster->setMigrationEnergy(1.66);
 	vector<double> energies = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(),
 			numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
 	cluster->setBindingEnergies(energies);
 	cluster->setConcentration(0.5);

 	// Set the diffusion factor, migration and binding energies based on the
 	// values from the tungsten benchmark for this problem for the second cluster
 	secondCluster->setDiffusionFactor(0.);
 	secondCluster->setMigrationEnergy(numeric_limits<double>::infinity());
 	energies = {numeric_limits<double>::infinity(), 0.432,
 			numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
 	secondCluster->setBindingEnergies(energies);
 	secondCluster->setConcentration(0.5);
 	// The flux can pretty much be anything except "not a number" (nan).
 	double flux = cluster->getTotalFlux(1000.0);
 	BOOST_TEST_MESSAGE("InterstitialClusterTester Message: \n" << "Total Flux is " << flux << "\n"
 			  << "   -Production Flux: " << cluster->getProductionFlux(1000.0) << "\n"
 			  << "   -Combination Flux: " << cluster->getCombinationFlux(1000.0) << "\n"
 			  << "   -Dissociation Flux: " << cluster->getDissociationFlux(1000.0) << "\n");
 	BOOST_REQUIRE_CLOSE(-2684., flux, 10);
 }

/**
 * This operation checks the reaction radius for VCluster.
 */
 BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	BOOST_TEST_MESSAGE("VClustertTester Message: BOOST_AUTO_TEST_CASE(checkReactionRadius): \n"
			<< "getReactionRadius needs to be fixed");
}
BOOST_AUTO_TEST_SUITE_END()

