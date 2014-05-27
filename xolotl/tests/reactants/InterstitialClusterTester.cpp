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
#include <InterstitialCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry = std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the InterstitialCluster.
 */BOOST_AUTO_TEST_SUITE(InterstitialCluster_testSuite)

/**
 * This operation checks the ability of the InterstitialCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	
	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();
	auto props = network->getProperties();
	
	// Prevent dissociation from being added to the connectivity array
	props["dissociationsEnabled"] = "false";
	
	// Check the reaction connectivity of the 4th interstitial cluster (4I)
	
	{
		// Get the connectivity array from the reactant
		auto reactant = (PSICluster *) network->get("I", 4);
		// Check the type name
		BOOST_REQUIRE_EQUAL("I",reactant->getType());
		auto reactionConnectivity = reactant->getConnectivity();
		
		// Check the connectivity for He, V, and I
		
		int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
			
			// V
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			
			// I
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			
			// HeV
			1, 1, 1, 1, 1, 0, 0, 0, 0,
			1, 1, 1, 1, 0, 0, 0, 0,
			1, 1, 1, 0, 0, 0, 0,
			1, 1, 0, 0, 0, 0,
			1, 1, 1, 1, 1,
			1, 1, 1, 1,
			1, 1, 1,
			1, 1,
			1,
			
			// HeI
			// Only a single-Interstitial cluster can react with a HeI cluster
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 1, 1,
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
  * This operation checks the InterstitialCluster get*Flux methods.
  */
 BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
 	// Local Declarations
 	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork();

 	// Get an I cluster with compostion 0,0,1.
 	auto cluster = (PSICluster *) network->get("I", 1);
 	// Get one that it combines with (I2)
 	auto secondCluster = (PSICluster *) network->get("I", 2);
 	// Set the diffusion factor, migration and binding energies based on the
 	// values from the tungsten benchmark for this problem.
 	cluster->setDiffusionFactor(2.13E+10);
 	cluster->setMigrationEnergy(0.013);
 	vector<double> energies = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(),
 			numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
 	cluster->setBindingEnergies(energies);
	cluster->setTemperature(1000.0);
 	cluster->setConcentration(0.5);

 	// Set the diffusion factor, migration and binding energies based on the
 	// values from the tungsten benchmark for this problem for the second cluster
 	secondCluster->setDiffusionFactor(1.065E+10);
 	secondCluster->setMigrationEnergy(0.013);
	secondCluster->setTemperature(1000.0);
 	energies = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(),
 			2.12, numeric_limits<double>::infinity()};
 	secondCluster->setBindingEnergies(energies);
 	secondCluster->setConcentration(0.5);
 	// The flux can pretty much be anything except "not a number" (nan).
 	double flux = cluster->getTotalFlux(1000.0);
 	BOOST_TEST_MESSAGE("InterstitialClusterTester Message: \n" << "Total Flux is " << flux << "\n"
 			  << "   -Production Flux: " << cluster->getProductionFlux(1000.0) << "\n"
 			  << "   -Combination Flux: " << cluster->getCombinationFlux(1000.0) << "\n"
 			  << "   -Dissociation Flux: " << cluster->getDissociationFlux(1000.0) << "\n");
 	BOOST_REQUIRE_CLOSE(-67088824870., flux, 10);
 }

/**
 * This operation checks the reaction radius for InterstitialCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {

	vector<shared_ptr<InterstitialCluster>> clusters;
	shared_ptr<InterstitialCluster> cluster;
	double expectedRadii[] = { 0.1578547805, 0.1984238001, 0.2268820159,
			0.2495375620, 0.2686693072, 0.2853926671, 0.3003469838,
			0.3139368664, 0.3264365165, 0.3380413550 };

	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<InterstitialCluster>(
				new InterstitialCluster(i, registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				.000001);
	}
}

BOOST_AUTO_TEST_SUITE_END()
