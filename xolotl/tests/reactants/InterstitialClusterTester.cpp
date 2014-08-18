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
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			1, 1, 1, 1, 1,
			1, 1, 1, 1,
			1, 1, 1,
			1, 1,
			1,
			
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

 	// Compute the rate constants that are needed for the flux
 	cluster->computeRateConstants(1000.0);
 	// The flux can pretty much be anything except "not a number" (nan).
 	double flux = cluster->getTotalFlux(1000.0);
 	BOOST_TEST_MESSAGE("InterstitialClusterTester Message: \n" << "Total Flux is " << flux << "\n"
 			  << "   -Production Flux: " << cluster->getProductionFlux(1000.0) << "\n"
 			  << "   -Combination Flux: " << cluster->getCombinationFlux(1000.0) << "\n"
 			  << "   -Dissociation Flux: " << cluster->getDissociationFlux(1000.0) << "\n"
 			  << "   -Emission Flux: " << cluster->getEmissionFlux(1000.0) << "\n");
 	BOOST_REQUIRE_CLOSE(-103424236193.0, flux, 0.1);
 }

 /**
  * This operation checks the HeCluster get*PartialDerivatives methods.
  */
 BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
 	// Local Declarations
 	// The vector of partial derivatives to compare with
 	double knownPartials[] = {-5.26951e+10, -5.54256e+10, 0.0, 0.0, -3.86349e+10,
 			-4.19101e+10, -2.90683e+11, 1.82094e+13, 5.1489e+12, 4.85134e+12, 0.0,
 			-3.86349e+10, -3.39657e+10, 0.0, 4.26504e+12};
 	// Get the simple reaction network
 	shared_ptr<ReactionNetwork> network = getSimpleReactionNetwork(3);

 	// Get an I cluster with compostion 0,0,1.
 	auto cluster = (PSICluster *) network->get("I", 1);
 	// Set the diffusion factor, migration and binding energies based on the
 	// values from the tungsten benchmark for this problem.
 	cluster->setDiffusionFactor(2.13E+10);
 	cluster->setMigrationEnergy(0.013);
 	vector<double> energies = {numeric_limits<double>::infinity(), numeric_limits<double>::infinity(),
 			numeric_limits<double>::infinity(), numeric_limits<double>::infinity()};
 	cluster->setBindingEnergies(energies);
	cluster->setTemperature(1000.0);
 	cluster->setConcentration(0.5);

 	// Compute the rate constants that are needed for the partial derivatives
 	cluster->computeRateConstants(1000.0);
 	// Get the vector of partial derivatives
 	auto partials = cluster->getPartialDerivatives(1000.0);

 	// Check the size of the partials
 	BOOST_REQUIRE_EQUAL(partials.size(), 15);

 	// Check all the values
 	for (int i = 0; i < partials.size(); i++) {
 		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
 	}
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
