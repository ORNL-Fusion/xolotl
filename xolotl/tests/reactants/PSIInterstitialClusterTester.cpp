/*
 * PSIClusterTester.cpp
 *
 *  Created on: May 6, 2013
 *      Author: Jay Jay Billings
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <PSICluster.h>
#include <PSIInterstitialCluster.h>
#include "SimpleReactionNetwork.h"
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the PSIInterstitialCluster.
 */
BOOST_AUTO_TEST_SUITE(PSIInterstitialCluster_testSuite)

/**
 * This operation checks the ability of the PSIInterstitialCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork();

	// Check the reaction connectivity of the 4th interstitial cluster (4I)
	// Get the connectivity array from the reactant
	auto reactant = (PSICluster *) network->get(Species::I, 4);

	// Check the type name
	BOOST_REQUIRE(ReactantType::I == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 0, 0, 0, 0,

			// D
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// T
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// V
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

			// I
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

			// Mixed
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
			0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

			// HeI
			1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
			0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the PSIInterstitialCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	auto network = getSimplePSIReactionNetwork();
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an I cluster with compostion 0,0,1.
	auto cluster = (PSICluster *) network->get(Species::I, 1);
	// Get one that it combines with (I2)
	auto secondCluster = (PSICluster *) network->get(Species::I, 2);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.13E+10);
	cluster->setMigrationEnergy(0.013);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(1.065E+10);
	secondCluster->setMigrationEnergy(0.013);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(9021773486621.2, flux, 0.1);

	return;
}

/**
 * This operation checks the PSIInterstitialCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] =
			{ -5.26951e+10, 0, 0, 0, 0, 0, -3.39657e+10, -3.86349e+10,
					-2.90683e+11, 1.82504e+13, -3.39657e+10, 0, 0, 0, 0, 0 };
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(2);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an I cluster with compostion 0,0,1.
	auto cluster = (PSICluster *) network->get(Species::I, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.13E+10);
	cluster->setMigrationEnergy(0.013);
	cluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the partials
	network->setTemperature(1000.0, 0);
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 16U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	return;
}

/**
 * This operation checks the reaction radius for PSIInterstitialCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the interstitial cluster
	shared_ptr<PSIInterstitialCluster> cluster;
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1578547805, 0.1984238001, 0.2268820159,
			0.2495375620, 0.2686693072, 0.2853926671, 0.3003469838,
			0.3139368664, 0.3264365165, 0.3380413550 };

	// Check all the values
	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<PSIInterstitialCluster>(
				new PSIInterstitialCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.000001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
