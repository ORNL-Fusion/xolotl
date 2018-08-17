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
#include "SimpleReactionNetwork.h"
#include <PSIHeInterstitialCluster.h>
#include <PSIHeCluster.h>
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
 * This suite is responsible for testing the PSIHeInterstitialCluster.
 */
BOOST_AUTO_TEST_SUITE(PSIHeInterstitialCluster_testSuite)

BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	// Create a simple reaction network and create a HeInterstitial cluster
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork(0);
	PSIHeInterstitialCluster cluster(4, 2, *(network.get()), registry);

	// Get the composition back
	auto composition = cluster.getComposition();

	// Check the composition is the created one
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 4);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 0);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 2);

	// Check if it is a mixed cluster
	BOOST_REQUIRE_EQUAL(cluster.isMixed(), true);

	return;
}

/**
 * This operation checks the ability of the PSIHeInterstitialCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimplePSIReactionNetwork();

	// Check the reaction connectivity of the HeI cluster
	// with 5He and 3I
	// Get the connectivity array from the reactant
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 5;
	composition[toCompIdx(Species::V)] = 0;
	composition[toCompIdx(Species::I)] = 3;
	auto reactant = (PSICluster *) network->get(ReactantType::HeI, composition);
	// Check the type name
	BOOST_REQUIRE(ReactantType::HeI == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 5);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 3);

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 0, 0, 0, 0, 0,

			// D
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// T
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// V
			1, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// I
			1, 1, 1, 0, 0, 0, 0, 0, 0, 0,

			// Mixed
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// HeI
			0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the ability of the PSIHeInterstitialCluster to compute
 * the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {
	// Local Declarations
	auto network = getSimplePSIReactionNetwork();
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an HeI cluster with compostion 1,0,1.
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 1;
	composition[toCompIdx(Species::V)] = 0;
	composition[toCompIdx(Species::I)] = 1;
	auto cluster = (PSICluster *) network->get(ReactantType::HeI, composition);
	// Get one that it combines with (I)
	auto secondCluster = (PSICluster *) network->get(Species::I, 1);
	// Set the diffusion factor and migration energy to arbitrary values
	cluster->setDiffusionFactor(1.5E+10);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the preprocessor
	secondCluster->setDiffusionFactor(2.13E+10);
	secondCluster->setMigrationEnergy(0.013);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	// Check the flux
	BOOST_REQUIRE_CLOSE(-16982855380.0, flux, 0.1);

	return;
}

/**
 * This operation checks the PSIHeInterstitialCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 3.24895e+12, 0, 0, 0, 0, 0, -2.58738e+10, 0, 0,
			0, 0, 0, 0, 0, 0, 0 };
	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(2);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Get an HeI cluster with compostion 2,0,1.
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 1;
	composition[toCompIdx(Species::V)] = 0;
	composition[toCompIdx(Species::I)] = 1;
	auto cluster = (PSICluster *) network->get(ReactantType::HeI, composition);
	// Set the diffusion factor and migration energy to arbitrary values
	cluster->setDiffusionFactor(1.5E+10);
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
 * This operation checks the reaction radius for PSIHeInterstitialCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the HeI clsuter
	shared_ptr<PSIHeInterstitialCluster> cluster;

	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1372650265, 0.1778340462, 0.2062922619,
			0.2289478080, 0.2480795532 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<PSIHeInterstitialCluster>(
				new PSIHeInterstitialCluster(1, i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.000001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
