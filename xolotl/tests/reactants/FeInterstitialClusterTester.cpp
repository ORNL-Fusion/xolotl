#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <FeCluster.h>
#include <FeInterstitialCluster.h>
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
 * This suite is responsible for testing the FeInterstitialCluster.
 */
BOOST_AUTO_TEST_SUITE(FeInterstitialCluster_testSuite)

/**
 * This operation checks the ability of the FeInterstitialCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimpleFeReactionNetwork();

	// Check the reaction connectivity of the first interstitial cluster (I)
	// Get the connectivity array from the reactant
	auto reactant = (FeCluster *) network->get(Species::I, 1);

	// Check the type name
	BOOST_REQUIRE(ReactantType::I == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 0, 0, 0,

			// V
			1, 1, 1, 1, 1, 1, 1, 1, 1,

			// I
			1,

			// HeV
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the FeInterstitialCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	auto network = getSimpleFeReactionNetwork();

	// Get an I cluster with compostion 0,0,1.
	auto cluster = (FeCluster *) network->get(Species::I, 1);
	// Get one that it combines with (V)
	auto secondCluster = (FeCluster *) network->get(Species::V, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(1.0e+11);
	cluster->setMigrationEnergy(0.34);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(1.0e+11);
	secondCluster->setMigrationEnergy(0.67);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	network->reinitializeNetwork();
	network->computeRateConstants(0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-1754900681, flux, 0.1);

	return;
}

/**
 * This operation checks the FeInterstitialCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 1.04524e-09, 0.0133446, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, -3.25273e+10, -3.67546e+10, -3.97199e+10, -4.20806e+10,
			-4.40741e+10, -4.58167e+10, -4.73749e+10, -4.8791e+10, -5.00934e+10,
			-8.738435e+06, -3.25273e+10, -3.25273e+10, -3.67546e+10,
			-3.67546e+10, 0.0 };
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(2);

	// Get an I cluster with compostion 0,0,1.
	auto cluster = (FeCluster *) network->get(Species::I, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(2.13E+10);
	cluster->setMigrationEnergy(0.013);
	cluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the partials
	network->setTemperature(1000.0, 0);
	network->reinitializeNetwork();
	network->computeRateConstants(0);
	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 23U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	return;
}

/**
 * This operation checks the reaction radius for FeInterstitialCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the interstitial cluster
	shared_ptr<FeInterstitialCluster> cluster;
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1413109 };

	// Check all the values
	for (int i = 1; i <= 1; i++) {
		cluster = shared_ptr<FeInterstitialCluster>(
				new FeInterstitialCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.0001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
