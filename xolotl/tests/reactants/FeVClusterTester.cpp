#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <FeCluster.h>
#include <FeVCluster.h>
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
 * This suite is responsible for testing the FeVCluster.
 */BOOST_AUTO_TEST_SUITE(FeVCluster_testSuite)

/**
 * This operation checks the ability of the FeVCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimpleFeReactionNetwork();

	// Get the connectivity array from the reactant for a vacancy cluster of size 2.
	auto reactant = (FeCluster *) network->get(Species::V, 2);
	// Check the type name
	BOOST_REQUIRE(ReactantType::V == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 0, 0, 0,

			// V
			1, 1, 1, 0, 0, 0, 0, 0, 0,

			// I
			1,

			// HeV
			0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the FeVCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	auto network = getSimpleFeReactionNetwork();

	// Get an V cluster with compostion 0,1,0.
	auto cluster = (FeCluster *) network->get(Species::V, 1);
	// Get one that it combines with (V2)
	auto secondCluster = (FeCluster *) network->get(Species::V, 2);
	// Set the diffusion factor, and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(1.0e+11);
	cluster->setMigrationEnergy(0.67);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(5.0e+10);
	secondCluster->setMigrationEnergy(0.62);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	network->reinitializeNetwork();
	network->computeRateConstants(0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(2013449798, flux, 0.1);

	return;
}

/**
 * This operation checks the FeVCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { -1.16486e+08, -1.22155e+08, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, -5.96809e+08, 4.4006e+09, -1.75765e+07, -8.03182e+07,
			-9.51941e+07, -1.02278e+08, -1.07081e+08, -1.10915e+08, 666012,
			-7.45991e+07, -7.38128e+07, -7.45991e+07, 8730.85, 2.07215, 0.0 };
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(2);

	// Get an V cluster with compostion 0,1,0.
	auto cluster = (FeCluster *) network->get(Species::V, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(1.0e+11);
	cluster->setMigrationEnergy(0.67);
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
 * This operation checks the reaction radius for FeVCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the vacancy clsuter
	shared_ptr<FeVCluster> cluster;

	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1413109, 0.17804059, 0.2038056, 0.224317088,
			0.241638258 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<FeVCluster>(
				new FeVCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.0001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()

