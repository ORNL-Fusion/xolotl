#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <FeCluster.h>
#include "SimpleReactionNetwork.h"
#include <FeHeCluster.h>
#include <FeHeVCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <math.h>
#include <limits>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the FeHeVCluster.
 */
BOOST_AUTO_TEST_SUITE(FeHeVCluster_testSuite)

BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	// Create a simple reaction network and create a HeV cluster
	shared_ptr<ReactionNetwork> network = getSimpleFeReactionNetwork(0);
	FeHeVCluster cluster(4, 5, *(network.get()), registry);

	// Get the composition back
	auto composition = cluster.getComposition();

	// Check the composition is the created one
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 4);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 5);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::I)], 0);

	// Check if it is a mixed cluster
	BOOST_REQUIRE_EQUAL(cluster.isMixed(), true);

	return;
}

/**
 * This operation checks the ability of the FeHeVCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimpleFeReactionNetwork();

	// Check the reaction connectivity of the FeHeV cluster
	// with 3He and 2V
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 2;
	composition[toCompIdx(Species::V)] = 2;
	composition[toCompIdx(Species::I)] = 0;
	auto reactant = (FeCluster *) network->get(ReactantType::HeV, composition);

	// Check the type name
	BOOST_REQUIRE(ReactantType::HeV == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the composition
	composition = reactant->getComposition();
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::He)], 2);
	BOOST_REQUIRE_EQUAL(composition[toCompIdx(Species::V)], 2);

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 0, 0, 0, 0, 0,

			// V
			1, 1, 0, 0, 0, 0, 0, 0, 0,

			// I
			1,

			// HeV
			0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the ability of the FeHeVCluster to compute the total flux.
 */
BOOST_AUTO_TEST_CASE(checkTotalFlux) {
	// Local Declarations
	auto network = getSimpleFeReactionNetwork();

	// Get an HeV cluster with compostion 2,1,0.
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 2;
	composition[toCompIdx(Species::V)] = 1;
	composition[toCompIdx(Species::I)] = 0;
	auto cluster = (FeCluster *) network->get(ReactantType::HeV, composition);
	// Get one that it combines with (He)
	auto secondCluster = (FeCluster *) network->get(Species::He, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the preprocessor.
	cluster->setDiffusionFactor(0.0);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	secondCluster->setDiffusionFactor(1.0e+11);
	secondCluster->setMigrationEnergy(0.06);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	network->reinitializeNetwork();
	network->computeRateConstants(0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-69104736618, flux, 0.1);

	return;
}

/**
 * This operation checks the FeHeVCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { 0.0, 7.2852e-13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			-1.77577, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.77577, 0.0,
			-3.05711e-07, 0.0, 4.93257e-08, 0.0 };
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(2);

	// Get an HeV cluster with compostion 2,1,0.
	IReactant::Composition composition;
	composition[toCompIdx(Species::He)] = 2;
	composition[toCompIdx(Species::V)] = 1;
	composition[toCompIdx(Species::I)] = 0;
	auto cluster = (FeCluster *) network->get(ReactantType::HeV, composition);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(0.0);
	cluster->setMigrationEnergy(numeric_limits<double>::infinity());
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
 * This operation checks the reaction radius for FeHeVCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the HeV cluster
	shared_ptr<FeHeVCluster> cluster;

	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.1413109, 0.17804059, 0.2038056, 0.224317088,
			0.24163825 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<FeHeVCluster>(
				new FeHeVCluster(1, i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.0001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
