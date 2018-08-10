#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <FeCluster.h>
#include "SimpleReactionNetwork.h"
#include <FeHeCluster.h>
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
 * This suite is responsible for testing the HeCluster.
 */
BOOST_AUTO_TEST_SUITE(HeCluster_testSuite)

/**
 * This operation checks the ability of the HeCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	shared_ptr<ReactionNetwork> network = getSimpleFeReactionNetwork();

	// Check the reaction connectivity of the 6th He reactant (numHe=6)
	// Get the connectivity array from the reactant
	auto reactant = (FeCluster *) network->get(Species::He, 6);

	// Check the type name
	BOOST_REQUIRE(ReactantType::He == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for He, V, and I
	int connectivityExpected[] = {
			// He
			1, 1, 1, 1, 1, 1, 1, 0,

			// V
			0, 0, 0, 0, 0, 0, 0, 0, 0,

			// I
			0,

			// HeV
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0,

			// temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	return;
}

/**
 * This operation checks the HeCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Local Declarations
	auto network = getSimpleFeReactionNetwork();

	// Get an He cluster with compostion 1,0,0.
	auto cluster = (FeCluster *) network->get(Species::He, 1);
	// Get one that it combines with (He2)
	auto secondCluster = (FeCluster *) network->get(Species::He, 2);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(1.0e+11);
	cluster->setMigrationEnergy(0.06);
	cluster->setConcentration(0.5);

	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem for the second cluster
	cluster->setDiffusionFactor(5.0e+10);
	cluster->setMigrationEnergy(0.06);
	secondCluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	network->reinitializeNetwork();
	network->computeRateConstants(0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);

	BOOST_REQUIRE_CLOSE(-236566273567.17, flux, 0.1);

	return;
}

/**
 * This operation checks the HeCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = { -1.50326e+12, -1.946346e+11, -1.990529e+11,
			-2.02802e+11, -2.059688e+11, -2.087369e+11, -2.112123e+11,
			3.26601367e+8, -1.382094e+11, -1.49712e+11, 0.0, 0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, -1.37276e+11, 1.47931e-12, -1.49695e+11, 1.60243e-12,
			0.0 };
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(2);

	// Get an He cluster with compostion 1,0,0.
	auto cluster = (FeCluster *) network->get(Species::He, 1);
	// Set the diffusion factor and migration energy based on the
	// values from the tungsten benchmark for this problem.
	cluster->setDiffusionFactor(1.0e+11);
	cluster->setMigrationEnergy(0.06);
	cluster->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
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
 * This operation checks the reaction radius for HeCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create a helium cluster
	shared_ptr<FeHeCluster> cluster;

	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.3, 0.321479647, 0.336547116, 0.3485423,
			0.3586717, 0.36752612, 0.3754438, 0.382639, 0.389257, 0.395401469 };

	// Check all the values
	for (int i = 1; i <= 10; i++) {
		cluster = shared_ptr<FeHeCluster>(
				new FeHeCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.0001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
