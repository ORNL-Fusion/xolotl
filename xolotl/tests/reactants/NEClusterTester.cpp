#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <xolotlPerf/xolotlPerf.h>
#include <xolotlPerf/dummy/DummyHandlerRegistry.h>
#include <NECluster.h>
#include <NEClusterReactionNetwork.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>
#include "SimpleReactionNetwork.h"

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the NECluster.
 */
BOOST_AUTO_TEST_SUITE (NECluster_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkDiffusionCoefficient) {
	// Get the simple reaction network
	auto network = getSimpleNEReactionNetwork(0);
	// Set a fission rate for the diffusion to work
	network->setFissionRate(8.0e-9);
	// Create a cluster
	NECluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature
	cluster.addGridPoints(1);

	// Check 3 different temperature for each regime
	cluster.setMigrationEnergy(0.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(600.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0064001319, 0.00001);
	BOOST_REQUIRE_CLOSE(600.0, cluster.getTemperature(0), 0.0001);
	cluster.setTemperature(1500.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.20003754747, 0.00001);
	cluster.setTemperature(2000.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 18.114870067, 0.00001);

	return;
}

/**
 * This operation tests the default values returned by select flux routines.
 */
BOOST_AUTO_TEST_CASE(checkDefaultFluxes) {
	// Get the simple reaction network
	auto network = getSimpleNEReactionNetwork(0);
	// Create a cluster
	NECluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature
	cluster.addGridPoints(1);

	// Check the default values of the fluxes
	BOOST_REQUIRE_CLOSE(cluster.getProductionFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getCombinationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getDissociationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getTotalFlux(0), 0.0, 1e-5);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
