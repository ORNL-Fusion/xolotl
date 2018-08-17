#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>
#include <FeCluster.h>
#include <FeClusterReactionNetwork.h>
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
 * This suite is responsible for testing the FeCluster.
 */
BOOST_AUTO_TEST_SUITE (FeCluster_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkDiffusionCoefficient) {
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(0);
	// Create a cluster
	FeCluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature and diffusion coef
	cluster.addGridPoints(1);

	// Check E_m = 0.0
	cluster.setMigrationEnergy(0.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), exp(0.0), 0.00001);
	BOOST_REQUIRE_CLOSE(1.0, cluster.getTemperature(0), 0.0001);

	// Make sure the diffusion coefficient is 0.0 if E_m is infinite
	cluster.setMigrationEnergy(numeric_limits<double>::infinity());
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0, 0.000001);

	// Make sure the diffusion coefficient is zero if the diffusion factor is zero
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(0.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0, 0.000001);

	// Make sure the diffusion coefficient is equal to the diffusion factor
	// if the temperature is infinite
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	cluster.setTemperature(numeric_limits<double>::infinity(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 1.0, 0.000001);

	// Throw something random in there to be certain
	cluster.setMigrationEnergy(0.013);
	cluster.setDiffusionFactor(1.08E10);
	cluster.setTemperature(1500.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 9766651101.800613,
			0.0000001);

	return;
}

/**
 * This operation tests the default values returned by select flux routines.
 */
BOOST_AUTO_TEST_CASE(checkDefaultFluxes) {
	// Get the simple reaction network
	auto network = getSimpleFeReactionNetwork(0);
	// Create a cluster
	FeCluster cluster(*(network.get()), registry);
	// Add a grid point for the rates
	cluster.addGridPoints(1);

	// Check the default values of the fluxes
	BOOST_REQUIRE_CLOSE(cluster.getProductionFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getCombinationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getDissociationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getEmissionFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getTotalFlux(0), 0.0, 1e-5);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
