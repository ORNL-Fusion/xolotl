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
#include <memory>
#include <typeinfo>
#include <limits>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the PSIClusterNetworkLoader. It
 * creates a string stream that contains each of the available PSICluster types
 * and checks that the loader returns a list with each type in it.
 */
BOOST_AUTO_TEST_SUITE (PSIClusterNetworkLoader_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkDiffusionCoefficient) {

	// Local Declarations
	double temperature = 1.0;
	PSICluster cluster(1);

	// Check E_m = 0.0
	cluster.setMigrationEnergy(0.0);
	cluster.setDiffusionFactor(1.0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(1.0),exp(0.0),0.00001);

	// Make sure the diffusion coefficient is 0.0 if E_m is infinite
	cluster.setMigrationEnergy(std::numeric_limits<double>::infinity());
	cluster.setDiffusionFactor(1.0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(1.0),0.0,0.00000);

	// Make sure the diffusion coefficient is zero if the diffusion factor is zero
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(0.0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(1.0),0.0,0.00000);

	// Make sure the diffusion coefficient is equal to the diffusion factor
	// if the temperature is infinite
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(1.0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(
					std::numeric_limits<double>::infinity()),1.0,0.00000);

	// Throw something random in there to be certain
	cluster.setMigrationEnergy(0.013);
	cluster.setDiffusionFactor(1.08E10);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(1500.0),9766651101.800613,0.0000001);

}

BOOST_AUTO_TEST_SUITE_END()

