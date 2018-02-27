#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <PSIDCluster.h>
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
 * This suite is responsible for testing the PSIDCluster.
 */
BOOST_AUTO_TEST_SUITE(PSIDCluster_testSuite)

/**
 * This operation checks the reaction radius for DCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create a deuterium cluster
	shared_ptr<PSIDCluster> cluster;

	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] =
			{ 0.075, 0.0809312, 0.0850918, 0.0884041, 0.0912012 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<PSIDCluster>(
				new PSIDCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
