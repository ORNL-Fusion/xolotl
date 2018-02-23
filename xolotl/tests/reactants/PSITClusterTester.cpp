#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <PSITCluster.h>
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
 * This suite is responsible for testing the TCluster.
 */
BOOST_AUTO_TEST_SUITE(PSITCluster_testSuite)

/**
 * This operation checks the reaction radius for TCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create a tritium cluster
	shared_ptr<PSITCluster> cluster;

	// Get the simple reaction network
	auto network = getSimplePSIReactionNetwork(0);

	// The vector of radii to compare with
	double expectedRadii[] = { 0.3, 0.3, 0.3, 0.3,
			0.3 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<PSITCluster>(
				new PSITCluster(i, *(network.get()), registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
