#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <DCluster.h>
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
 * This suite is responsible for testing the DCluster.
 */
BOOST_AUTO_TEST_SUITE(DCluster_testSuite)

/**
 * This operation checks the reaction radius for DCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create a deuterium cluster
	shared_ptr<DCluster> cluster;

	// The vector of radii to compare with
	double expectedRadii[] = { 0.012407, 0.0156319, 0.017894, 0.0196949,
			0.0212157 };

	// Check all the values
	for (int i = 1; i <= 5; i++) {
		cluster = shared_ptr<DCluster>(new DCluster(i, registry));
		BOOST_REQUIRE_CLOSE(expectedRadii[i - 1], cluster->getReactionRadius(),
				0.001);
	}

	return;
}

BOOST_AUTO_TEST_SUITE_END()
