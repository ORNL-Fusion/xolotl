#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "FeFitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the FeFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FeFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Specify the number of grid points that will be used
	int nGridpts = 21;
	// Specify the step size between grid points
	double step = 1.05;

    auto testFitFlux = make_shared<FeFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.116361, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[5], 0.098120, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[15], 0.0, 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
