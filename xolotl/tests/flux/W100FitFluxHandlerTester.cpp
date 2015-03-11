#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W100FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W100FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.476819, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.225961, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.097220, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkHeFluence) {
	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Check that the fluence is 0.0 at the beginning
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 0.0);

	// Increment the helium fluence
	testFitFlux->incrementHeFluence(1.0e-8);
	
	// Check that the fluence is not 0.0 anymore
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 1.0e-8);

	return;
}

BOOST_AUTO_TEST_CASE(checkHeFlux) {
	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Set the factor to change the helium flux
    testFitFlux->setHeFlux(2.5);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

    // Check the value of the helium flux
    BOOST_REQUIRE_EQUAL(testFitFlux->getHeFlux(), 2.5);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 1.192047, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.564902, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.243050, 0.01);

	return;
}


BOOST_AUTO_TEST_SUITE_END()
