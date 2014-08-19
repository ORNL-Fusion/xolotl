#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "WFitFluxHandler.h"

using namespace std;
using namespace xolotlSolver;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

    auto testFitFlux = std::make_shared<xolotlSolver::WFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create a composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	std::vector<double> x = {1.25, 0.0, 0.0};

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);

	BOOST_TEST_MESSAGE( "\nWFitFluxHandlerTester Message: \n"
						<< "incidentFlux = " << testFlux << " with position "
						<< "(" << x[0] << "," << x[1] << "," << x[2] << ") "
						<< "at time = " << currTime << "\n");
	BOOST_REQUIRE_CLOSE(testFlux, 0.449045, 0.01);

}

BOOST_AUTO_TEST_CASE(checkHeFluence) {

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

    auto testFitFlux = std::make_shared<xolotlSolver::WFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create the composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	std::vector<double> x = {1.25, 0.0, 0.0};

	// Check the flux
	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);
	BOOST_REQUIRE_CLOSE(testFlux, 0.449045, 0.01);

	// Set the maximum helium fluence value
	testFitFlux->setMaxHeFluence(8.0e-09);
	BOOST_REQUIRE_EQUAL(testFitFlux->getMaxHeFluence(), 8.0e-09);

	// Increment the helium fluence
	testFitFlux->incrementHeFluence(1.0e-8);
	auto fluence1 = testFitFlux->getHeFluence();
	// Increment the helium fluence again, although the helium fluence value should not change
	// because it already reached the maximum
	testFitFlux->incrementHeFluence(1.0e-7);
	// Check to make sure the helium fluence value did not increment a second time
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), fluence1);
	// Since the maximum helium fluence value has been exceeded, an incident flux of 0 should be returned
	BOOST_REQUIRE_EQUAL(testFitFlux->getIncidentFlux(compVec, x, 1), 0.0);

}

BOOST_AUTO_TEST_CASE(checkHeFlux) {

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

    auto testFitFlux = std::make_shared<xolotlSolver::WFitFluxHandler>();
    // Set the factor to change the Helium flux
    testFitFlux->setHeFlux(2.5);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

    BOOST_REQUIRE_EQUAL(testFitFlux->getHeFlux(), 2.5);

	// Create a composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	std::vector<double> x = {1.25, 0.0, 0.0};

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);
	BOOST_REQUIRE_CLOSE(testFlux, 2.5 * 0.449045, 0.01);
}


BOOST_AUTO_TEST_SUITE_END()
