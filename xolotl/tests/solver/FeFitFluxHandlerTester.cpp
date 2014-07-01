#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "FeFitFluxHandler.h"

using namespace std;
using namespace xolotlSolver;

/**
 * The test suite is responsible for testing the FeFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FeFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {

	// Specify the number of grid points that will be used
	int nGridpts = 15;
	// Specify the step size between grid points
	double step = 0.5714285714285714;

    auto testFitFlux = std::make_shared<xolotlSolver::FeFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(nGridpts, step);

	// Create a composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	// Create a vector representing the position of the cluster
	std::vector<double> x = {2.857142857142857, 0.0, 0.0};
	std::vector<double> x1 = {4, 0.0, 0.0};

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);
	auto testFlux1 = testFitFlux->getIncidentFlux(compVec, x1, 1);

	BOOST_TEST_MESSAGE( "\n" << "\nFeFitFluxHandlerTester Message: \n"
						<< "incidentFlux = " << testFlux << " with composition "
						<< "(" << compVec[0] << "," << compVec[1] << "," << compVec[2] << "), "
						<< "at position " << "(" << x[0] << "," << x[1] << "," << x[2] << "), "
						<< "at time = " << currTime << "\n");
	BOOST_REQUIRE_EQUAL(testFlux, 1.0394067746526634);

	BOOST_TEST_MESSAGE( "\n" << "incidentFlux = " << testFlux1 << " with composition "
						<< "(" << compVec[0] << "," << compVec[1] << "," << compVec[2] << "), "
						<< "at position " << "(" << x[0] << "," << x[1] << "," << x[2] << "), "
						<< "at time = " << currTime << "\n");
	BOOST_REQUIRE_EQUAL(testFlux1, 0.9887042553323879);

}


BOOST_AUTO_TEST_SUITE_END()
