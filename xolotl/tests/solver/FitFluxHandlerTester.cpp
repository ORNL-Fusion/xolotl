#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "FitFluxHandler.h"

using namespace std;
using namespace xolotlSolver;

/**
 * The test suite is responsible for testing the FitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {

    auto testFitFlux = std::make_shared<xolotlSolver::FitFluxHandler>();

	// Create a composition vector
	std::vector<int> compVec = {1, 0, 0};
	// Create a time
	double currTime = 1.0;

	std::vector<double> x = {0.0, 1.142857142857143, 0.0};
	// x is a gridpoint in PetscSolver, RHSFunction
	//x=1.142857142857143=8/7 where x = xi * hx with xi=1, hx=1.142857142857143
	double fitFunction = 0.0006 * x[1] * x[1] * x[1] - 0.0087 * x[1] * x[1] + 0.0300 * x[1];
	//fitFunction = 0.02381807580174927

	auto testFlux = testFitFlux->getIncidentFlux(compVec, x, 1);

	BOOST_TEST_MESSAGE( "\n" << "\nFitFluxHandlerTester Message: \n"
						<< "incidentFlux = " << testFlux << " with composition "
						<< "(" << x[0] << "," << x[1] << "," << x[2] << ") "
						<< "at time = " << currTime << "\n"
						<< "fitFunction = " << fitFunction << "\n");
	BOOST_REQUIRE_EQUAL(testFlux, fitFunction);

}

BOOST_AUTO_TEST_SUITE_END()
