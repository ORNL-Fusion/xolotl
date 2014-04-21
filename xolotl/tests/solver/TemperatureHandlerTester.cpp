#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "TemperatureHandler.h"

using namespace std;
using namespace xolotlSolver;

/*
 * The test suite is responsible for testing the TemperatureHandler.
 */
BOOST_AUTO_TEST_SUITE (TemperatureHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(check_getTemperature) {

	auto testTemp = std::make_shared<xolotlSolver::TemperatureHandler>();

	// Create a time
	double currTime = 1.0;

	std::vector<double> x = {0.0, 1.142857142857143, 0.0};
	// x is a gridpoint in PetscSolver, RHSFunction
	//x=1.142857142857143=8/7 where x = xi * hx with xi=1, hx=1.142857142857143

	auto temp = testTemp->getTemperature(x, 1);

	BOOST_TEST_MESSAGE( "\n" << "\nTemperatureHandlerTester Message: \n"
						<< "temperature = " << temp << " at "
						<< "(" << x[0] << "," << x[1] << "," << x[2] << ") "
						<< "at time = " << currTime << "\n");
	BOOST_REQUIRE_EQUAL(temp, 1000.0);

}

BOOST_AUTO_TEST_SUITE_END()
