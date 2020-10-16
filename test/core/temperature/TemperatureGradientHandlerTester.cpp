#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <mpi.h>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/temperature/GradientHandler.h>

using namespace std;
using namespace xolotl;
using namespace core;

/**
 * The test suite is responsible for testing the GradientHandler.
 */
BOOST_AUTO_TEST_SUITE(TemperatureGradientHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(check_getTemperature)
{
	MPI_Init(NULL, NULL);
	// Create the temperature handler
	auto testTemp = make_shared<temperature::GradientHandler>(1000.0, 900.0);

	// Create a time
	double currTime = 1.0;

	// Create a position
	plsm::SpaceVector<double, 3> x = {0.0, 0.0, 0.0};

	// Get the temperature
	double temp = testTemp->getTemperature(x, currTime);
	BOOST_REQUIRE_CLOSE(temp, 1000.0, 0.001);

	// Get the temperature at a different location
	x[0] = 0.6;
	temp = testTemp->getTemperature(x, currTime);
	BOOST_REQUIRE_CLOSE(temp, 940.0, 0.001);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
