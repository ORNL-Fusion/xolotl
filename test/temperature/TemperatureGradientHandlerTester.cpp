#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <TemperatureGradientHandler.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the TemperatureGradientHandler.
 */
BOOST_AUTO_TEST_SUITE (TemperatureGradientHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(check_getTemperature) {
	// Create the temperature handler
	auto testTemp = make_shared<TemperatureGradientHandler>(1000.0, 900.0);

	// Create a time
	double currTime = 1.0;

	// Create a position
	Point<3> x = { 0.0, 0.0, 0.0 };

	// Get the temperature
	double temp = testTemp->getTemperature(x, currTime);
	BOOST_REQUIRE_CLOSE(temp, 1000.0, 0.001);

	// Get the temperature at a different location
	x[0] = 0.6;
	temp = testTemp->getTemperature(x, currTime);
	BOOST_REQUIRE_CLOSE(temp, 940.0, 0.001);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
