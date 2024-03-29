#define BOOST_TEST_MODULE Regression

#include <string>

#include <boost/test/included/unit_test.hpp>

#include <xolotl/perf/os/OSTimer.h>

using namespace std;
using namespace xolotl;

/**
 * This suite is responsible for testing the OSTimer.
 */
BOOST_AUTO_TEST_SUITE(OSTimer_testSuite)

BOOST_AUTO_TEST_CASE(checkTiming)
{
	auto tester = perf::os::OSTimer();
	double sleepSeconds = 2.0;

	// Simulate some computation/communication with a sleep of known duration.
	// Time the duration of the operation.
	tester.start();
	sleep(sleepSeconds);
	tester.stop();

	// Require that the value of this Timer is within 3% of the
	// duration of the sleep.
	BOOST_REQUIRE_CLOSE(sleepSeconds, tester.getValue(), 0.03);
	BOOST_REQUIRE_EQUAL("s", tester.getUnits());
}

BOOST_AUTO_TEST_CASE(accumulate)
{
	auto tester = perf::os::OSTimer();

	const unsigned int sleepSeconds = 2;

	tester.start();
	sleep(sleepSeconds);
	tester.stop();
	tester.start();
	sleep(sleepSeconds);
	tester.stop();

	double timerValue = tester.getValue();
	double expValue = 2 * sleepSeconds; // we had two sleep intervals
	BOOST_REQUIRE_CLOSE(expValue, timerValue, 0.03);
}

BOOST_AUTO_TEST_CASE(reset)
{
	auto tester = perf::os::OSTimer();

	const unsigned int sleepSeconds = 2;

	tester.start();
	sleep(sleepSeconds);
	tester.stop();
	tester.reset();
	BOOST_REQUIRE_EQUAL(tester.getValue(), 0.0);
	tester.start();
	sleep(sleepSeconds);
	tester.stop();

	double timerValue = tester.getValue();
	double expValue = sleepSeconds; // should only represent last sleep interval
	BOOST_REQUIRE_CLOSE(expValue, timerValue, 0.03);
}

BOOST_AUTO_TEST_SUITE_END()
