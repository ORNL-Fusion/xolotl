#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <stdio.h>
#include <stdlib.h>

#include <string>

#include <boost/test/included/unit_test.hpp>

#include <xolotl/perf/EventCounter.h>

using namespace std;
using namespace xolotl::perf;

/**
 * This suite is responsible for testing the EventCounter.
 */
BOOST_AUTO_TEST_SUITE(EventCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkInitialValue)
{
	auto tester = EventCounter();

	BOOST_TEST_MESSAGE("\n"
		<< "EventCounter Message: \n"
		<< "tester.getValue() " << tester.getValue() << "\n");

	// Require that the value of this EventCounter is 0
	BOOST_REQUIRE_EQUAL(0U, tester.getValue());
}

BOOST_AUTO_TEST_CASE(checkCounting)
{
	auto tester = EventCounter();

	for (int i = 0; i < 3; i++) {
		// increment the EventCounter
		tester.increment();
	}

	BOOST_TEST_MESSAGE("\n"
		<< "EventCounter Message: \n"
		<< "tester.getValue() = " << tester.getValue() << "\n");

	// Require that the value of this EventCounter is 3
	BOOST_REQUIRE_EQUAL(3U, tester.getValue());
}

BOOST_AUTO_TEST_SUITE_END()
