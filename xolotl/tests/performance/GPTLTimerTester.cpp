#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <papi.h>
#include <gptl.h>
#include <GPTLTimer.h>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the GPTLTimer.
 */
BOOST_AUTO_TEST_SUITE (GPTLTimer_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	GPTLinitialize();

	GPTLTimer tester("test");

	BOOST_TEST_MESSAGE( "\n" << "GPTLTimer Message: \n" << "tester.getName() = " << tester.getName() << "\n"
				  );

	//Require that the name of this GPTLTimer is "test"
	BOOST_REQUIRE_EQUAL("test", tester.getName());

    GPTLfinalize();
}


BOOST_AUTO_TEST_CASE(checkTiming) {

	GPTLinitialize();
	GPTLTimer tester("test");
	double sleepSeconds = 2.0;

	//Output the version of PAPI that is being used
	BOOST_TEST_MESSAGE("\n" << "PAPI_VERSION = " << PAPI_VERSION_MAJOR(PAPI_VERSION) << "."
			  << PAPI_VERSION_MINOR(PAPI_VERSION) << "." << PAPI_VERSION_REVISION(PAPI_VERSION) << "\n");


	double wall, usr, sys;
	double wallStart, wallStop;

	//start the timer
	tester.start();

	// GPTLstamp is used here to get the wallclock timestamp when the timer is started
	GPTLstamp(&wall, &usr, &sys);
	BOOST_TEST_MESSAGE("\n" << "Started timer at: " << "wallclock = " << wall << "s");
	wallStart = wall;

	sleep(sleepSeconds);

	//stop the timer
	tester.stop();

	// GPTLstamp is used here to get the wallclock timestamp when the timer is stopped
	GPTLstamp(&wall, &usr, &sys);
	BOOST_TEST_MESSAGE("\n" << "Stopped timer at: " << "wallclock = " << wall << "s");
	wallStop = wall;

	//Output the difference between the wallclock timestamps when the timer was started and stopped
	BOOST_TEST_MESSAGE("\n" << "Difference between wallclock at stop and start: " << wallStop << "s - " << wallStart
			<< "s = " << wallStop - wallStart << "s");

	BOOST_TEST_MESSAGE( "\n" << "GPTLTimer Message: \n" << "tester.getName() = " << tester.getName() << "\n"
			  << "tester.getValue() = " << tester.getValue() << "s" << "\n"
			  << "tester.getValue() - " << sleepSeconds << "s = " << tester.getValue()-sleepSeconds << "s");

	//Require that the value of this GPTLTimer is within 3% of the value of sleepSeconds
	BOOST_REQUIRE_CLOSE(sleepSeconds, tester.getValue(),0.03);

	BOOST_REQUIRE_EQUAL("s", tester.getUnits());

    GPTLfinalize();
}

BOOST_AUTO_TEST_SUITE_END()
