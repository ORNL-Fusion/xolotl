#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyHardwareCounter.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

const std::vector<HardwareQuantities> test_hardwareQuantities = {L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,
			BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

/**
 * This suite is responsible for testing the DummyHardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (DummyHardwareCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	DummyHardwareCounter tester("test",test_hardwareQuantities);

	BOOST_REQUIRE_EQUAL("unused", tester.getName());
}

BOOST_AUTO_TEST_CASE(check_getValues) {

	DummyHardwareCounter tester("test",test_hardwareQuantities);

	BOOST_TEST_MESSAGE("\n" << "DummyHardwareCounter Message: \n"
							<< "tester.getValues().size() = " << tester.getValues().size());

	BOOST_REQUIRE_EQUAL(0, tester.getValues().size());

}

BOOST_AUTO_TEST_CASE(check_getHardwareQuantities) {

	DummyHardwareCounter tester("test",test_hardwareQuantities);

	BOOST_TEST_MESSAGE("\n" << "DummyHardwareCounter Message: \n"
							<< "tester.getHardwareQuantities().size() = " << tester.getHardwareQuantities().size());

	BOOST_REQUIRE_EQUAL(0, tester.getHardwareQuantities().size());

}

BOOST_AUTO_TEST_SUITE_END()

