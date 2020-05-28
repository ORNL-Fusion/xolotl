#define BOOST_TEST_MODULE Regression

#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <boost/test/included/unit_test.hpp>
#include "papi.h"
#include "xolotlPerf/papi/PAPIHardwareCounter.h"
#include "tests/performance/perfTestConfig.h"

using namespace xolotlPerf;

const IHardwareCounter::SpecType test_ctrSpec = {
		IHardwareCounter::Instructions, IHardwareCounter::Cycles,
		IHardwareCounter::FPOps, IHardwareCounter::L3CacheMisses };

// Normally, PAPI would be initialized by the HandlerRegistry.
// Since our purpose is to test the Timer class and not the registry,
// we recreate the initialization explicitly.
bool InitializePAPI(void) {
	bool ret = true;

	if (!PAPI_is_initialized()) {
		int papiVersion = PAPI_library_init(PAPI_VER_CURRENT);
		if (papiVersion != PAPI_VER_CURRENT) {
			BOOST_TEST_MESSAGE(
					"PAPI library version mismatch: asked for" << PAPI_VER_CURRENT << ", got " << papiVersion);
			ret = false;
		}
	}
	return ret;
}

/**
 * This suite is responsible for testing the HardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (HardwareCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	bool papiInitialized = InitializePAPI();
	BOOST_REQUIRE_EQUAL(papiInitialized, true);

	PAPIHardwareCounter tester("test", test_ctrSpec);

	//Output the version of PAPI that is being used
	BOOST_TEST_MESSAGE(
			"\n" << "PAPI_VERSION = " << PAPI_VERSION_MAJOR(PAPI_VERSION) << "." << PAPI_VERSION_MINOR(PAPI_VERSION) << "." << PAPI_VERSION_REVISION(PAPI_VERSION) << "\n");

	BOOST_REQUIRE_EQUAL("test", tester.getName());
}

BOOST_AUTO_TEST_CASE(check_getSpecification) {

	bool papiInitialized = InitializePAPI();
	BOOST_REQUIRE_EQUAL(papiInitialized, true);

	PAPIHardwareCounter tester("test", test_ctrSpec);

	BOOST_TEST_MESSAGE(
			"\n" << "PAPIHardwareCounter Message: \n" << "ctrSpec = ");
	const IHardwareCounter::SpecType& ctrSpec = tester.getSpecification();
	for (auto iter = ctrSpec.begin(); iter != ctrSpec.end(); ++iter) {
		BOOST_TEST_MESSAGE(" " << tester.getCounterName(*iter) << " ");
	}

	BOOST_REQUIRE_EQUAL(test_ctrSpec.size(), ctrSpec.size());
	if (test_ctrSpec.size() == ctrSpec.size()) {
		for (unsigned i = 0; i < test_ctrSpec.size(); i++) {
			BOOST_REQUIRE_EQUAL(test_ctrSpec[i], ctrSpec[i]);
		}
	}
}

BOOST_AUTO_TEST_CASE(check_getValues) {
	bool papiInitialized = InitializePAPI();
	BOOST_REQUIRE_EQUAL(papiInitialized, true);

	PAPIHardwareCounter tester("test", test_ctrSpec);

	const unsigned int nMultiplies = 1000;
	tester.start();
	double a = 2;
	for (unsigned int i = 0; i < nMultiplies; i++) {
		a *= ((double) (i + 1));
	}
	tester.stop();

	// Use variable 'a' so that compiler can't optimize the above loop away
	BOOST_TEST_MESSAGE(
			"\nPAPIHardwareCounter test:\n" << "produced value: " << a << '\n');

	// Output the counts we measured.
	auto testVals = tester.getValues();
	BOOST_TEST_MESSAGE("produced " << testVals.size() << " values");
	BOOST_REQUIRE_EQUAL(testVals.size(), test_ctrSpec.size());

	std::ostringstream mstr;
	for (unsigned int i = 0; i < testVals.size(); ++i) {
		mstr << " " << tester.getCounterName(test_ctrSpec[i]) << ": "
				<< testVals[i] << '\n';
	}
	BOOST_TEST_MESSAGE("PAPIHardwareCounter measured:\n" << mstr.str());

	// Verify we got what we expected.
	// It isn't easy to verify instructions or cycles - these
	// are compiler- and processor-dependent.
	// However, we know how many floating point operations we
	// did in our test loop.
	// NOTE: the indices into testVals must change if you change
	// the test_ctrSpec specification.
	std::ostringstream vstr;
	vstr << "We believe:\n" << "* this test program was compiled "
#if defined(XOLOTL_TEST_HWCTR_DEBUGEXP)
			<< "without"
#else
			<< "with"
#endif // defined(XOLOTL_HWCTR_DEBUGEXP)
			<< " compiler optimizations enabled.\n"
			<< "* the hardware counter should measure approximately "
			<< XOLOTL_HWCTR_EXPVAL << " FP ops +/-" << XOLOTL_HWCTR_TOL
			<< "%\nwhen executing a multiplication statement " << nMultiplies
			<< " times in a loop.";
	BOOST_TEST_MESSAGE(vstr.str());

	BOOST_REQUIRE_CLOSE(testVals[2], (double)XOLOTL_HWCTR_EXPVAL,
			XOLOTL_HWCTR_TOL);
}

BOOST_AUTO_TEST_SUITE_END()
