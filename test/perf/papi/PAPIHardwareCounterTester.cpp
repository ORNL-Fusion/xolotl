#define BOOST_TEST_MODULE Regression

#include <papi.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

#include <boost/test/included/unit_test.hpp>

#include <xolotl/perf/papi/PAPIHardwareCounter.h>
#include <xolotl/test/perf/config.h>

using namespace xolotl::perf;
using namespace papi;

struct PAPIFixture
{
	PAPIFixture()
	{
		if (!PAPI_is_initialized()) {
			int papiVersion = PAPI_library_init(PAPI_VER_CURRENT);
			if (papiVersion != PAPI_VER_CURRENT) {
				BOOST_TEST_MESSAGE("PAPI library version mismatch: asked for"
					<< PAPI_VER_CURRENT << ", got " << papiVersion);
				throw std::runtime_error("Unable to construct PAPI fixture");
			}
		}
	}
};

BOOST_GLOBAL_FIXTURE(PAPIFixture);

IHardwareCounter::SpecType
getCtrSpec()
{
	using SpecInfo = std::pair<IHardwareCounter::CounterSpec, int>;
	std::vector<SpecInfo> specs = {
		{IHardwareCounter::Instructions, PAPI_TOT_INS},
		{IHardwareCounter::Cycles, PAPI_TOT_CYC},
		{IHardwareCounter::FPOps, PAPI_FP_OPS},
		{IHardwareCounter::L3CacheMisses, PAPI_L3_TCM}};

	IHardwareCounter::SpecType ret;

	for (const auto& spec : specs) {
		if (PAPI_query_event(spec.second) == PAPI_OK) {
			ret.push_back(spec.first);
		}
	}

	return ret;
}

/**
 * This suite is responsible for testing the HardwareCounter.
 */
BOOST_AUTO_TEST_SUITE(HardwareCounter_testSuite)

BOOST_AUTO_TEST_CASE(check_getSpecification)
{
	PAPIHardwareCounter tester(getCtrSpec());

	BOOST_TEST_MESSAGE("\n"
		<< "PAPIHardwareCounter Message: \n"
		<< "ctrSpec = ");
	const IHardwareCounter::SpecType& ctrSpec = tester.getSpecification();
	for (auto iter = ctrSpec.begin(); iter != ctrSpec.end(); ++iter) {
		BOOST_TEST_MESSAGE(" " << tester.getCounterName(*iter) << " ");
	}

	auto testCtrSpec = getCtrSpec();
	BOOST_REQUIRE_EQUAL(testCtrSpec.size(), ctrSpec.size());
	if (testCtrSpec.size() == ctrSpec.size()) {
		for (unsigned i = 0; i < testCtrSpec.size(); i++) {
			BOOST_REQUIRE_EQUAL(testCtrSpec[i], ctrSpec[i]);
		}
	}
}

BOOST_AUTO_TEST_CASE(check_getValues)
{
	auto testCtrSpec = getCtrSpec();
	PAPIHardwareCounter tester(testCtrSpec);

	const unsigned int nMultiplies = 1000;
	tester.start();
	double a = 2;
	for (unsigned int i = 0; i < nMultiplies; i++) {
		a *= ((double)(i + 1));
	}
	tester.stop();

	// Use variable 'a' so that compiler can't optimize the above loop away
	BOOST_TEST_MESSAGE("\nPAPIHardwareCounter test:\n"
		<< "produced value: " << a << '\n');

	// Output the counts we measured.
	auto testVals = tester.getValues();
	BOOST_TEST_MESSAGE("produced " << testVals.size() << " values");
	BOOST_REQUIRE_EQUAL(testVals.size(), testCtrSpec.size());

	std::ostringstream mstr;
	for (unsigned int i = 0; i < testVals.size(); ++i) {
		mstr << " " << tester.getCounterName(testCtrSpec[i]) << ": "
			 << testVals[i] << '\n';
	}
	BOOST_TEST_MESSAGE("PAPIHardwareCounter measured:\n" << mstr.str());

	// if FP ops are not supported, quit here
	if (std::find(testCtrSpec.begin(), testCtrSpec.end(),
			IHardwareCounter::FPOps) == testCtrSpec.end()) {
		return;
	}

	// Verify we got what we expected.
	// It isn't easy to verify instructions or cycles - these
	// are compiler- and processor-dependent.
	// However, we know how many floating point operations we
	// did in our test loop.
	// NOTE: the indices into testVals must change if you change
	// the testCtrSpec specification.
	std::ostringstream vstr;
	vstr << "We believe:\n"
		 << "* this test program was compiled "
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

	BOOST_REQUIRE_CLOSE(
		testVals[2], (double)XOLOTL_HWCTR_EXPVAL, XOLOTL_HWCTR_TOL);
}

BOOST_AUTO_TEST_SUITE_END()
