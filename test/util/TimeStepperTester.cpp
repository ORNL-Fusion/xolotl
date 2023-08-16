#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/util/GrowthFactorStepSequence.h>
#include <xolotl/util/TimeStepper.h>

using namespace xolotl::util;

BOOST_AUTO_TEST_SUITE(StepSequence_testSuite)

BOOST_AUTO_TEST_CASE(growthFactor)
{
	auto stepper =
		TimeStepper(std::make_unique<GrowthFactorStepSequence>(1.0, 10.0, 1.3),
			0.0, 200.0, 50);
	std::vector<double> vals;
	for (stepper.start(); stepper; ++stepper) {
		vals.push_back(stepper.currentTime());
	}
	BOOST_REQUIRE_CLOSE_FRACTION(0.0, vals[0], 0.00001);
	BOOST_REQUIRE_CLOSE_FRACTION(1.0, vals[1], 0.00001);
	BOOST_REQUIRE_CLOSE_FRACTION(2.3, vals[2], 0.00001);
	BOOST_REQUIRE_CLOSE_FRACTION(200.0, vals[26], 0.00001);
}

BOOST_AUTO_TEST_SUITE_END()
