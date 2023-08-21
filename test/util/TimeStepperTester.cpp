#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/util/GrowthFactorStepSequence.h>
#include <xolotl/util/TimeStepper.h>

using namespace xolotl::util;

BOOST_AUTO_TEST_SUITE(StepSequence_testSuite)

BOOST_AUTO_TEST_CASE(growthFactor)
{
	// Reach end time first
	auto seq = std::make_unique<GrowthFactorStepSequence>(1.0, 10.0, 1.3);
	auto stepper = TimeStepper(
		std::make_unique<GrowthFactorStepSequence>(*seq), 0.0, 200.0, 50);
	std::vector<double> vals;
	for (stepper.start(); stepper; ++stepper) {
		vals.push_back(stepper.currentTime());
	}
	BOOST_REQUIRE(vals.size() == 27);
	BOOST_REQUIRE_CLOSE_FRACTION(0.0, vals[0], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(1.0, vals[1], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(2.3, vals[2], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(200.0, vals[26], 1e-6);

	// Reach max steps first
	stepper = TimeStepper(
		std::make_unique<GrowthFactorStepSequence>(*seq), 10.0, 200.0, 20);
	vals.clear();
	for (stepper.start(); stepper; ++stepper) {
		vals.push_back(stepper.currentTime());
	}
	BOOST_REQUIRE(vals.size() == 21);
	BOOST_REQUIRE_CLOSE_FRACTION(10.0, vals[0], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(11.0, vals[1], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(12.3, vals[2], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(152.015, vals[20], 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
