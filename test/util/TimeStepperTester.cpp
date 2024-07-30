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
	std::vector<double> timeVals;
	std::vector<double> stepVals;
	for (stepper.start(); stepper; ++stepper) {
		timeVals.push_back(stepper.currentTime());
		stepVals.push_back(stepper.currentTimeStepSize());
	}
	BOOST_REQUIRE_EQUAL(timeVals.size(), 26);
	BOOST_REQUIRE_CLOSE_FRACTION(0.0, timeVals[0], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(1.0, timeVals[1], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(2.3, timeVals[2], 1e-6);
	BOOST_REQUIRE(200.0 <= stepper.currentTime());

	auto agreeTol = std::numeric_limits<double>::epsilon();

	BOOST_REQUIRE_CLOSE_FRACTION(timeVals[0], stepper.timeAtStep(0), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(timeVals[1], stepper.timeAtStep(1), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(timeVals[3], stepper.timeAtStep(3), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		timeVals[24], stepper.timeAtStep(24), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[0], stepper.timeStepSizeAtStep(0), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[1], stepper.timeStepSizeAtStep(1), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[3], stepper.timeStepSizeAtStep(3), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[24], stepper.timeStepSizeAtStep(24), agreeTol);

	// Reach max steps first
	stepper = TimeStepper(
		std::make_unique<GrowthFactorStepSequence>(*seq), 10.0, 200.0, 20);
	timeVals.clear();
	for (stepper.start(); stepper; ++stepper) {
		timeVals.push_back(stepper.currentTime());
		stepVals.push_back(stepper.currentTimeStepSize());
	}
	BOOST_REQUIRE(timeVals.size() == 21);
	BOOST_REQUIRE_CLOSE_FRACTION(10.0, timeVals[0], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(11.0, timeVals[1], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(12.3, timeVals[2], 1e-6);
	BOOST_REQUIRE_CLOSE_FRACTION(152.015, timeVals[20], 1e-6);

	BOOST_REQUIRE_CLOSE_FRACTION(timeVals[0], stepper.timeAtStep(0), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(timeVals[1], stepper.timeAtStep(1), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(timeVals[3], stepper.timeAtStep(3), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		timeVals[19], stepper.timeAtStep(19), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[0], stepper.timeStepSizeAtStep(0), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[1], stepper.timeStepSizeAtStep(1), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[3], stepper.timeStepSizeAtStep(3), agreeTol);
	BOOST_REQUIRE_CLOSE_FRACTION(
		stepVals[19], stepper.timeStepSizeAtStep(19), agreeTol);

	stepper = TimeStepper(
		std::make_unique<GrowthFactorStepSequence>(1.0e-12, 1.0, 1.1), 0.0,
		5000.0, 10000);
	for (stepper.start(); stepper; ++stepper) { }
	BOOST_REQUIRE(stepper.currentStep() == 5280);
	BOOST_REQUIRE(stepper.currentTime() >= 5000);
}

BOOST_AUTO_TEST_SUITE_END()
