#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/util/GrowthFactorStepSequence.h>

using namespace xolotl::util;

BOOST_AUTO_TEST_SUITE(StepSequence_testSuite)

BOOST_AUTO_TEST_CASE(growthFactor)
{
	auto seq = GrowthFactorStepSequence(1.0, 10.0, 1.2);
	seq.start();
	std::vector<double> vals;
	for (std::size_t i = 0; i < 20; ++i) {
		vals.push_back(seq.current());
		seq.step();
	}
	BOOST_REQUIRE_CLOSE_FRACTION(1.0, vals[0], 0.00001);
	BOOST_REQUIRE_CLOSE_FRACTION(1.2, vals[1], 0.00001);
	BOOST_REQUIRE_CLOSE_FRACTION(1.44, vals[2], 0.00001);
	BOOST_REQUIRE_CLOSE_FRACTION(10.0, vals[13], 0.00001);
}

BOOST_AUTO_TEST_SUITE_END()
