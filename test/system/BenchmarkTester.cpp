#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <xolotl/test/SystemTestCase.h>
using xolotl::test::SystemTestCase;

BOOST_AUTO_TEST_SUITE(Benchmark_testSuite)

BOOST_AUTO_TEST_CASE(NE_9999DOF_plsm)
{
	SystemTestCase{"NE_9999DOF_plsm"}.run();
}

BOOST_AUTO_TEST_CASE(NE_20079DOF_plsm)
{
	SystemTestCase{"NE_20079DOF_plsm"}.run();
}

BOOST_AUTO_TEST_SUITE_END()
