#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <SystemTestCase.h>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Benchmark_testSuite)

BOOST_AUTO_TEST_CASE(NE_9999DOF_plsm)
{
    xolotl::test::SystemTestCase{"NE_9999DOF_plsm"}();
}

BOOST_AUTO_TEST_CASE(NE_20079DOF_plsm)
{
    xolotl::test::SystemTestCase{"NE_20079DOF_plsm"}();
}

BOOST_AUTO_TEST_SUITE_END()
