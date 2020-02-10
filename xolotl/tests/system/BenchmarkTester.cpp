#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <SystemTest.h>

BOOST_AUTO_TEST_SUITE(Benchmark_testSuite)

BOOST_AUTO_TEST_CASE(NE_9999DOF_plsm) {
    testUtils::runSystemTestCase("NE_9999DOF_plsm");
}

BOOST_AUTO_TEST_SUITE_END()
