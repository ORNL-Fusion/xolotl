#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <SystemTest.h>

BOOST_AUTO_TEST_SUITE(System_testSuite)

BOOST_AUTO_TEST_CASE(system_NE_1) {
    testUtils::runSystemTestCase("system_NE_1");
}

BOOST_AUTO_TEST_CASE(system_NE_2) {
    testUtils::runSystemTestCase("system_NE_2");
}

BOOST_AUTO_TEST_SUITE_END()
