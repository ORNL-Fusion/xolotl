#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <SystemTest.h>

BOOST_AUTO_TEST_SUITE(System_testSuite)

BOOST_AUTO_TEST_CASE(system_NE_0) {
    testUtils::runSystemTestCase("system_NE_0");
}

BOOST_AUTO_TEST_CASE(system_NE_1) {
	// 0D + heterogeneous nucleation
    testUtils::runSystemTestCase("system_NE_1");
}

BOOST_AUTO_TEST_CASE(system_NE_2) {
	// 0D + grouping + re-solution
    testUtils::runSystemTestCase("system_NE_2");
}

BOOST_AUTO_TEST_CASE(system_NE_3) {
	// 1D + temperature gradient
    testUtils::runSystemTestCase("system_NE_3");
}

BOOST_AUTO_TEST_CASE(system_NE_4) {
	// 2D
    testUtils::runSystemTestCase("system_NE_4");
}

BOOST_AUTO_TEST_CASE(system_NE_5) {
	// 3D
    testUtils::runSystemTestCase("system_NE_5");
}

BOOST_AUTO_TEST_CASE(system_PSI_1) {
	// 1D + HeVI + advection + modifiedTM + attenuation + surface + reflective bulk
    testUtils::runSystemTestCase("system_PSI_1");
}

BOOST_AUTO_TEST_CASE(system_PSI_2) {
    BOOST_REQUIRE(testUtils::copyFile("tridyn.dat"));
	// 1D + HeDTVI + grouping + heat
    testUtils::runSystemTestCase("system_PSI_2");
}

BOOST_AUTO_TEST_SUITE_END()
