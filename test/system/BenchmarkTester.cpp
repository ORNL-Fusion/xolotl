#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <Kokkos_Core.hpp>
using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

#include <xolotl/test/MPITestUtils.h>
#include <xolotl/test/SystemTestCase.h>
using xolotl::test::SystemTestCase;
using xolotl::test::getMPICommSize;

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_SUITE(Benchmark_testSuite)

BOOST_AUTO_TEST_CASE(benchmark_NE_9999DOF)
{
    if (getMPICommSize() > 1) {
        return;
    }
	SystemTestCase{"benchmark_NE_9999DOF"}.withTimer().run();
}

BOOST_AUTO_TEST_CASE(benchmark_NE_20079DOF)
{
    if (getMPICommSize() > 1) {
        return;
    }
	SystemTestCase{"benchmark_NE_20079DOF"}.withTimer().run();
}

BOOST_AUTO_TEST_SUITE_END()
