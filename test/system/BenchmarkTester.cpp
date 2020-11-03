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

BOOST_AUTO_TEST_SUITE(Benchmark)

BOOST_AUTO_TEST_CASE(NE_1)
{
    if (getMPICommSize() > 1) {
        return;
    }
    // 0D, 9999 DOF ungrouped
	SystemTestCase{"benchmark_NE_1"}.withTimer().run();
}

BOOST_AUTO_TEST_CASE(NE_2)
{
    if (getMPICommSize() > 1) {
        return;
    }
    // 0D, 20079 DOF grouped
	SystemTestCase{"benchmark_NE_2"}.withTimer().run();
}

BOOST_AUTO_TEST_CASE(NE_3)
{
	if (getMPICommSize() > 400) {
		return;
    // 3D
	}
	SystemTestCase{"benchmark_NE_3"}.withTimer().run();
}

BOOST_AUTO_TEST_CASE(PSI_1)
{
	if (getMPICommSize() > 10) {
		return;
	}
	SystemTestCase::copyFile("tridyn_benchmark_PSI_1.dat");
	// 1D + HeDTVI + grouping + heat
	SystemTestCase{"benchmark_PSI_1"}.withTimer().run();
}

BOOST_AUTO_TEST_CASE(PSI_2)
{
	if (getMPICommSize() > 20) {
		return;
	}
	SystemTestCase::copyFile("tridyn_benchmark_PSI_2.dat");
	// 1D + HeDVI + advection + modifiedTM + attenuation + surface + reflective
	SystemTestCase{"benchmark_PSI_2"}.withTimer().run();
}

BOOST_AUTO_TEST_CASE(PSI_3)
{
	if (getMPICommSize() > 20) {
		return;
	}
	// 1D + HeVI + advection + modifiedTM + attenuation + surface + reflective +
	// reduced matrix method bulk
	SystemTestCase{"benchmark_PSI_3"}.withTimer().run();
}

BOOST_AUTO_TEST_SUITE_END()
