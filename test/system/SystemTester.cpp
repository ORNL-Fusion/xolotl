#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>

#include <Kokkos_Core.hpp>
using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

#include <xolotl/test/MPITestUtils.h>
#include <xolotl/test/SystemTestCase.h>
using xolotl::test::getMPICommSize;
using xolotl::test::SystemTestCase;

BOOST_GLOBAL_FIXTURE(MPIFixture);

using xolotl::test::SystemTestOptions;
BOOST_GLOBAL_FIXTURE(SystemTestOptions);

BOOST_AUTO_TEST_SUITE(System_testSuite)

BOOST_AUTO_TEST_CASE(system_NE_0)
{
	if (getMPICommSize() > 1) {
		return;
	}
	SystemTestCase{"system_NE_0"}.tolerance(1.0e-9).run();
}

BOOST_AUTO_TEST_CASE(system_NE_1)
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D + heterogeneous nucleation
	SystemTestCase{"system_NE_1"}.run();
}

BOOST_AUTO_TEST_CASE(system_NE_2)
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D + grouping + re-solution
	SystemTestCase{"system_NE_2"}.run();
}

BOOST_AUTO_TEST_CASE(system_NE_3)
{
	if (getMPICommSize() > 25) {
		return;
	}
	// 1D + temperature gradient
	SystemTestCase{"system_NE_3"}.run();
}

BOOST_AUTO_TEST_CASE(system_NE_4)
{
	if (getMPICommSize() > 400) {
		return;
    // 3D
	}
	SystemTestCase{"system_NE_4"}.run();
}

BOOST_AUTO_TEST_CASE(system_NE_5)
{
	if (getMPICommSize() > 100) {
		return;
	}
	// 2D
	SystemTestCase{"system_NE_5"}.run();
}

BOOST_AUTO_TEST_CASE(system_PSI_1)
{
	if (getMPICommSize() > 20) {
		return;
	}
	// 1D + HeVI + advection + modifiedTM + attenuation + surface + reflective
	// bulk
	SystemTestCase{"system_PSI_1"}.run();
}

BOOST_AUTO_TEST_CASE(system_PSI_2)
{
	if (getMPICommSize() > 10) {
		return;
	}
	SystemTestCase::copyFile("tridyn_system_PSI_2.dat");
	// 1D + HeDTVI + grouping + heat
	SystemTestCase{"system_PSI_2"}.run();
}

BOOST_AUTO_TEST_CASE(system_PSI_3)
{
	if (getMPICommSize() > 20) {
		return;
	}
	SystemTestCase::copyFile("tridyn_system_PSI_3.dat");
	// 1D + HeDVI + advection + modifiedTM + attenuation + surface + reflective
	SystemTestCase{"system_PSI_3"}.run();
}

BOOST_AUTO_TEST_CASE(system_PSI_4)
{
	if (getMPICommSize() > 20) {
		return;
	}
	SystemTestCase::copyFile("tridyn_system_PSI_4.dat");
	// 1D + HeTVI + advection + modifiedTM + attenuation + surface + reflective
	SystemTestCase{"system_PSI_4"}.run();
}

BOOST_AUTO_TEST_CASE(system_PSI_5)
{
	if (getMPICommSize() > 20) {
		return;
	}
	// 1D + HeVI + advection + modifiedTM + attenuation + surface + reflective +
	// reduced matrix method bulk
	SystemTestCase{"system_PSI_5"}.run();
}

BOOST_AUTO_TEST_CASE(system_Fe_1)
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D
	SystemTestCase{"system_Fe_1", "bubble_36.dat"}.run();
}

BOOST_AUTO_TEST_CASE(system_Alloy_1)
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D
	SystemTestCase{"system_Alloy_1", "Alloy.dat"}.run();
}

BOOST_AUTO_TEST_SUITE_END()
