#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;

#include <xolotl/test/KokkosFixture.h>

BOOST_GLOBAL_FIXTURE(KokkosFixture);

#include <xolotl/test/MPITestUtils.h>
#include <xolotl/test/SystemTestCase.h>
using xolotl::test::getMPICommSize;
using xolotl::test::SystemTestCase;

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_SUITE(System)

BOOST_AUTO_TEST_CASE_WITH_DECOR(AZr_0, * utf::label("0D"))
{
	if (getMPICommSize() > 1) {
		return;
	}
	SystemTestCase{"system_AZr_0", "AlphaZr.dat"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_0, * utf::label("0D"))
{
	if (getMPICommSize() > 1) {
		return;
	}
	SystemTestCase{"system_NE_0"}.tolerance(1.0e-9).run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_1, * utf::label("0D"))
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D + heterogeneous nucleation
	SystemTestCase{"system_NE_1"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_2, * utf::label("0D"))
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D + grouping + re-solution
	SystemTestCase{"system_NE_2"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_3, * utf::label("1D"))
{
	if (getMPICommSize() > 25) {
		return;
	}
	// 1D + temperature gradient
	SystemTestCase{"system_NE_3"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_4, * utf::label("2D"))
{
	if (getMPICommSize() > 100) {
		return;
	}
	// 2D
	SystemTestCase{"system_NE_4"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_1, * utf::label("1D"))
{
	if (getMPICommSize() > 20) {
		return;
	}
	// 1D + HeVI + advection + modifiedTM + attenuation + surface + reflective
	// bulk
	SystemTestCase{"system_PSI_1"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_2, * utf::label("1D"))
{
	if (getMPICommSize() > 20) {
		return;
	}
	SystemTestCase::copyFile("tridyn_system_PSI_2.dat");
	// 1D + HeTVI + advection + modifiedTM + attenuation + surface + reflective
	SystemTestCase{"system_PSI_2"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_3, * utf::label("1D"))
{
	if (getMPICommSize() > 10) {
		return;
	}
	// 1D + HeVI + pulsed flux + sink + I grouping + surface + reflective
	// bulk
	SystemTestCase{"system_PSI_3"}.tolerance(1.0e-7).run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_4, * utf::label("1D"))
{
	if (getMPICommSize() > 100) {
		return;
	}
	// 1D + He + heat
	SystemTestCase{"system_PSI_4"}.tolerance(1.0e-4).run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(Fe_1, * utf::label("0D"))
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D
	SystemTestCase{"system_Fe_1", "bubble_36.dat"}.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(Alloy_1, * utf::label("0D"))
{
	if (getMPICommSize() > 1) {
		return;
	}
	// 0D
	SystemTestCase{"system_Alloy_1", "Alloy.dat"}.run();
}

BOOST_AUTO_TEST_SUITE_END()
