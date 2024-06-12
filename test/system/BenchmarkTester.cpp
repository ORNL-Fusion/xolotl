#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test_suite.hpp>
namespace utf = boost::unit_test;

#include <xolotl/test/KokkosFixture.h>

BOOST_GLOBAL_FIXTURE(KokkosFixture);

#include <xolotl/test/MPITestUtils.h>
#include <xolotl/test/SystemTestCase.h>
using xolotl::test::getMPICommSize;
using xolotl::test::SystemTestCase;

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_SUITE(Benchmark)

BOOST_AUTO_TEST_CASE_WITH_DECOR(AZr_1, *utf::label("0D"))
{
	// 0D, 1000 in each direction, grouped
	SystemTestCase{"benchmark_AZr_1", "AlphaZr.dat"}
		.mpiLimits(1, 1)
		.withTimer()
		.run();
}

// BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_1, *utf::label("0D"))
//{
//	// 0D, 10000 DOF ungrouped
//	SystemTestCase{"benchmark_NE_1"}.mpiLimits(1, 1).withTimer().run();
// }
//
// BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_2, *utf::label("0D"))
//{
//	// 0D, 2010 DOF grouped
//	SystemTestCase{"benchmark_NE_2"}.mpiLimits(1, 1).withTimer().run();
// }
//
// BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_3, *utf::label("3D"))
//{
//	// 3D
//	SystemTestCase{"benchmark_NE_3"}.mpiLimits(1, 3969).withTimer().run();
// }
//
// BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_4, *utf::label("2D"))
//{
//	// 2D, longer
//	SystemTestCase{"benchmark_NE_4"}.mpiLimits(16, 15876).withTimer().run();
// }
//
// BOOST_AUTO_TEST_CASE_WITH_DECOR(NE_5, *utf::label("0D"))
//{
//	// 0D, grouped, re-solution
//	SystemTestCase{"benchmark_NE_5"}
//		.mpiLimits(1, 1)
//		.tolerance(5.0e-9)
//		.withTimer()
//		.run();
// }

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_1, *utf::label("1D"))
{
	// 1D + HeV + 4e25 flux W100
	SystemTestCase{"benchmark_PSI_1"}
		.tolerance(1.0e-5)
		.mpiLimits(1, 100)
		.withTimer()
		.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_2, *utf::label("1D"))
{
	// 1D + + HeV + 5e27 flux W100 + bursting
	SystemTestCase{"benchmark_PSI_2"}.mpiLimits(1, 1).withTimer().run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_3, *utf::label("1D"))
{
	// 1D + HeV + case g + bursting
	SystemTestCase{"benchmark_PSI_3"}.mpiLimits(32, 32).withTimer().run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_4, *utf::label("1D"))
{
	SystemTestCase::copyFile("tridyn_benchmark_PSI_4.dat");
	// 1D + ITER_He
	SystemTestCase{"benchmark_PSI_4"}
		.mpiLimits(32, 32)
		.tolerance(5.0e-10)
		.withTimer()
		.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_5, *utf::label("1D"))
{
	SystemTestCase::copyFile("tridyn_benchmark_PSI_5.dat");
	// 1D + ITER_BPO
	SystemTestCase{"benchmark_PSI_5"}
		.mpiLimits(32, 32)
		.tolerance(5.0e-10)
		.withTimer()
		.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_7, *utf::label("1D"))
{
	// 1D + pulsed
	SystemTestCase{"benchmark_PSI_7"}
		.tolerance(1.0e-8)
		.mpiLimits(4, 25)
		.withTimer()
		.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_8, *utf::label("1D"))
{
	SystemTestCase::copyFile("flux_benchmark_PSI_8.dat");
	SystemTestCase::copyFile("tridyn_benchmark_PSI_8.dat");
	// 1D + PISCES + varying flux
	SystemTestCase{"benchmark_PSI_8"}
		.tolerance(1.0e-6)
		.mpiLimits(1, 20)
		.withTimer()
		.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_9, *utf::label("1D"))
{
	SystemTestCase::copyFile("temp_benchmark_PSI_9.dat");
	SystemTestCase::copyFile("tridyn_benchmark_PSI_9.dat");
	// 1D + PISCES + varying temperature
	SystemTestCase{"benchmark_PSI_9"}
		.tolerance(1.0e-6)
		.mpiLimits(1, 20)
		.withTimer()
		.run();
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(PSI_10, *utf::label("1D"))
{
	// 1D + reduced jacobian
	SystemTestCase{"benchmark_PSI_10"}.mpiLimits(1, 20).withTimer().run();
}

BOOST_AUTO_TEST_SUITE_END()
