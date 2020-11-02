#pragma once

#include <mpi.h>

namespace xolotl
{
namespace test
{
inline int
getMPIRank()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

inline int
getMPICommSize()
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	return size;
}
} // namespace test
} // namespace xolotl

// Do not put this in a namespace - it causes macros like
// BOOST_GLOBAL_FIXTURE to break.

struct MPIFixture
{
	MPIFixture(void)
	{
		BOOST_TEST_MESSAGE("Initializing MPI");
		auto& mts = boost::unit_test::framework::master_test_suite();
		MPI_Init(&mts.argc, &mts.argv);
		BOOST_TEST_MESSAGE("Done initializing MPI");
	}

	~MPIFixture(void)
	{
		BOOST_TEST_MESSAGE("Finalizing MPI");
		MPI_Finalize();
		BOOST_TEST_MESSAGE("Done finalizing MPI");
	}
};
