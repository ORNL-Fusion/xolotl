#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/modified/SoretDiffusionHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/ConfOptions.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>
#include <xolotl/util/MPIUtils.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace modified;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the SoretDiffusionHandler.
 */
BOOST_AUTO_TEST_SUITE(SoretDiffusionHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal part of the Jacobian,
 * and the compute diffusion methods.
 */
BOOST_AUTO_TEST_CASE(checkDiffusion)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=1 0 0 0 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	util::mpiInit(cl.argc, cl.argv);
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	// Create a grid
	std::vector<double> grid;
	std::vector<double> temperatures;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double)l);
		temperatures.push_back(2000.0 - 200.0 * l);
	}

	// Create the network
	using NetworkType =
		network::PSIReactionNetwork<network::PSIFullSpeciesList>;
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe = opts.getMaxImpurity();
	NetworkType::AmountType maxD = opts.getMaxD();
	NetworkType::AmountType maxT = opts.getMaxT();
	NetworkType network({maxHe, maxD, maxT, maxV, maxI}, grid.size(), opts);
	// Get its size
	const int dof = network.getDOF();

	// Create the diffusion handler
	SoretDiffusionHandler soretHandler;

	// Create coordinate entries list
	std::vector<RowColPair> soretEntries;

	// Initialize it
	soretHandler.initialize(network, soretEntries, grid, 0);

	// Test which cluster diffuses
	BOOST_REQUIRE_EQUAL(soretEntries.size(), 1);
	BOOST_REQUIRE_EQUAL(soretEntries[0][0], 0);
	BOOST_REQUIRE_EQUAL(soretEntries[0][1], 0);

	// The size parameter in the x direction
	double hx = 1.0;

	// The arrays of concentration
	test::DOFView concentration("concentration", 3, dof);
	test::DOFView newConcentration("newConcentration", 3, dof);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {3, dof}),
		KOKKOS_LAMBDA(int i, int n) {
			auto id = static_cast<double>(i * dof + n);
			concentration(i, n) = id * id;
		});

	// Set the temperatures
	network.setTemperatures(temperatures, grid);

	// Get the offset for the grid point in the middle
	// Supposing the 3 grid points are laid-out as follow:
	// 0 | 1 | 2
	auto concOffset = subview(concentration, 1, Kokkos::ALL);
	auto updatedConcOffset = subview(newConcentration, 1, Kokkos::ALL);

	// Fill the concVector with the pointer to the middle, left, and right grid
	// points
	using ConcSubView = Kokkos::View<const double*>;
	Kokkos::Array<ConcSubView, 3> concVector;
	concVector[0] = concOffset; // middle
	concVector[1] = subview(concentration, 0, Kokkos::ALL); // left
	concVector[2] = subview(concentration, 2, Kokkos::ALL); // right

	// Compute the diffusion at this grid point
	soretHandler.computeDiffusion(network,
		core::StencilConcArray{concVector.data(), 3}, updatedConcOffset, hx, hx,
		1);

	// Check the new values of updatedConcOffset
	auto hUpdatedConcOffset =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(hUpdatedConcOffset[0], 27632823604, 0.01);

	// Initialize the indices and values to set in the Jacobian
	int nDiff = 1;
	Kokkos::View<double*> values("values", 3 * nDiff);

	// Compute the partial derivatives for the diffusion a the grid point 1
	auto valid = soretHandler.computePartialsForDiffusion(network,
		core::StencilConcArray{concVector.data(), 3}, values, hx, hx, 0);

	auto hVals = create_mirror_view_and_copy(Kokkos::HostSpace{}, values);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(valid, true);

	// Check some values
	BOOST_REQUIRE_CLOSE(hVals[0], -1523706589, 0.01);
	BOOST_REQUIRE_CLOSE(hVals[1], -8153185774, 0.01);
	BOOST_REQUIRE_CLOSE(hVals[2], 8153185774, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
