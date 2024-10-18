#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/diffusion/DummyDiffusionHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/ConfOptions.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>
#include <xolotl/util/MPIUtils.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace diffusion;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the DummyDiffusionHandler.
 */
BOOST_AUTO_TEST_SUITE(DummyDiffusionHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal part of the Jacobian,
 * and the compute diffusion methods for the dummy handler.
 */
BOOST_AUTO_TEST_CASE(checkDiffusion)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 1 0" << std::endl;
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
		temperatures.push_back(1000.0);
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
	DummyDiffusionHandler diffusionHandler(opts.getMigrationThreshold());

	// Initialize it
	std::vector<core::RowColPair> idPairs;
	diffusionHandler.initialize(network, idPairs);

	// Check the total number of diffusing clusters, here 0
	BOOST_REQUIRE_EQUAL(diffusionHandler.getNumberOfDiffusing(), 0);

	BOOST_REQUIRE(idPairs.empty());

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

	// Set the temperature to 1000 K to initialize the diffusion coefficients
	network.setTemperatures(temperatures, grid);

	// Get the offset for the grid point in the middle
	// Supposing the 3 grid points are laid-out as follow:
	// 0 | 1 | 2
	auto concOffset = subview(concentration, 1, Kokkos::ALL);
	auto updatedConcOffset = subview(newConcentration, 1, Kokkos::ALL);

	// Fill the concVector with the pointer to the middle, left, and right grid
	// points
	// Kokkos::Array<Kokkos::View<const double*>, 3> concVector;
	Kokkos::Array<Kokkos::View<const double*>, 3> concVector;
	concVector[0] = concOffset; // middle
	concVector[1] = subview(concentration, 0, Kokkos::ALL); // left
	concVector[2] = subview(concentration, 2, Kokkos::ALL); // right

	// Compute the diffusion at this grid point
	diffusionHandler.computeDiffusion(network,
		StencilConcArray{concVector.data(), 3}, updatedConcOffset, hx, hx, 1,
		1);

	// Check the new values of updatedConcOffset
	auto updatedConcOffsetMirror =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[0], 0.0, 0.01); // Does not diffuse
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[1], 0.0, 0.01); // Does not diffuse
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[3], 0.0, 0.01); // Does not diffuse
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[5], 0.0, 0.01); // Does not diffuse
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[7], 0.0, 0.01); // Does not diffuse

	// Don't even test the Jacobian because the number of diffusing cluster is 0

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
