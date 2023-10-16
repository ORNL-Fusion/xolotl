#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/diffusion/Diffusion3DHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/Options.h>
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
 * This suite is responsible for testing the Diffusion3DHandler.
 */
BOOST_AUTO_TEST_SUITE(Diffusion3DHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal part of the Jacobian,
 * and the compute diffusion methods.
 */
BOOST_AUTO_TEST_CASE(checkDiffusion)
{
	// Create the option to create a network
	xolotl::options::Options opts;
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
	Diffusion3DHandler diffusionHandler(opts.getMigrationThreshold());

	// Create a collection of advection handlers
	std::vector<advection::IAdvectionHandler*> advectionHandlers;

	// Initialize it
	std::vector<core::RowColPair> idPairs;
	diffusionHandler.initialize(network, idPairs);
	diffusionHandler.initializeDiffusionGrid(
		advectionHandlers, grid, 5, 0, 3, 1.0, 0, 3, 1.0, 0);

	// Check the total number of diffusing clusters
	BOOST_REQUIRE_EQUAL(diffusionHandler.getNumberOfDiffusing(), 8);

	// The step size in the x direction
	double hx = 1.0;
	// The size parameter in the y direction
	double sy = 1.0;
	// The size parameter in the z direction
	double sz = 1.0;

	// The arrays of concentration
	test::DOFView concentration("concentration", 27, dof);
	test::DOFView newConcentration("newConcentration", 27, dof);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {27, dof}),
		KOKKOS_LAMBDA(int i, int n) {
			auto id = static_cast<double>(i * dof + n);
			concentration(i, n) = id * id / 10.0;
		});

	// Set the temperature to 1000K to initialize the diffusion coefficients
	network.setTemperatures(temperatures, grid);

	// Get the offset for the grid point in the middle
	// Supposing the 27 grid points are laid-out as follow (a cube!):
	// 6 | 7 | 8    15 | 16 | 17    24 | 25 | 26
	// 3 | 4 | 5    12 | 13 | 14    21 | 22 | 23
	// 0 | 1 | 2    9  | 10 | 11    18 | 19 | 20
	//   front         middle           back
	auto concOffset = subview(concentration, 13, Kokkos::ALL);
	auto updatedConcOffset = subview(newConcentration, 13, Kokkos::ALL);

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// top, front, and back grid points
	Kokkos::Array<Kokkos::View<const double*>, 7> concVector;
	concVector[0] = concOffset; // middle
	concVector[1] = subview(concentration, 12, Kokkos::ALL); // left
	concVector[2] = subview(concentration, 14, Kokkos::ALL); // right
	concVector[3] = subview(concentration, 10, Kokkos::ALL); // bottom
	concVector[4] = subview(concentration, 16, Kokkos::ALL); // top
	concVector[5] = subview(concentration, 4, Kokkos::ALL); // front
	concVector[6] = subview(concentration, 22, Kokkos::ALL); // back

	// Compute the diffusion at this grid point
	diffusionHandler.computeDiffusion(network,
		StencilConcArray{concVector.data(), 7}, updatedConcOffset, hx, hx, 0,
		sy, 1, sz, 1);

	// Check the new values of updatedConcOffset
	auto updatedConcOffsetMirror =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[1], 3.3744e+13, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[3], 1.6526e+13, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[5], 6.6489e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[7], 8.7794e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[9], 6.5338e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[11], 1.6182e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[13], 2.5353e+10, 0.01);
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[15], 0.0, 0.01); // He_8 does not diffuse
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[0], 2.6578e+09, 0.01);

	// Initialize the indices and values to set in the Jacobian
	int nDiff = diffusionHandler.getNumberOfDiffusing();
	auto val = Kokkos::View<double*>("val", 7 * nDiff);

	// Compute the partial derivatives for the diffusion a the grid point 1
	diffusionHandler.computePartialsForDiffusion(
		network, val, hx, hx, 0, sy, 1, sz, 1);

	// Check some values
	auto valMirror = create_mirror_view_and_copy(Kokkos::HostSpace{}, val);
	BOOST_REQUIRE_CLOSE(valMirror[0], -3031876, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[5], 505312, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[12], 6415444736, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[20], 3141913616, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[26], 1264105042, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
