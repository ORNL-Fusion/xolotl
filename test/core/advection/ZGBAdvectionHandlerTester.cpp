#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/advection/ZGBAdvectionHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/Util.h>
#include <xolotl/util/MPIUtils.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace advection;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the ZGBAdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(ZGBAdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection)
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

	// Create the advection handler and initialize it with a sink at
	// 2nm in the Z direction
	ZGBAdvectionHandler advectionHandler;
	std::vector<core::RowColPair> idPairs;
	advectionHandler.initialize(network, idPairs);
	advectionHandler.setLocation(2.0);
	advectionHandler.setDimension(3);

	// Check if grid points are on the sink
	plsm::SpaceVector<double, 3> pos0{0.1, 3.0, 0.0};
	plsm::SpaceVector<double, 3> pos1{2.0, 2.0, 2.0};
	BOOST_REQUIRE_EQUAL(advectionHandler.isPointOnSink(pos0), false);
	BOOST_REQUIRE_EQUAL(advectionHandler.isPointOnSink(pos1), true);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 7);

	// Set the size parameters
	double hx = 1.0;
	double hy = 0.5;
	double hz = 2.0;

	// The arrays of concentration
	test::DOFView concentration("concentration", 27, dof);
	test::DOFView newConcentration("newConcentration", 27, dof);

	// Initialize their values
	Kokkos::parallel_for(
		Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {27, dof}),
		KOKKOS_LAMBDA(int i, int n) {
			auto id = static_cast<double>(i * dof + n);
			concentration(i, n) = id * id;
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

	// Set the grid position
	plsm::SpaceVector<double, 3> gridPosition{hx, hy, hz};

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(network, gridPosition,
		StencilConcArray{concVector.data(), 7}, updatedConcOffset, hx, hx, 0,
		hy, 1, hz, 1);

	// Check the new values of updatedConcOffset
	auto updatedConcOffsetMirror =
		create_mirror_view_and_copy(Kokkos::HostSpace{}, updatedConcOffset);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[1], 5.4795e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[3], 5.0808e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[5], 6.2076e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[7], 1.0761e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[9], 1.4932e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[11], 5.5940e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffsetMirror[13], 1.5789e+10, 0.01);
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[0], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(
		updatedConcOffsetMirror[15], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	auto val = Kokkos::View<double*>("val", 7 * nAdvec);

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(
		network, val, gridPosition, hx, hx, 0, hy, 1, hz, 1);

	// Check values
	auto valMirror = create_mirror_view_and_copy(Kokkos::HostSpace{}, val);
	BOOST_REQUIRE_CLOSE(valMirror[0], 3.76893e+06, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[1], 3.76893e+06, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[2], 3.45234e+06, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[3], 3.45234e+06, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[4], 4.16701e+06, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[5], 4.16701e+06, 0.01);
	BOOST_REQUIRE_CLOSE(valMirror[6], 7.13648e+06, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 0);
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 1); // z

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
