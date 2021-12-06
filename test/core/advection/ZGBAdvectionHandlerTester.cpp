#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/advection/ZGBAdvectionHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/Options.h>
#include <xolotl/test/CommandLine.h>
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
	network.syncClusterDataOnHost();
	network.getSubpaving().syncZones(plsm::onHost);
	// Get its size
	const int dof = network.getDOF();

	// Create ofill
	network::IReactionNetwork::SparseFillMap ofill;

	// Create the advection handler and initialize it with a sink at
	// 2nm in the Z direction
	ZGBAdvectionHandler advectionHandler;
	advectionHandler.initialize(network, ofill);
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
	double concentration[27 * dof];
	double newConcentration[27 * dof];

	// Initialize their values
	for (int i = 0; i < 27 * dof; i++) {
		concentration[i] = (double)i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	network.setTemperatures(temperatures, grid);
	network.syncClusterDataOnHost();

	// Get pointers
	double* conc = &concentration[0];
	double* updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 27 grid points are laid-out as follow (a cube!):
	// 6 | 7 | 8    15 | 16 | 17    24 | 25 | 26
	// 3 | 4 | 5    12 | 13 | 14    21 | 22 | 23
	// 0 | 1 | 2    9  | 10 | 11    18 | 19 | 20
	//   front         middle           back
	double* concOffset = conc + 13 * dof;
	double* updatedConcOffset = updatedConc + 13 * dof;

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// top, front, and back grid points
	double* concVector[7]{};
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 12 * dof; // left
	concVector[2] = conc + 14 * dof; // right
	concVector[3] = conc + 10 * dof; // bottom
	concVector[4] = conc + 16 * dof; // top
	concVector[5] = conc + 4 * dof; // front
	concVector[6] = conc + 22 * dof; // back

	// Set the grid position
	plsm::SpaceVector<double, 3> gridPosition{hx, hy, hz};

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(network, gridPosition, concVector,
		updatedConcOffset, hx, hx, 0, hy, 1, hz, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], 5.4795e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], 5.0808e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], 6.2076e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 1.0761e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], 1.4932e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[11], 5.5940e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[13], 1.5789e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[15], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	IdType indices[nAdvec];
	double val[7 * nAdvec];
	// Get the pointer on them for the compute advection method
	IdType* indicesPointer = &indices[0];
	double* valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(network, valPointer,
		indicesPointer, gridPosition, hx, hx, 0, hy, 1, hz, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 1);
	BOOST_REQUIRE_EQUAL(indices[1], 3);
	BOOST_REQUIRE_EQUAL(indices[2], 5);
	BOOST_REQUIRE_EQUAL(indices[3], 7);
	BOOST_REQUIRE_EQUAL(indices[4], 9);
	BOOST_REQUIRE_EQUAL(indices[5], 11);
	BOOST_REQUIRE_EQUAL(indices[6], 13);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], 3.76893e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 3.76893e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 3.45234e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 3.45234e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 4.16701e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 4.16701e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[6], 7.13648e+06, 0.01);

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
