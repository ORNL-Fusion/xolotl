#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <mpi.h>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/advection/XGBAdvectionHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/Options.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace advection;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

/**
 * This suite is responsible for testing the XGBAdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(XGBAdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection)
{
	// Create the option to create a network
	xolotl::options::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 1 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char** argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	// Initialize MPI
	MPI_Init(&argc, &argv);
	opts.readParams(argc, argv);

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
	// 2nm in the X direction
	XGBAdvectionHandler advectionHandler;
	advectionHandler.initialize(network, ofill);
	advectionHandler.setLocation(2.0);
	advectionHandler.setDimension(2);

	// Check if grid points are on the sink
	plsm::SpaceVector<double, 3> pos0{0.1, 3.0, 0.0};
	plsm::SpaceVector<double, 3> pos1{2.0, 2.0, 0.0};
	BOOST_REQUIRE_EQUAL(advectionHandler.isPointOnSink(pos0), false);
	BOOST_REQUIRE_EQUAL(advectionHandler.isPointOnSink(pos1), true);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 7);

	// Set the size parameters
	double hx = 1.0;
	double hy = 0.5;

	// The arrays of concentration
	double concentration[9 * dof];
	double newConcentration[9 * dof];

	// Initialize their values
	for (int i = 0; i < 9 * dof; i++) {
		concentration[i] = (double)i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	network.setTemperatures(temperatures);
	network.syncClusterDataOnHost();

	// Get pointers
	double* conc = &concentration[0];
	double* updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 9 grid points are laid-out as follow:
	// 6 | 7 | 8
	// 3 | 4 | 5
	// 0 | 1 | 2
	double* concOffset = conc + 4 * dof;
	double* updatedConcOffset = updatedConc + 4 * dof;

	// Fill the concVector with the pointer to the middle, left, right, bottom,
	// and top grid points
	double** concVector = new double*[5];
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 3 * dof; // left
	concVector[2] = conc + 5 * dof; // right
	concVector[3] = conc + 1 * dof; // bottom
	concVector[4] = conc + 7 * dof; // top

	// Set the grid position
	plsm::SpaceVector<double, 3> gridPosition{hx, hy, 0.0};

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(
		network, gridPosition, concVector, updatedConcOffset, hx, hx, 0, hy, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -5.5382e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], -5.3677e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], -6.8446e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], -1.2366e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], -1.7859e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[11], -6.9542e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[13], -2.0377e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[15], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int indices[nAdvec];
	double val[3 * nAdvec];
	// Get the pointer on them for the compute advection method
	int* indicesPointer = &indices[0];
	double* valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(
		network, valPointer, indicesPointer, gridPosition, hx, hx, 0, hy, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 1);
	BOOST_REQUIRE_EQUAL(indices[1], 3);
	BOOST_REQUIRE_EQUAL(indices[2], 5);
	BOOST_REQUIRE_EQUAL(indices[3], 7);
	BOOST_REQUIRE_EQUAL(indices[4], 9);
	BOOST_REQUIRE_EQUAL(indices[5], 11);
	BOOST_REQUIRE_EQUAL(indices[6], 13);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -1.20606e+08, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 7.53788e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -1.10475e+08, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 6.90469e+06, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], -1.33344e+08, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 8.33401e+06, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], -1); // x
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 0);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
