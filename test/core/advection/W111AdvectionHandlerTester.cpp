#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/advection/W111AdvectionHandler.h>
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
 * This suite is responsible for testing the W111AdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(W111AdvectionHandler_testSuite)

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

	// Create a collection of advection handlers
	std::vector<IAdvectionHandler*> advectionHandlers;

	// Create the advection handler and initialize it
	W111AdvectionHandler advectionHandler;
	advectionHandler.initialize(network, ofill);
	advectionHandler.initializeAdvectionGrid(advectionHandlers, grid, 3, 0);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 7);

	// Check the clusters in ofill
	BOOST_REQUIRE_EQUAL(ofill[1][0], 1);
	BOOST_REQUIRE_EQUAL(ofill[3][0], 3);
	BOOST_REQUIRE_EQUAL(ofill[5][0], 5);
	BOOST_REQUIRE_EQUAL(ofill[7][0], 7);
	BOOST_REQUIRE_EQUAL(ofill[9][0], 9);
	BOOST_REQUIRE_EQUAL(ofill[11][0], 11);
	BOOST_REQUIRE_EQUAL(ofill[13][0], 13);

	// Set the size parameter in the x direction
	double hx = 1.0;

	// Create the arrays of concentration
	double concentration[3 * dof];
	double newConcentration[3 * dof];

	// Initialize their values
	for (int i = 0; i < 3 * dof; i++) {
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
	double* concOffset = conc + dof;
	double* updatedConcOffset = updatedConc + dof;

	// Fill the concVector with the pointer to the middle, left, and right grid
	// points
	double* concVector[3]{};
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * dof; // right

	// Set the grid position
	plsm::SpaceVector<double, 3> gridPosition{hx, 0.0, 0.0};

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(
		network, gridPosition, concVector, updatedConcOffset, hx, hx, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -2.0171e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], -2.2012e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], -2.8037e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], -2.6928e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[9], -1.0761e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[11], -3.6701e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[13], -1.0429e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[15], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	IdType indices[nAdvec];
	double val[2 * nAdvec];
	// Get the pointer on them for the compute advection method
	IdType* indicesPointer = &indices[0];
	double* valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(
		network, valPointer, indicesPointer, gridPosition, hx, hx, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 1);
	BOOST_REQUIRE_EQUAL(indices[1], 3);
	BOOST_REQUIRE_EQUAL(indices[2], 5);
	BOOST_REQUIRE_EQUAL(indices[3], 7);
	BOOST_REQUIRE_EQUAL(indices[4], 9);
	BOOST_REQUIRE_EQUAL(indices[5], 11);
	BOOST_REQUIRE_EQUAL(indices[6], 13);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -815207266.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 50950454.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -700039625.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 43752477.0, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 1); // x
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 0);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
