#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <Diffusion2DHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <Options.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the Diffusion2DHandler.
 */
BOOST_AUTO_TEST_SUITE(Diffusion2DHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal part of the Jacobian,
 * and the compute diffusion methods.
 */
BOOST_AUTO_TEST_CASE(checkDiffusion) {

	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 1 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	// Initialize MPI for HDF5
	MPI_Init(&argc, &argv);
	opts.readParams(argc, argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Create the network
	auto network = loader.generate(opts);
	// Get its size
	const int dof = network->getDOF();
	// Initialize the rates
	network->addGridPoints(3);

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l);
	}

	// Create the diffusion handler
	Diffusion2DHandler diffusionHandler;

	// Create a collection of advection handlers
	std::vector<IAdvectionHandler *> advectionHandlers;

	// Create ofill
	xolotlCore::IReactionNetwork::SparseFillMap ofill;

	// Initialize it
	diffusionHandler.initializeOFill(*network, ofill);
	diffusionHandler.initializeDiffusionGrid(advectionHandlers, grid, 5, 0, 3,
			1.0, 0);

	// Check the total number of diffusing clusters
	BOOST_REQUIRE_EQUAL(diffusionHandler.getNumberOfDiffusing(), 8);

	// The step size in the x direction
	double hx = 1.0;
	// The size parameter in the y direction
	double sy = 1.0;

	// The arrays of concentration
	double concentration[9 * dof];
	double newConcentration[9 * dof];

	// Initialize their values
	for (int i = 0; i < 9 * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	network->setTemperature(1000.0, 0);
	network->setTemperature(1000.0, 1);
	network->setTemperature(1000.0, 2);

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 9 grid points are laid-out as follow:
	// 6 | 7 | 8
	// 3 | 4 | 5
	// 0 | 1 | 2
	double *concOffset = conc + 4 * dof;
	double *updatedConcOffset = updatedConc + 4 * dof;

	// Fill the concVector with the pointer to the middle, left, right, bottom, and top grid points
	double **concVector = new double*[5];
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 3 * dof; // left
	concVector[2] = conc + 5 * dof; // right
	concVector[3] = conc + 1 * dof; // bottom
	concVector[4] = conc + 7 * dof; // top

	// Compute the diffusion at this grid point
	diffusionHandler.computeDiffusion(*network, concVector, updatedConcOffset,
			hx, hx, 0, sy, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 4.632e+13, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], 2.2685e+13, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], 9.1268e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], 1.2051e+13, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], 8.9687e+12, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], 2.2213e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], 3.4801e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // Does not diffuse
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 3.6484e+09, 0.01);

	// Initialize the indices and values to set in the Jacobian
	int nDiff = diffusionHandler.getNumberOfDiffusing();
	int indices[nDiff];
	double val[5 * nDiff];
	// Get the pointer on them for the compute diffusion method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the diffusion a the grid point 1
	diffusionHandler.computePartialsForDiffusion(*network, valPointer,
			indicesPointer, hx, hx, 0, sy, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);
	BOOST_REQUIRE_EQUAL(indices[4], 4);
	BOOST_REQUIRE_EQUAL(indices[5], 5);
	BOOST_REQUIRE_EQUAL(indices[6], 6);
	BOOST_REQUIRE_EQUAL(indices[7], 8);

	// Check some values
	BOOST_REQUIRE_CLOSE(val[0], -2.56618e+10, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 6.41544e+09, 0.01);
	BOOST_REQUIRE_CLOSE(val[8], 3.14191e+09, 0.01);
	BOOST_REQUIRE_CLOSE(val[12], 1.26411e+09, 0.01);
	BOOST_REQUIRE_CLOSE(val[15], -6.67657e+09, 0.01);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
