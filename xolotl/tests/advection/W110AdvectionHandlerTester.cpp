#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <W110AdvectionHandler.h>
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
 * This suite is responsible for testing the W110AdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(W110AdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 1 0" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array
	opts.readParams(argv);

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

	// Create ofill
	xolotlCore::IReactionNetwork::SparseFillMap ofill;

	// Create a collection of advection handlers
	std::vector<IAdvectionHandler *> advectionHandlers;

	// Create the advection handler and initialize it
	W110AdvectionHandler advectionHandler;
	advectionHandler.initialize(*network, ofill);
	advectionHandler.initializeAdvectionGrid(advectionHandlers, grid);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 7);

	// Check the clusters in ofill
	BOOST_REQUIRE_EQUAL(ofill[0][0], 0);
	BOOST_REQUIRE_EQUAL(ofill[1][0], 1);
	BOOST_REQUIRE_EQUAL(ofill[2][0], 2);
	BOOST_REQUIRE_EQUAL(ofill[3][0], 3);
	BOOST_REQUIRE_EQUAL(ofill[4][0], 4);
	BOOST_REQUIRE_EQUAL(ofill[5][0], 5);
	BOOST_REQUIRE_EQUAL(ofill[6][0], 6);

	// Set the size parameter in the x direction
	double hx = 1.0;

	// Create the arrays of concentration
	double concentration[3 * dof];
	double newConcentration[3 * dof];

	// Initialize their values
	for (int i = 0; i < 3 * dof; i++) {
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
	double *concOffset = conc + dof;
	double *updatedConcOffset = updatedConc + dof;

	// Fill the concVector with the pointer to the middle, left, and right grid points
	double **concVector = new double*[3];
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * dof; // right

	// Set the grid position
	NDPoint<3> gridPosition { hx, 0.0, 0.0 };

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(*network, gridPosition, concVector,
			updatedConcOffset, hx, hx, 1, 1);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], -5.56329e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -4.93646e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], -1.00995e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], -1.36081e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], -6.08651e+11, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], -1.85328e+10, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], -3.54365e+09, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int indices[nAdvec];
	double val[2 * nAdvec];
	// Get the pointer on them for the compute advection method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(*network, valPointer,
			indicesPointer, gridPosition, hx, hx, 1, 1);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);
	BOOST_REQUIRE_EQUAL(indices[4], 4);
	BOOST_REQUIRE_EQUAL(indices[5], 5);
	BOOST_REQUIRE_EQUAL(indices[6], 6);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -205476899.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 12842306.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -161884163.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 10117760.0, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 1); //x
	BOOST_REQUIRE_EQUAL(stencil[1], 0);
	BOOST_REQUIRE_EQUAL(stencil[2], 0);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
