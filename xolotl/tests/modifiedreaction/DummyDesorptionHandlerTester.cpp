#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <DummyDesorptionHandler.h>
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
 * This suite is responsible for testing the DummyDesorptionHandler.
 */
BOOST_AUTO_TEST_SUITE(DummyDesorptionHandler_testSuite)

/**
 * Method checking the initialization and the compute desorption methods.
 */
BOOST_AUTO_TEST_CASE(checkDesorption) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 1 0 2 6" << std::endl;
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
	network->addGridPoints(1);
	network->setTemperature(1000.0, 0);

	// Suppose we have a grid with 13 grip points and distance of
	// 0.1 nm between grid points
	std::vector<double> grid;
	for (int l = 0; l < 13; l++) {
		grid.push_back((double) l * 0.1);
	}
	// Set the surface position
	int surfacePos = 0;

	// Create the desorption handler
	DummyDesorptionHandler desorptionHandler;

	// Initialize it
	desorptionHandler.initialize(grid);
	desorptionHandler.initializeIndex1D(surfacePos, *network, grid);

	// The arrays of concentration
	double concentration[13 * dof];
	double newConcentration[13 * dof];

	// Initialize their values
	for (int i = 0; i < 13 * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the second grid point
	double *concOffset = conc + dof;
	double *updatedConcOffset = updatedConc + dof;

	// Compute the desorption at the second grid point
	desorptionHandler.computeDesorption(concOffset, updatedConcOffset, 1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[14], 0.0, 0.01); // D_1

	// Initialize the indices and values to set in the Jacobian
	int nD = desorptionHandler.getNumberOfDesorbing();
	int indices[nD];
	double val[nD];
	// Get the pointer on them for the compute desorption method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the desorption at the grid point 1
	int nDesorbing = desorptionHandler.computePartialsForDesorption(valPointer,
			indicesPointer, 1, 0);

	// Verify that no cluster is undergoing desorption
	BOOST_REQUIRE_EQUAL(nDesorbing, 0);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
