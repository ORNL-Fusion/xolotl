#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <DesorptionHandler.h>
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
 * This suite is responsible for testing the DesorptionHandler.
 */
BOOST_AUTO_TEST_SUITE(DesorptionHandler_testSuite)

/**
 * Method checking the initialization and the compute desorption methods.
 */
BOOST_AUTO_TEST_CASE(checkModifiedDesorption) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 1 0 2 6" << std::endl;
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

	// Suppose we have a grid with 13 grip points and distance of
	// 0.1 nm between grid points
	int nGrid = 13;
	// Initialize the rates
	network->addGridPoints(nGrid);
	std::vector<double> grid;
	for (int l = 0; l < nGrid; l++) {
		grid.push_back((double) l * 0.1);
		network->setTemperature(1200.0, l);
	}
	// Set the surface position
	int surfacePos = 0;

	// Create the desorption handler
	DesorptionHandler desorptionHandler;

	// Initialize it
	desorptionHandler.setSolutionEnergy(1.0);
	desorptionHandler.setTemperature(1200.0);
	desorptionHandler.initialize(13);
	desorptionHandler.initializeIndex1D(surfacePos, *network, 13, 0);

	// The arrays of concentration
	double concentration[nGrid * dof];
	double newConcentration[nGrid * dof];

	// Initialize their values
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the first grid point
	double *concOffset = conc + 1 * dof;
	double *updatedConcOffset = updatedConc + 1 * dof;

	// Putting the concentrations in the network so that the rate for
	// desorption is computed correctly
	network->updateConcentrationsFromArray(concOffset);

	// Compute the desorption at the first grid point
	desorptionHandler.computeDesorption(concOffset, updatedConcOffset, 1, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[14], -9.84040153e-14, 0.01); // D_1

	// Get the offset for the ninth grid point
	concOffset = conc + 8 * dof;
	updatedConcOffset = updatedConc + 8 * dof;

	// Putting the concentrations in the network so that the rate for
	// desorption is computed correctly
	network->updateConcentrationsFromArray(concOffset);

	// Compute the desorption at the ninth grid point
	desorptionHandler.computeDesorption(concOffset, updatedConcOffset, 8);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[14], 0.0, 0.01); // D_1

	// Initialize the indices and values to set in the Jacobian
	int nD = desorptionHandler.getNumberOfDesorbing();
	int indices[nD];
	double val[nD];
	// Get the pointer on them for the compute desorption method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the desorption at the grid point 8
	int nDesorbing = desorptionHandler.computePartialsForDesorption(valPointer,
			indicesPointer, 8);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(nDesorbing, 0);

	// Change the grid point
	concOffset = conc + 1 * dof;
	updatedConcOffset = updatedConc + 1 * dof;
	network->updateConcentrationsFromArray(concOffset);

	// Compute the partial derivatives for the desorption at the grid point 1
	nDesorbing = desorptionHandler.computePartialsForDesorption(valPointer,
			indicesPointer, 1, 0);

	// Check values
	BOOST_REQUIRE_EQUAL(nDesorbing, 1);
	BOOST_REQUIRE_EQUAL(indices[0], 14);
	BOOST_REQUIRE_CLOSE(val[0], -6.809966e-18, 0.01); // D_1

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
