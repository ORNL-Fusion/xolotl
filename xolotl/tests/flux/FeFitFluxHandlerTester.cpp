#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "FeFitFluxHandler.h"
#include <mpi.h>
#include <HDF5NetworkLoader.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the FeFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FeFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = loader.load().get();
	// Get its size
	const int dof = network->getDOF();

	// Create an empty grid because we want 0D
	std::vector<double> grid;
	// Specify the surface position
	int surfacePos = 0;

	// Create the iron flux handler
	auto testFitFlux = make_shared<FeFitFluxHandler>();
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// The array of concentration
	double newConcentration[dof];

	// Initialize their values
	for (int i = 0; i < dof; i++) {
		newConcentration[i] = 0.0;
	}

	// The pointer to the grid point we want
	double *updatedConc = &newConcentration[0];
	double *updatedConcOffset = updatedConc;

	// Update the concentrations
	testFitFlux->computeIncidentFlux(currTime, updatedConcOffset, 0, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(newConcentration[0], 1.49e-05, 0.01); // I
	BOOST_REQUIRE_CLOSE(newConcentration[1], 0.0, 0.01); // I_2
	BOOST_REQUIRE_CLOSE(newConcentration[6], 2.11e-11, 0.01); // He
	BOOST_REQUIRE_CLOSE(newConcentration[14], 9.91e-06, 0.01); // V
	BOOST_REQUIRE_CLOSE(newConcentration[24], 1.51e-06, 0.01); // V_2
	BOOST_REQUIRE_CLOSE(newConcentration[39], 2.60e-07, 0.01); // V_3
	BOOST_REQUIRE_CLOSE(newConcentration[58], 1.58e-07, 0.01); // V_4
	BOOST_REQUIRE_CLOSE(newConcentration[79], 6.29e-08, 0.01); // V_5
	BOOST_REQUIRE_CLOSE(newConcentration[215], 3.16e-08, 0.01); // V_9

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
