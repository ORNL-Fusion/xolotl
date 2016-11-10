#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W211FitFluxHandler.h"
#include <mpi.h>
#include <HDF5NetworkLoader.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the W211FitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W211FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = loader.load().get();

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the W211 flux handler
	auto testFitFlux = make_shared<W211FitFluxHandler>();
	// Set the flux amplitude
	testFitFlux->setFluxAmplitude(1.0);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(network, surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.454047, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.249781, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.096172, 0.01);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
