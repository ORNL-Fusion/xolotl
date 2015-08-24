#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <W100FitFluxHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W100FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkGetIncidentFlux) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

	// Create the flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Set the factor to change the helium flux
    testFitFlux->setFluxAmplitude(1.0);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(network, nGridpts, step);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.476819, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.225961, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.097220, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkFluxIndex) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

	// Create the flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Set the factor to change the helium flux
    testFitFlux->setFluxAmplitude(1.0);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(network, nGridpts, step);

    // Check the value of the index of the cluster for the flux
    BOOST_REQUIRE_EQUAL(testFitFlux->getIncidentFluxClusterIndex(), 0);

	return;
}

BOOST_AUTO_TEST_CASE(checkFluence) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

	// Create the flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Set the factor to change the helium flux
    testFitFlux->setFluxAmplitude(1.0);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(network, nGridpts, step);

	// Check that the fluence is 0.0 at the beginning
	BOOST_REQUIRE_EQUAL(testFitFlux->getFluence(), 0.0);

	// Increment the fluence
	testFitFlux->incrementFluence(1.0e-8);
	
	// Check that the fluence is not 0.0 anymore
	BOOST_REQUIRE_EQUAL(testFitFlux->getFluence(), 1.0e-8);

	return;
}

BOOST_AUTO_TEST_CASE(checkFluxAmplitude) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

	// Create the flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Set the factor to change the helium flux
    testFitFlux->setFluxAmplitude(1.0);
    // Set the factor to change the helium flux
    testFitFlux->setFluxAmplitude(2.5);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(network, nGridpts, step);

    // Check the value of the helium flux
    BOOST_REQUIRE_EQUAL(testFitFlux->getFluxAmplitude(), 2.5);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 1.192047, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.564902, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.243050, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkTimeProfileFlux) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 1.25;

	// Create a file with a time profile for the flux
	// First column with the time and the second with
	// the amplitude (in He/nm2/s) at that time.
	std::ofstream writeFluxFile("fluxFile.dat");
	writeFluxFile << "0.0 1000.0 \n"
			"1.0 4000.0 \n"
			"2.0 2000.0 \n"
			"3.0 3000.0 \n"
			"4.0 0.0";
	writeFluxFile.close();

	// Create the flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the time profile for the flux handler
    testFitFlux->initializeTimeProfile("fluxFile.dat");
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(network, nGridpts, step);

	// Create a time
	double currTime = 0.5;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 1192.047, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 564.902, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 243.050, 0.01);
	// Check the value of the helium flux
    BOOST_REQUIRE_EQUAL(testFitFlux->getFluxAmplitude(), 2500.0);

    // Change the current time
    currTime = 3.5;

	// Get the flux vector
	testFluxVec = testFitFlux->getIncidentFluxVec(currTime);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 715.228, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 338.941, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 145.830, 0.01);
	// Check the value of the helium flux
    BOOST_REQUIRE_EQUAL(testFitFlux->getFluxAmplitude(), 1500.0);

    // Remove the created file
    std::string tempFile = "fluxFile.dat";
    std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
