#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <PSIClusterReactionNetwork.h>
#include <DummyHandlerRegistry.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>
#include <Options.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HDF5Utils
 */
BOOST_AUTO_TEST_SUITE(HDF5Utils_testSuite)

/**
 * Method checking the writing and reading of the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkIO) {
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

	// Create the options needed to load the network
	Options opts;
	// Load the network
	auto network = loader.load(opts);

	// Get the size of the network
	int networkSize = network->size();
	// Set the time step number
	int timeStep = 0;
	// Initialize the HDF5 file
	HDF5Utils::initializeFile("test.h5");

	// Set the number of grid points and step size
	int nGrid = 5;
	double stepSize = 0.5;
	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;
	// Set the surface information
	int iSurface = 3;
	double nInter = 1.0, previousFlux = 0.1;
	// Write the header in the HDF5 file
	HDF5Utils::fillHeader(nGrid, stepSize);

	// Write the network in the HDF5 file
	HDF5Utils::fillNetwork(filename);

	// Finalize the HDF5 file
	HDF5Utils::finalizeFile();

	// Open it again to add the concentrations
	HDF5Utils::openFile("test.h5");

	// Add the concentration sub group
	HDF5Utils::addConcentrationSubGroup(timeStep, currentTime, previousTime,
			currentTimeStep);

	// Write the surface position
	HDF5Utils::writeSurface1D(timeStep, iSurface, nInter, previousFlux);

	// Add the concentration dataset
	int length = 5;
	int gridPoint = 0;
	HDF5Utils::addConcentrationDataset(length, gridPoint);

	// Create a vector of concentration for one grid point
	double concArray[length][2];
	// Fill it
	for (int i = 0; i < length; i++) {
		// Fill the concArray
		concArray[i][0] = (double) i;
		concArray[i][1] = (double) i * 10.0 - 5.0;
	}

	// Write the concentrations in the HDF5 file
	HDF5Utils::fillConcentrations(concArray, gridPoint);

	// Close the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	// Read the header of the written file
	int nx = 0, ny = 0, nz = 0;
	double hx = 0.0, hy = 0.0, hz = 0.0;
	HDF5Utils::readHeader("test.h5", nx, hx, ny, hy, nz, hz);
	// Check the obtained values
	BOOST_REQUIRE_EQUAL(nx, nGrid);
	BOOST_REQUIRE_CLOSE(hx, stepSize, 0.0001);
	BOOST_REQUIRE_EQUAL(ny, 0);
	BOOST_REQUIRE_CLOSE(hy, 0.0, 0.0001);
	BOOST_REQUIRE_EQUAL(nz, 0);
	BOOST_REQUIRE_CLOSE(hz, 0.0, 0.0001);

	// Read the times
	double t = 0.0, dt = 0.0;
	HDF5Utils::readTimes("test.h5", 0, t, dt);
	BOOST_REQUIRE_CLOSE(t, currentTime, 0.0001);
	BOOST_REQUIRE_CLOSE(dt, currentTimeStep, 0.0001);
	double previousReadTime = HDF5Utils::readPreviousTime("test.h5", 0);
	BOOST_REQUIRE_CLOSE(previousReadTime, previousTime, 0.0001);

	// Read the surface position
	int surfacePos = HDF5Utils::readSurface1D("test.h5", timeStep);
	BOOST_REQUIRE_EQUAL(surfacePos, iSurface);
	double nInterstitial = HDF5Utils::readNInterstitial1D("test.h5", timeStep);
	BOOST_REQUIRE_CLOSE(nInterstitial, nInter, 0.0001);
	double previousIFlux = HDF5Utils::readPreviousIFlux1D("test.h5", timeStep);
	BOOST_REQUIRE_CLOSE(previousIFlux, previousFlux, 0.0001);

	// Read the network of the written file
	auto networkVector = HDF5Utils::readNetwork("test.h5");
	// Get all the reactants
	auto const& reactants = network->getAll();
	// Check the network vector
	for (IReactant& it : reactants) {
		// Get the i-th reactant in the network
		auto& reactant = (PSICluster&) it;
		int id = reactant.getId() - 1;
		// Get the corresponding line from the HDF5 file
		auto line = networkVector.at(id);

		// Check the composition
		auto& composition = reactant.getComposition();
		BOOST_REQUIRE_EQUAL((int ) line[0],
				composition[toCompIdx(Species::He)]);
		BOOST_REQUIRE_EQUAL((int ) line[1], composition[toCompIdx(Species::D)]);
		BOOST_REQUIRE_EQUAL((int ) line[2], composition[toCompIdx(Species::T)]);
		BOOST_REQUIRE_EQUAL((int ) line[3], composition[toCompIdx(Species::V)]);
		BOOST_REQUIRE_EQUAL((int ) line[4], composition[toCompIdx(Species::I)]);

		// Check the formation energy
		auto formationEnergy = reactant.getFormationEnergy();
		BOOST_REQUIRE_EQUAL(line[5], formationEnergy);

		// Check the migration energy
		double migrationEnergy = reactant.getMigrationEnergy();
		BOOST_REQUIRE_EQUAL(line[6], migrationEnergy);

		// Check the diffusion factor
		double diffusionFactor = reactant.getDiffusionFactor();
		BOOST_REQUIRE_EQUAL(line[7], diffusionFactor);
	}

	// If the HDF5 file contains initial concentrations
	int tempTimeStep = -2;
	if (HDF5Utils::hasConcentrationGroup("test.h5", tempTimeStep)) {
		// Read the concentrations at the given grid point
		auto returnedVector = HDF5Utils::readGridPoint("test.h5", tempTimeStep,
				gridPoint);

		// Check the size of the vector
		BOOST_REQUIRE_EQUAL(returnedVector.size(), length);
		// Check the values
		for (unsigned int i = 0; i < returnedVector.size(); i++) {
			BOOST_REQUIRE_CLOSE(returnedVector.at(i).at(0),
					concArray[i][0], 0.0001);
			BOOST_REQUIRE_CLOSE(returnedVector.at(i).at(1),
					concArray[i][1], 0.0001);
		}
	}
}

/**
 * Method checking the writing and reading of the surface position specifically
 * in the case of a 2D grid.
 */
BOOST_AUTO_TEST_CASE(checkSurface2D) {
	// Initialize the HDF5 file
	HDF5Utils::initializeFile("test.h5");

	// Set the number of grid points and step size
	int nGrid = 5;
	double stepSize = 0.5;
	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;
	// Write the header in the HDF5 file
	HDF5Utils::fillHeader(nGrid, stepSize);

	// Finalize the HDF5 file
	HDF5Utils::finalizeFile();

	// Open it again to add the concentrations
	HDF5Utils::openFile("test.h5");

	// Set the time step number
	int timeStep = 0;

	// Add the concentration sub group
	HDF5Utils::addConcentrationSubGroup(timeStep, currentTime, previousTime,
			currentTimeStep);

	// Set the surface information in 2D
	std::vector<int> iSurface = { 2, 3, 2, 0, 5 };
	std::vector<double> nInter = { 0.0, 0.0, 0.5, 0.6, 0.5 };
	std::vector<double> previousFlux = { 0.0, 0.1, 3.0, -1.0, 5.0 };

	// Write the surface position
	HDF5Utils::writeSurface2D(timeStep, iSurface, nInter, previousFlux);

	// Close the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	// Read the surface position
	auto surfacePos = HDF5Utils::readSurface2D("test.h5", timeStep);
	// Check all the values
	for (int i = 0; i < surfacePos.size(); i++) {
		BOOST_REQUIRE_EQUAL(surfacePos[i], iSurface[i]);
	}

	// Read the interstitial quantity
	auto nInterstitial = HDF5Utils::readNInterstitial2D("test.h5", timeStep);
	// Check all the values
	for (int i = 0; i < nInterstitial.size(); i++) {
		BOOST_REQUIRE_CLOSE(nInterstitial[i], nInter[i], 0.0001);
	}

	// Read the interstitial flux
	auto previousIFlux = HDF5Utils::readPreviousIFlux2D("test.h5", timeStep);
	// Check all the values
	for (int i = 0; i < previousIFlux.size(); i++) {
		BOOST_REQUIRE_CLOSE(previousIFlux[i], previousFlux[i], 0.0001);
	}
}

/**
 * Method checking the writing and reading of the surface position specifically
 * in the case of a 3D grid.
 */
BOOST_AUTO_TEST_CASE(checkSurface3D) {
	// Initialize the HDF5 file
	HDF5Utils::initializeFile("test.h5");

	// Set the number of grid points and step size
	int nGrid = 5;
	double stepSize = 0.5;
	// Set the time information
	double currentTime = 0.0001;
	double previousTime = 0.00001;
	double currentTimeStep = 0.000001;
	// Write the header in the HDF5 file
	HDF5Utils::fillHeader(nGrid, stepSize);

	// Finalize the HDF5 file
	HDF5Utils::finalizeFile();

	// Open it again to add the concentrations
	HDF5Utils::openFile("test.h5");

	// Set the time step number
	int timeStep = 0;

	// Add the concentration sub group
	HDF5Utils::addConcentrationSubGroup(timeStep, currentTime, previousTime,
			currentTimeStep);

	// Set the surface information in 2D
	std::vector<std::vector<int> > iSurface = { { 2, 4, 1, 0, 5 }, { 2, 3, 2, 0,
			5 }, { 6, 1, 2, 3, 2 } };
	std::vector<std::vector<double> > nInter = { { 0.0, 0.0, 0.0, 0.0, 0.0 }, {
			2.0, 3.0, 2.0, 0.0, 0.5 }, { 0.0, 0.0, 0.0, 0.0, 0.0 } };
	std::vector<std::vector<double> > previousFlux = {
			{ 0.0, 0.0, 0.0, 0.0, 0.0 }, { -2.0, 3.0, 2.0, 0.0, -0.5 }, { 0.0,
					0.0, 0.0, 0.0, 0.0 } };

	// Write the surface position
	HDF5Utils::writeSurface3D(timeStep, iSurface, nInter, previousFlux);

	// Close the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	// Read the surface position
	auto surfacePos = HDF5Utils::readSurface3D("test.h5", timeStep);
	// Check all the values
	for (int i = 0; i < surfacePos.size(); i++) {
		for (int j = 0; j < surfacePos[0].size(); j++) {
			BOOST_REQUIRE_EQUAL(surfacePos[i][j], iSurface[i][j]);
		}
	}
	auto nInterstitial = HDF5Utils::readNInterstitial3D("test.h5", timeStep);
	// Check all the values
	for (int i = 0; i < nInterstitial.size(); i++) {
		for (int j = 0; j < nInterstitial[0].size(); j++) {
			BOOST_REQUIRE_CLOSE(nInterstitial[i][j], nInter[i][j], 0.0001);
		}
	}
	auto previousIFlux = HDF5Utils::readPreviousIFlux3D("test.h5", timeStep);
	// Check all the values
	for (int i = 0; i < previousIFlux.size(); i++) {
		for (int j = 0; j < previousIFlux[0].size(); j++) {
			BOOST_REQUIRE_CLOSE(previousIFlux[i][j], previousFlux[i][j],
					0.0001);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
